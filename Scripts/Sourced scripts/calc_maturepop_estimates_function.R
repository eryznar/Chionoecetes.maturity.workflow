# PURPOSE: function to generate ogives, SAM, and mature/immature bio/abund from sdmTMB maturity models
# Author: Emily Ryznar

calc_maturepop_estimates <- function(model, crab_data, years, species, region,
                                     district, size_1mm, size_min, size_max,
                                     fill_missing_years, output, n_sim) {
  
  outputs <- list()
  
  #-----------------------------
  # Prep specimen data -> sub1
  #-----------------------------
  sub1 <- crab_data$specimen %>%
          filter(
            REGION == region,
            SPECIES == species,
            SHELL_CONDITION == 2,
            SEX == 1,
            !(YEAR == 2025 & SPECIES == "SNOW" & SIZE == 172.5)) %>%
          mutate(
            BIN_5MM = cut_width(SIZE_1MM, width = 5, center = 2.5,
                                closed = "left", dig.lab = 4),
            BIN2 = BIN_5MM) %>%
          separate(BIN2, sep = ",", into = c("LOWER", "UPPER")) %>%
          mutate(
            LOWER = as.numeric(sub('.', '', LOWER)),
            UPPER = as.numeric(gsub('.$', '', UPPER)),
            SIZE_5MM = (UPPER + LOWER) / 2) %>%
          dplyr::select(!c(BIN_5MM, LOWER, UPPER)) %>%
          mutate(
            YEAR_F = as.factor(YEAR),
            YEAR_SCALED = scale(YEAR)) %>%
          st_as_sf(coords = c("LONGITUDE", "LATITUDE"),
                   crs = "+proj=longlat +datum=WGS84") %>%
          st_transform(crs = "+proj=utm +zone=2") %>%
          cbind(st_coordinates(.)) %>%
          as.data.frame() %>%
          mutate(
            LATITUDE = Y / 1000,
            LONGITUDE = X / 1000) 
  
  #-----------------------------
  # District logic
  #-----------------------------
  if (species == "TANNER" && is.null(district)) {
    # TANNER + NULL: keep E/W and add ALL copy
    sub1 <- dplyr::bind_rows(
      sub1,
      sub1 %>% dplyr::mutate(DISTRICT = "ALL")
    )
  } else if (species == "TANNER" && !is.null(district) && district == "ALL") {
    # TANNER + ALL: combine E166 and W166 into ALL only
    sub1 <- sub1 %>%
      dplyr::filter(DISTRICT %in% c("E166", "W166")) %>%
      dplyr::mutate(DISTRICT = "ALL")
  } else if (species == "SNOW" && is.null(district)) {
    # SNOW + NULL: leave all districts
    # sub1 unchanged
  } else if (!is.null(district)) {
    # any species with a specified district: filter to that / those
    sub1 <- sub1 %>%
      dplyr::filter(DISTRICT %in% district)
  }
  
  
  #-----------------------------
  # Do you want the model to fill in missing chela years? 
  #-----------------------------
  sub1 <- if (fill_missing_years) {
    sub1 # Yes
  } else {
    dplyr::filter(sub1, YEAR %in% unique(model$data$YEAR)) # No
  }
 
  #-----------------------------
  # Specify outputs
  #-----------------------------
  all_opts <- c("ogives", "SAM", "cpue", "bioabund")
  if (missing(output) || is.null(output)) {
    output <- all_opts
  } else {
    output <- match.arg(output, choices = all_opts, several.ok = TRUE)
  }
  
  #------------------------------
  # Specify missing chela years
  #------------------------------
  missing_chela <- data.frame(SPECIES  = c("SNOW", "SNOW", "SNOW", "SNOW",
                                           "TANNER", "TANNER"),
                              YEAR = c(2008, 2012, 2014, 2016,
                                       2013, 2015))
  
  #------------------------------------------------------------------
  # sdmTMB simulations
  #------------------------------------------------------------------
  set.seed(1)
  message("Simulating from sdmTMB model")
  pmat.sim <- predict(model, sub1, type = "response", nsim = n_sim)
  
  #------------------------------------------------------------------
  # OGIVES / SAM with SF-at-size uncertainty
  #------------------------------------------------------------------
  if (any(c("ogives", "SAM") %in% output)) {
    
    # 1. Size-bin SF summaries from original design weights
    SF_var <- sub1 %>%
      group_by(YEAR, REGION, DISTRICT, SIZE_5MM, STATION_ID) %>%
      summarise(
        SF_haul = sum(SAMPLING_FACTOR, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      group_by(YEAR, REGION, DISTRICT, SIZE_5MM) %>%
      summarise(
        SF_sum   = sum(SF_haul, na.rm = TRUE),   # total SF at size
        SF_sd    = sd(SF_haul, na.rm = TRUE),    # SD across hauls
        n_hauls  = dplyr::n(),
        .groups  = "drop"
      ) %>%
      mutate(
        SF_sum  = tidyr::replace_na(SF_sum, 0),
        SF_sd   = tidyr::replace_na(SF_sd,  0),
        n_hauls = ifelse(is.na(n_hauls) | n_hauls <= 0, 1, n_hauls),
        SF_se   = SF_sd / sqrt(n_hauls)         # approx SE for SF_sum
      )
    
    # join SF summaries to each crab
    sub1_sf <- sub1 %>%
      left_join(
        SF_var %>% dplyr::select(YEAR, REGION, DISTRICT, SIZE_5MM, SF_sum, SF_se),
        by = c("YEAR", "REGION", "DISTRICT", "SIZE_5MM")
      )
    
    B    <- 25
    nsim <- ncol(pmat.sim)
    
    # SAM helper
    get_sam <- function(size, p) {
      o <- order(size)
      size <- size[o]; p <- p[o]
      if (all(p < 0.5) || all(p > 0.5)) return(NA_real_)
      i_upper <- which(p >= 0.5)[1]
      i_lower <- i_upper - 1
      if (is.na(i_upper) || i_upper <= 1) return(NA_real_)
      size_lower <- size[i_lower]; size_upper <- size[i_upper]
      prop_lower <- p[i_lower];    prop_upper <- p[i_upper]
      size_lower + ((0.5 - prop_lower) / (prop_upper - prop_lower)) *
        (size_upper - size_lower)
    }
    
    ogive_draws_list <- vector("list", nsim)
    SAM_draws_list   <- vector("list", nsim)
    
    # outer loop over sdmTMB sims
    for (s in seq_len(nsim)) {
      message(paste("SF draws for sdmTMB sim", s, "of", nsim))
      
      pmat_s <- pmat.sim[, s]
      sub1_p <- cbind(sub1_sf, pmat = pmat_s)
      
      ogive_inner <- vector("list", B)
      SAM_inner   <- vector("list", B)
      
      # inner loop: SF-at-size draws
      for (b in seq_len(B)) {
        
        # draw SF total for each size bin
        SF_draw_bins <- SF_var %>%
          mutate(
            SF_total_draw = pmax(
              rnorm(n(), mean = SF_sum, sd = SF_se),
              0
            )
          ) %>%
          dplyr::select(YEAR, REGION, DISTRICT, SIZE_5MM, SF_total_draw)
        
        # allocate SF_total_draw back to crabs using original SAMPLING_FACTOR
        sub_draw <- sub1_p %>%
          left_join(
            SF_draw_bins,
            by = c("YEAR", "REGION", "DISTRICT", "SIZE_5MM")
          ) %>%
          group_by(YEAR, REGION, DISTRICT, SIZE_5MM) %>%
          mutate(
            w_rel = ifelse(
              sum(SAMPLING_FACTOR, na.rm = TRUE) > 0,
              SAMPLING_FACTOR / sum(SAMPLING_FACTOR, na.rm = TRUE),
              0
            ),
            SF_draw = w_rel * SF_total_draw
          ) %>%
          ungroup()
        
        ogive_b <- sub_draw %>%
          group_by(YEAR, SPECIES, REGION, DISTRICT, SIZE_5MM) %>%
          summarise(
            denom = sum(SF_draw, na.rm = TRUE),
            num   = sum(pmat * SF_draw, na.rm = TRUE),
            p_b   = ifelse(denom > 0, num / denom, 0),
            .groups = "drop"
          )
        
        ogive_inner[[b]] <- ogive_b
        
        SAM_b <- ogive_b %>%
          group_by(YEAR, SPECIES, REGION, DISTRICT) %>%
          summarise(
            SAM = get_sam(SIZE_5MM, p_b),
            .groups = "drop"
          ) %>%
          mutate(sim = s, inner_draw = b)
        
        SAM_inner[[b]] <- SAM_b
      }
      
      ogive_draws_list[[s]] <- dplyr::bind_rows(ogive_inner)
      SAM_draws_list[[s]]   <- dplyr::bind_rows(SAM_inner)
    }
    
    ogive_all <- dplyr::bind_rows(ogive_draws_list)
    SAM_all   <- dplyr::bind_rows(SAM_draws_list)
    
    #------------------------
    # Ogive summary
    #------------------------
    if ("ogives" %in% output) {
      message("Summarizing ogives")
      
      ogives <- ogive_all %>%
        group_by(YEAR, SPECIES, REGION, DISTRICT, SIZE_5MM) %>%
        summarise(
          PROP_MATURE = mean(p_b, na.rm = TRUE),
          VAR_TOTAL        = var(p_b,  na.rm = TRUE),
          PROP_MATURE_LOGIT = {
            p_clip <- pmin(pmax(p_b, 1e-6), 1 - 1e-6)
            mean(qlogis(p_clip), na.rm = TRUE)
          },
          PROP_MATURE_LOGIT_SD = {
            p_clip <- pmin(pmax(p_b, 1e-6), 1 - 1e-6)
            sd(qlogis(p_clip), na.rm = TRUE)
          },
          .groups = "drop"
        ) %>%
        mutate(
          PROP_MATURE = ifelse(is.nan(PROP_MATURE), 0, PROP_MATURE),
          VAR_TOTAL        = replace_na(VAR_TOTAL, 0),
          VAR_TOTAL        = pmax(VAR_TOTAL, 0),
          PROP_MATURE_SD   = sqrt(VAR_TOTAL),
          PROP_MATURE = pmin(pmax(PROP_MATURE, 0), 1),
          hi_raw = PROP_MATURE + 1.96 * PROP_MATURE_SD,
          lo_raw = PROP_MATURE - 1.96 * PROP_MATURE_SD,
          PROP_MATURE_CI_hi = pmin(1, hi_raw),
          PROP_MATURE_CI_lo = pmax(0, lo_raw),
          PROP_MATURE_SD = (PROP_MATURE_CI_hi - PROP_MATURE) / 1.96
        ) %>%
        left_join(
          SF_var %>%
            dplyr::select(YEAR, REGION, DISTRICT, SIZE_5MM, SF_sum),
          by = c("YEAR", "REGION", "DISTRICT", "SIZE_5MM")
        ) %>%
        mutate(
          SF          = replace_na(SF_sum, 0),
          NUM_MATURE  = SF * PROP_MATURE,
          NUM_IMMATURE= SF - NUM_MATURE,
          TOTAL_CRAB  = SF,
          # Did the model interpolate missing chela years?
          INTERPOLATED = if_else(
            YEAR %in% missing_chela$YEAR[missing_chela$SPECIES == unique(SPECIES)],
            "Y", "N")) %>%
        filter(!(YEAR == 2025 & SPECIES == "SNOW" & SIZE_5MM == 172.5)) %>%# dropping likely miscoded Tanner until this is corrected in the specimen file
        dplyr::select(YEAR, SPECIES, REGION, DISTRICT, SIZE_5MM, PROP_MATURE, PROP_MATURE_SD, 
                      PROP_MATURE_LOGIT, PROP_MATURE_LOGIT_SD, NUM_MATURE, NUM_IMMATURE, TOTAL_CRAB, INTERPOLATED)
      
      outputs$ogives <- ogives
    }
    
    #------------------------
    # SAM summary
    #------------------------
    if ("SAM" %in% output) {
      message("Summarizing SAM")
      
      SAM <- SAM_all %>%
        group_by(YEAR, SPECIES, REGION, DISTRICT) %>%
        summarise(
          VAR_TOTAL = var(SAM,  na.rm = TRUE),
          SAM_SD    = sqrt(VAR_TOTAL),
          SAM  = mean(SAM, na.rm = TRUE),
          .groups   = "drop"
        ) %>%
        full_join(
          expand.grid(YEAR = seq(min(.$YEAR, na.rm = TRUE),
                                max(.$YEAR, na.rm = TRUE)),
                      DISTRICT = unique(.$DISTRICT)),
        ) %>%
        # Did the model interpolate missing chela years?
        mutate(INTERPOLATED = if_else(
          YEAR %in% missing_chela$YEAR[missing_chela$SPECIES == unique(SPECIES)],
          "Y", "N"))
      
      outputs$SAM <- SAM
    }
  } # end if(ogives/SAM)
  
  #------------------------------------------------------------------
  # CPUE 
  #------------------------------------------------------------------
  if ("cpue" %in% output) {
    nsim <- ncol(pmat.sim)
    cpue.df <- data.frame()
    
    for (ii in 1:nsim) {
      message(paste0("Calculating cpue sim ", ii))
      
      fit.sim <- pmat.sim[, ii]
      
      crab_data$specimen <- cbind(sub1, fit.sim) %>%
        filter(!(SPECIES == "TANNER" & DISTRICT == "ALL")) %>% # adding this in until crabpack can handle EBS-wide for tanner
        rename(PROP_MATURE = fit.sim) %>%
        mutate(
          SAMPLING_FACTOR_MATURE   = SAMPLING_FACTOR * PROP_MATURE,
          SAMPLING_FACTOR_IMMATURE = SAMPLING_FACTOR - SAMPLING_FACTOR_MATURE
        )
      
      crab_data$specimen <- crab_data$specimen %>%
        mutate(SAMPLING_FACTOR = SAMPLING_FACTOR_MATURE)
      
      cpue_sim_mature <- suppressMessages(
        crabpack::calc_cpue(
          crab_data = crab_data, species = species,
          size_min = size_min, size_max = size_max, sex = "male", district = unique(crab_data$specimen$DISTRICT),
          bin_1mm = isTRUE(size_1mm),
          shell_condition = "new_hardshell"
        )
      ) %>% mutate(CATEGORY = "Mature male")
      
      crab_data$specimen <- crab_data$specimen %>%
        mutate(SAMPLING_FACTOR = SAMPLING_FACTOR_IMMATURE)
      
      cpue_sim_immature <- suppressMessages(
        crabpack::calc_cpue(
          crab_data = crab_data, species = species,
          size_min = size_min, size_max = size_max, sex = "male", district = unique(crab_data$specimen$DISTRICT),
          bin_1mm = isTRUE(size_1mm),
          shell_condition = "new_hardshell"
        )
      ) %>% mutate(CATEGORY = "Immature male")
      
      cpue_sim <- rbind(cpue_sim_mature, cpue_sim_immature)
      cpue.df  <- rbind(cpue.df, cpue_sim %>% mutate(sim = ii))
    }
    
    cpue.df2 <- cpue.df %>%
      group_by(SPECIES, YEAR, REGION, STATION_ID, LATITUDE, LONGITUDE,
               SHELL_TEXT, DISTRICT, STRATUM, TOTAL_AREA, CATEGORY) %>%
      summarise(
        NSIM = n(),
        COUNT_MEAN    = mean(COUNT,   na.rm = TRUE),
        CPUE_MEAN     = mean(CPUE,    na.rm = TRUE),
        CPUE_MT_MEAN  = mean(CPUE_MT, na.rm = TRUE),
        CPUE_LBS_MEAN = mean(CPUE_LBS,na.rm = TRUE),
        COUNT_VAR     = var(COUNT,   na.rm = TRUE),
        CPUE_VAR      = var(CPUE,    na.rm = TRUE),
        CPUE_MT_VAR   = var(CPUE_MT, na.rm = TRUE),
        CPUE_LBS_VAR  = var(CPUE_LBS,na.rm = TRUE),
        COUNT_SD      = sqrt(COUNT_VAR),
        CPUE_SD       = sqrt(CPUE_VAR),
        CPUE_MT_SD    = sqrt(CPUE_MT_VAR),
        CPUE_LBS_SD   = sqrt(CPUE_LBS_VAR),
        COUNT_CI      = 1.96 * COUNT_SD,
        CPUE_CI       = 1.96 * CPUE_SD,
        CPUE_MT_CI    = 1.96 * CPUE_MT_SD,
        CPUE_LBS_CI   = 1.96 * CPUE_LBS_SD,
        .groups       = "drop"
      )
    
    mature_cpue <- cpue.df2 %>%
      dplyr::select(!c(CPUE_VAR, CPUE_MT_VAR, CPUE_LBS_VAR,
                       CPUE_SD, CPUE_MT_SD, CPUE_LBS_SD,
                       COUNT_VAR, COUNT_SD, NSIM,
                       CPUE_CI, COUNT_CI, CPUE_MT_CI, CPUE_LBS_CI)) %>%
      rename(
        COUNT     = COUNT_MEAN,
        CPUE      = CPUE_MEAN,
        CPUE_MT   = CPUE_MT_MEAN,
        CPUE_LBS  = CPUE_LBS_MEAN
      ) %>%
      identity() %>%
      full_join(
        expand.grid(YEAR = seq(min(.$YEAR, na.rm = TRUE),
                               max(.$YEAR, na.rm = TRUE)),
                    DISTRICT = unique(.$DISTRICT)),
      ) %>%
      mutate(
        #COUNT_CI      = 1.96 * COUNT_SD,
        # CPUE_CI       = 1.96 * CPUE_SD,
        # CPUE_MT_CI    = 1.96 * CPUE_MT_SD,
        # CPUE_LBS_CI   = 1.96 * CPUE_LBS_SD,
        
        # ensure CI widths are non‑negative
        #COUNT_CI    = pmax(COUNT_CI,    0),
        # CPUE_CI    = pmax(CPUE_CI,    0),
        # CPUE_MT_CI   = pmax(CPUE_MT_CI ,   0),
        # CPUE_LBS_CI  = pmax(CPUE_LBS_CI,  0),
        
        # clamp lower bounds at 0
        # COUNT_lo    = pmax(COUNT   - COUNT_CI,   0),
        # COUNT_hi    = COUNT   + COUNT_CI,
        # 
        # CPUE_lo    = pmax(CPUE   - CPUE_CI,   0),
        # CPUE_hi    = CPUE   + CPUE_CI,
        
        # CPUE_MT_lo   = pmax(CPUE_MT  - CPUE_MT_CI,  0),
        # CPUE_MT_hi   = CPUE_MT  + CPUE_MT_CI,
        # 
        # CPUE_LBS_lo  = pmax(CPUE_LBS - CPUE_LBS_CI, 0),
        # CPUE_LBS_hi  = CPUE_LBS + CPUE_LBS_CI,
        # Did the model interpolate missing chela years?
        INTERPOLATED = if_else(
          YEAR %in% missing_chela$YEAR[missing_chela$SPECIES == unique(SPECIES)],
          "Y", "N"))
    
    outputs$mature_cpue <- mature_cpue
  }
  
  #------------------------------------------------------------------
  # BIOABUND 
  #------------------------------------------------------------------
  if ("bioabund" %in% output) {
    nsim <- ncol(pmat.sim)
    bioabund.df <- data.frame()
    
    for (ii in 1:nsim) {
      message(paste0("Calculating bioabund sim ", ii))
      
      fit.sim <- pmat.sim[, ii]
      
      crab_data$specimen <- cbind(sub1, fit.sim) %>%
        filter(!(SPECIES == "TANNER" & DISTRICT == "ALL")) %>% # adding this in until crabpack can handle EBS-wide for tanner
        rename(PROP_MATURE = fit.sim) %>%
        mutate(
          SAMPLING_FACTOR_MATURE   = SAMPLING_FACTOR * PROP_MATURE,
          SAMPLING_FACTOR_IMMATURE = SAMPLING_FACTOR - SAMPLING_FACTOR_MATURE
        )
      
      crab_data$specimen <- crab_data$specimen %>%
        mutate(SAMPLING_FACTOR = SAMPLING_FACTOR_MATURE)
      
      bioabund_sim_mature <- suppressMessages(
        crabpack::calc_bioabund(
          crab_data = crab_data, species = species,
          size_min = size_min, size_max = size_max,
          sex = "male", shell_condition = "new_hardshell"
        )
      ) %>% mutate(CATEGORY = "Mature male")
      
      crab_data$specimen <- crab_data$specimen %>%
        mutate(SAMPLING_FACTOR = SAMPLING_FACTOR_IMMATURE)
      
      bioabund_sim_immature <- suppressMessages(
        crabpack::calc_bioabund(
          crab_data = crab_data, species = species,
          size_min = size_min, size_max = size_max,
          sex = "male", shell_condition = "new_hardshell"
        )
      ) %>% mutate(CATEGORY = "Immature male")
      
      bioabund_sim <- rbind(bioabund_sim_mature, bioabund_sim_immature)
      bioabund.df  <- rbind(bioabund.df, bioabund_sim %>% mutate(sim = ii))
    }
    
    bioabund.df2 <- bioabund.df %>%
      group_by(YEAR, SPECIES, DISTRICT, CATEGORY) %>%
      summarise(
        NSIM = n(),
        ABUNDANCE_MEAN    = mean(ABUNDANCE,   na.rm = TRUE),
        BIOMASS_MT_MEAN   = mean(BIOMASS_MT,  na.rm = TRUE),
        BIOMASS_LBS_MEAN  = mean(BIOMASS_LBS, na.rm = TRUE),
        VAR_ABUNDANCE_between   = var(ABUNDANCE,   na.rm = TRUE),
        VAR_BIOMASS_MT_between  = var(BIOMASS_MT,  na.rm = TRUE),
        VAR_BIOMASS_LBS_between = var(BIOMASS_LBS, na.rm = TRUE),
        VAR_ABUNDANCE_within    = mean((ABUNDANCE_CI   / 1.96)^2, na.rm = TRUE),
        VAR_BIOMASS_MT_within   = mean((BIOMASS_MT_CI  / 1.96)^2, na.rm = TRUE),
        VAR_BIOMASS_LBS_within  = mean((BIOMASS_LBS_CI / 1.96)^2, na.rm = TRUE),
        ABUNDANCE_VAR    = VAR_ABUNDANCE_within   + VAR_ABUNDANCE_between,
        BIOMASS_MT_VAR   = VAR_BIOMASS_MT_within  + VAR_BIOMASS_MT_between,
        BIOMASS_LBS_VAR  = VAR_BIOMASS_LBS_within + VAR_BIOMASS_LBS_between,
        ABUNDANCE_SD     = sqrt(ABUNDANCE_VAR),
        BIOMASS_MT_SD    = sqrt(BIOMASS_MT_VAR),
        BIOMASS_LBS_SD   = sqrt(BIOMASS_LBS_VAR),
        ABUNDANCE_CI     = 1.96 * ABUNDANCE_SD,
        BIOMASS_MT_CI    = 1.96 * BIOMASS_MT_SD,
        BIOMASS_LBS_CI   = 1.96 * BIOMASS_LBS_SD,
        .groups          = "drop"
      )
    
    mature_bioabund <- bioabund.df2 %>%
      dplyr::select(!c(
        VAR_ABUNDANCE_between, VAR_ABUNDANCE_within,
        VAR_BIOMASS_MT_between, VAR_BIOMASS_MT_within,
        VAR_BIOMASS_LBS_between, VAR_BIOMASS_LBS_within,
        ABUNDANCE_VAR, BIOMASS_MT_VAR, BIOMASS_LBS_VAR,
        ABUNDANCE_SD, BIOMASS_MT_SD, BIOMASS_LBS_SD, NSIM)) %>%
      rename(
        ABUNDANCE   = ABUNDANCE_MEAN,
        BIOMASS_MT  = BIOMASS_MT_MEAN,
        BIOMASS_LBS = BIOMASS_LBS_MEAN
      ) %>%
      right_join(
        expand.grid(
          YEAR     = seq(min(.$YEAR, na.rm = TRUE),
                         max(.$YEAR, na.rm = TRUE)),
          CATEGORY = c("Mature male", "Immature male")
        ),
        by = c("YEAR", "CATEGORY")
      ) %>%
      filter(ABUNDANCE > 0) %>% # crabpack codes missing chela years as zero, filtering them here so they are assigned NA in next step
      full_join(
        expand.grid(YEAR = seq(min(.$YEAR, na.rm = TRUE),
                               max(.$YEAR, na.rm = TRUE)),
                    DISTRICT = unique(.$DISTRICT),
                    CATEGORY = unique(.$CATEGORY)),
      ) %>%
      mutate(
        # ensure CI widths are non‑negative
        ABUNDANCE_CI    = pmax(ABUNDANCE_CI,    0),
        BIOMASS_MT_CI   = pmax(BIOMASS_MT_CI,   0),
        BIOMASS_LBS_CI  = pmax(BIOMASS_LBS_CI,  0),
        
        # clamp lower bounds at 0
        ABUNDANCE_lo    = pmax(ABUNDANCE   - ABUNDANCE_CI,   0),
        ABUNDANCE_hi    = ABUNDANCE   + ABUNDANCE_CI,
        
        BIOMASS_MT_lo   = pmax(BIOMASS_MT  - BIOMASS_MT_CI,  0),
        BIOMASS_MT_hi   = BIOMASS_MT  + BIOMASS_MT_CI,
        
        BIOMASS_LBS_lo  = pmax(BIOMASS_LBS - BIOMASS_LBS_CI, 0),
        BIOMASS_LBS_hi  = BIOMASS_LBS + BIOMASS_LBS_CI,
        # Did the model interpolate missing chela years?
        INTERPOLATED = if_else(
            YEAR %in% missing_chela$YEAR[missing_chela$SPECIES == unique(SPECIES)],
            "Y", "N"))

    outputs$mature_bioabund <- mature_bioabund
  }
  
  return(outputs)
}
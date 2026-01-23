# PURPOSE: function to generate ogives, SAM, and mature/immature bio/abund from sdmTMB maturity models

# Author: Emily Ryznar

# NOTES:


# CREATE FUNCTION ----
# Create function to simulate from sdmTMB model, and calculate ogives, SAM, and mature bioabund with uncertainty
calc_maturepop_estimates <- function(model, crab_data, years, species, output){
  
  # list to store all requested outputs
  outputs <- list()
  gc()
  
  # filter specimen data by year and transform to sdmTMB coordinates
  crab_data$specimen %>%
    filter(YEAR %in% years, SPECIES == species, SHELL_CONDITION == 2, SEX == 1) %>%
    mutate(BIN_5MM = cut_width(SIZE_1MM, width = 5, center = 2.5, closed = "left", dig.lab = 4),
           BIN2 = BIN_5MM) %>%
    separate(BIN2, sep = ",", into = c("LOWER", "UPPER")) %>%
    mutate(LOWER = as.numeric(sub('.', '', LOWER)),
           UPPER = as.numeric(gsub('.$', '', UPPER)),
           SIZE_5MM = (UPPER + LOWER)/2) %>%
    dplyr::select(!c(BIN_5MM, LOWER, UPPER)) %>%
    mutate(YEAR_F = as.factor(YEAR),
           YEAR_SCALED = scale(YEAR)) %>%
    st_as_sf(., coords = c("LONGITUDE", "LATITUDE"), crs = "+proj=longlat +datum=WGS84") %>%
    st_transform(., crs = "+proj=utm +zone=2") %>%
    cbind(st_coordinates(.)) %>%
    as.data.frame(.) %>%
    mutate(LATITUDE = Y/1000, # scale to km so values don't get too large
           LONGITUDE = X/1000) %>%
    filter(YEAR %in% unique(model$data$YEAR)) -> sub1
  
  # Specify if output = "", all outputs are produced
  all_opts <- c("ogives", "SAM", "cpue", "bioabund")
  if (missing(output) || is.null(output)) {
    output <- all_opts
  } else {
    output <- match.arg(output, choices = all_opts, several.ok = TRUE)
  }
  
  
 # SIMULATE FROM SDMTMB JOINT PRECISION MATRIX TO ESTIMATE MODEL UNCERTAINTY ----
  set.seed(1)
  message("Simulating from sdmTMB model")
  pmat.sim <- predict(model, sub1, type = "response", nsim = 300)
  
  gc()
  
 # OGIVES/SAM ----
  if(any(c("ogives", "SAM") %in% output)){
  # Calculate haul_level variability for abundance-at-size (SAMPLING FACTOR) to account for its variance below
  SF_var <- sub1 %>%
    group_by(YEAR, DISTRICT, SIZE_5MM, STATION_ID) %>%
    summarise(
      SF_haul = sum(SAMPLING_FACTOR, na.rm = TRUE),  # haul-level SF at size
    ) %>%
    group_by(YEAR, DISTRICT, SIZE_5MM) %>%
    summarise(
      SF_mean = mean(SF_haul, na.rm = TRUE),
      SF_sd   = sd(SF_haul,   na.rm = TRUE),
      n_hauls = n()
    ) %>%
    mutate(
      SF_mean = replace_na(SF_mean, 0), # Make sure SF_mean/SF_sd don’t become NA for single or empty groupings:
      SF_sd   = replace_na(SF_sd, 0)
    )
  
  # Attach SF_mean / SF_sd back to the row-level data (sub1) and set the working SAMPLING_FACTOR equal to the mean
  sub1_sf <- sub1 %>%
    left_join(
      SF_var,
      by = c("YEAR", "DISTRICT", "SIZE_5MM")
    ) %>%
    mutate(
      # Replace (or define) SAMPLING_FACTOR as the mean haul-level SF for that YEAR × DISTRICT × SIZE_5MM
      SAMPLING_FACTOR = SF_mean
    )
  
  # Specify inner draws for SF perturbations per sdmTMB simulation
  B <- 100  
  
  # Specify number of outer sdmTMB simulations
  nsim <- ncol(pmat.sim)
  
  # function to compute SAM from a discrete ogive
  get_sam <- function(size, p) {
    o <- order(size)
    size <- size[o]; p <- p[o]
    
    # if curve never crosses 0.5, return NA
    if (all(p < 0.5) || all(p > 0.5)) return(NA_real_)
    
    # first index where p >= 0.5
    i_upper <- which(p >= 0.5)[1]
    i_lower <- i_upper - 1
    
    # guard for edge cases
    if (is.na(i_upper) || i_upper <= 1) return(NA_real_)
    
    size_lower <- size[i_lower]; size_upper <- size[i_upper]
    prop_lower <- p[i_lower];    prop_upper <- p[i_upper]
    
    # linear interpolation between the two size bins bracketing p = 0.5
    size_lower + ((0.5 - prop_lower) / (prop_upper - prop_lower)) *
      (size_upper - size_lower)
  }
  
  # Containers for to store results from nsim * B loops (outer loop = nsim = model simulations; inner loop = B = SF draws)
  ogive_draws_list <- vector("list", nsim)
  
  SAM_draws_list   <- vector("list", nsim)
  
  # Loop over nsim sdmTMB predictive simulations (outer loop)
  for (s in seq_len(nsim)) {
    message(paste("SF draws for sim", (1:nsim)[s]))
    
    # maturity probabilities for sim s, aligned with rows of sub1_sf
    pmat_s <- pmat.sim[, s]
    
    #Specify containers to store results from B SF draws (inner loop)
    ogive_inner <- vector("list", B)
    SAM_inner   <- vector("list", B)
    
    # Loop to perturb SF and recompute ogives/SAM B times
    for (b in seq_len(B)) {
      # Draw a perturbed SF for each record (row of sub1_sf)
      # mean = SF_mean (attached as SAMPLING_FACTOR), sd = SF_sd
      sf_draw <- suppressWarnings(rnorm(
        nrow(sub1_sf),
        mean = sub1_sf$SAMPLING_FACTOR,
        sd   = sub1_sf$SF_sd
      ))
      
      # Truncate negative draws to zero
      sf_draw[sf_draw < 0] <- 0
      
      # Compute SF-weighted ogive for SF draw b
      ogive_b <- sub1_sf %>%
        mutate(pmat = pmat_s, SF_draw = sf_draw) %>%
        group_by(YEAR, SPECIES, DISTRICT, SIZE_5MM) %>%
        summarise(
          denom = sum(SF_draw),
          num   = sum(pmat * SF_draw),
          p_b   = ifelse(denom > 0, num / denom, 0),
          .groups = "drop"
        )
      
      # store ogive 
      ogive_inner[[b]] <- ogive_b
      
      # Compute SAM for SF draw b
      SAM_b <- ogive_b %>%
        group_by(YEAR, SPECIES, DISTRICT) %>%
        summarise(
          SAM = get_sam(SIZE_5MM, p_b),
          .groups = "drop"
        ) %>%
        mutate(sim = s, inner_draw = b)
      
      # store SAM 
      SAM_inner[[b]] <- SAM_b
      
    } # close SF perturbation loop
    
    # Bind all inner draws for sdmTMB simulation s
    ogive_draws_list[[s]] <- bind_rows(ogive_inner)
    SAM_draws_list[[s]]   <- bind_rows(SAM_inner)
    
  } # close model simulation loop
  
  # Bind across all sdmTMB simulations
  ogive_all <- bind_rows(ogive_draws_list)
  SAM_all   <- bind_rows(SAM_draws_list)
  
  ## Ogive summary
  if ("ogives" %in% output) {
    message("Summarizing ogives")
    
    ogives <- ogive_all %>%
      group_by(YEAR, SPECIES, DISTRICT, SIZE_5MM) %>%
      summarise(
        PROP_MATURE_mean = mean(p_b, na.rm = TRUE),
        VAR_total        = var(p_b,  na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        # if all draws were NA, set maturity to 0 rather than NaN
        PROP_MATURE_mean = ifelse(is.nan(PROP_MATURE_mean), 0, PROP_MATURE_mean),
        VAR_total        = replace_na(VAR_total, 0),
        PROP_MATURE_sd   = sqrt(VAR_total),
        PROP_MATURE_lo   = PROP_MATURE_mean - 1.96 * PROP_MATURE_sd,
        PROP_MATURE_hi   = PROP_MATURE_mean + 1.96 * PROP_MATURE_sd
      )
    
    outputs$ogives <- ogives
  }
  
  ## SAM summary
  if ("SAM" %in% output) {
    message("Summarizing SAM")
    
    SAM <- SAM_all %>%
      group_by(YEAR, SPECIES, DISTRICT) %>%
      summarise(
        SAM_mean  = mean(SAM, na.rm = TRUE),
        VAR_total = var(SAM,  na.rm = TRUE),
        SAM_sd    = sqrt(VAR_total),
        SAM_lo    = SAM_mean - 1.96 * SAM_sd,
        SAM_hi    = SAM_mean + 1.96 * SAM_sd
      ) %>%
      right_join(data.frame(YEAR = seq(min(.$YEAR), max(.$YEAR))),
                 by = "YEAR")
    
    outputs$SAM <- SAM
  }
} # close "if" statement for SAM and ogives

 # CPUE ----
  if("cpue" %in% output){
  # Calculate mature/immature cpue
    nsim <- ncol(pmat.sim)
    cpue.df <- data.frame()
    
    for(ii in 1:nsim){
      message(paste0("Calculating cpue sim ", ii))
      
      fit.sim <- pmat.sim[,ii]
      
      # replace PROP_MATURE with each model simulation draw
      crab_data$specimen <- cbind(sub1, fit.sim) %>%
        rename(PROP_MATURE = fit.sim) %>%
        mutate(SAMPLING_FACTOR_MATURE = SAMPLING_FACTOR * PROP_MATURE,
               SAMPLING_FACTOR_IMMATURE = SAMPLING_FACTOR - SAMPLING_FACTOR_MATURE)
      
      # calculate cpue for each simulation
      crab_data$specimen <- crab_data$specimen %>%
        mutate(SAMPLING_FACTOR = SAMPLING_FACTOR_MATURE) # specifying sampling factor as mat_sf so crabpack recognizes
      # Mature
      cpue_sim_mature <-  crabpack::calc_cpue(crab_data = crab_data, species = species,
                                                      size_min = NULL, size_max = NULL,  sex = "male",
                                                      shell_condition = "new_hardshell") %>%
        mutate(CATEGORY = "Mature male")
      
      # Immature
      crab_data$specimen <- crab_data$specimen %>%
        mutate(SAMPLING_FACTOR = SAMPLING_FACTOR_IMMATURE) # specifying sampling factor as mat_sf so crabpack recognizes
      
      cpue_sim_immature <-  crabpack::calc_cpue(crab_data = crab_data, species = species,
                                                        size_min = NULL, size_max = NULL,  sex = "male",
                                                        shell_condition = "new_hardshell") %>%
        mutate(CATEGORY = "Immature male")
      
      # Bind
      cpue_sim <- rbind(cpue_sim_mature, cpue_sim_immature)
      
      cpue.df <- rbind(cpue.df,  cpue_sim %>% mutate(sim = ii))
    }
    
    # Now propagate both model and survey uncertainty into biomass/abundance
    cpue.df2 <-
      cpue.df %>%
      group_by(SPECIES, YEAR, REGION, STATION_ID, LATITUDE, LONGITUDE, SHELL_TEXT, DISTRICT, STRATUM, TOTAL_AREA, CATEGORY) %>%
      summarise(
        NSIM = n(),
        
        ## Means across maturity draws
        CPUE_MEAN    = mean(CPUE,   na.rm = TRUE),
        CPUE_MT_MEAN   = mean(CPUE_MT,  na.rm = TRUE),
        CPUE_LBS_MEAN   = mean(CPUE_LBS,  na.rm = TRUE),
        
        ## Between-simulation variance (model uncertainty), no expansion uncertainty from db
        CPUE_VAR   = var(CPUE,   na.rm = TRUE),
        CPUE_MT_VAR  = var(CPUE_MT,  na.rm = TRUE),
        CPUE_LBS_VAR = var(CPUE_LBS, na.rm = TRUE),
        
        ## SDs and 95% CIs
        CPUE_SD     = sqrt(CPUE_VAR),
        CPUE_MT_SD    = sqrt(CPUE_MT_VAR),
        CPUE_LBS_SD   = sqrt(CPUE_LBS_VAR),
        
        CPUE_CI   = 1.96 * CPUE_SD,
        CPUE_MT_CI  = 1.96 * CPUE_MT_SD,
        CPUE_LBS_CI = 1.96 * CPUE_LBS_SD,
      
      )
    
    mature_cpue <-
      cpue.df2 %>%
      dplyr::select(!c(CPUE_VAR,
                       CPUE_MT_VAR,
                       CPUE_LBS_VAR,
                       CPUE_SD, CPUE_MT_SD, CPUE_LBS_SD)) %>%
      rename(CPUE      = CPUE_MEAN,
             CPUE_MT     = CPUE_MT_MEAN,
             CPUE_LBS    = CPUE_LBS_MEAN) %>%
      mutate(
        Estimator     = "sdmTMB",
        # Assign NAs for years with no snow chela data
        CPUE = case_when(YEAR %in% c(2008, 2012, 2014, 2016, 2020) & SPECIES == "SNOW" ~ NA,
                              TRUE ~ CPUE),
        CPUE_CI = case_when(YEAR %in% c(2008, 2012, 2014, 2016, 2020)  & SPECIES == "SNOW" ~ NA,
                                 TRUE ~ CPUE_CI),
        CPUE_MT = case_when(YEAR %in% c(2008, 2012, 2014, 2016, 2020) & SPECIES == "SNOW" ~ NA,
                               TRUE ~ CPUE_MT),
        CPUE_MT_CI = case_when(YEAR %in% c(2008, 2012, 2014, 2016, 2020)  & SPECIES == "SNOW" ~ NA,
                                  TRUE ~ CPUE_MT_CI),
        CPUE_LBS = case_when(YEAR %in% c(2008, 2012, 2014, 2016, 2020)  & SPECIES == "SNOW" ~ NA,
                                TRUE ~ CPUE_LBS),
        CPUE_LBS_CI = case_when(YEAR %in% c(2008, 2012, 2014, 2016, 2020) & SPECIES == "SNOW" ~ NA,
                                   TRUE ~ CPUE_LBS_CI),
        # Assign NAs for years with no Tanner data by district
        CPUE = case_when(YEAR %in% c(2011, 2013, 2015, 2020) & SPECIES == "TANNER" & DISTRICT == "E166" ~ NA,
                              YEAR %in% c(2013, 2015, 2020) & SPECIES == "TANNER" & DISTRICT == "W166" ~ NA,
                              TRUE ~ CPUE),
        CPUE_CI = case_when(YEAR %in% c(2011, 2013, 2015, 2020) & SPECIES == "TANNER" & DISTRICT == "E166" ~ NA,
                                 YEAR %in% c(2013, 2015, 2020) & SPECIES == "TANNER" & DISTRICT == "W166" ~ NA,
                                 TRUE ~ CPUE_CI),
        CPUE_MT = case_when(YEAR %in% c(2011, 2013, 2015, 2020) & SPECIES == "TANNER" & DISTRICT == "E166" ~ NA,
                               YEAR %in% c(2013, 2015, 2020) & SPECIES == "TANNER" & DISTRICT == "W166" ~ NA,
                               TRUE ~ CPUE_MT),
        CPUE_MT_CI = case_when(YEAR %in% c(2011, 2013, 2015, 2020) & SPECIES == "TANNER" & DISTRICT == "E166" ~ NA,
                                  YEAR %in% c(2013, 2015, 2020) & SPECIES == "TANNER" & DISTRICT == "W166" ~ NA,
                                  TRUE ~ CPUE_MT_CI),
        CPUE_LBS = case_when(YEAR %in% c(2011, 2013, 2015, 2020) & SPECIES == "TANNER" & DISTRICT == "E166" ~ NA,
                                YEAR %in% c(2013, 2015, 2020) & SPECIES == "TANNER" & DISTRICT == "W166" ~ NA,
                                TRUE ~ CPUE_LBS),
        CPUE_LBS_CI = case_when(YEAR %in% c(2011, 2013, 2015, 2020) & SPECIES == "TANNER" & DISTRICT == "E166" ~ NA,
                                   YEAR %in% c(2013, 2015, 2020) & SPECIES == "TANNER" & DISTRICT == "W166" ~ NA,
                                   TRUE ~ CPUE_LBS_CI))
    
    outputs$mature_cpue <- mature_cpue
    
    #write.csv(mature_cpue, paste0("./Maturity data processing/Output/", species, "_maturecpue.csv"))
  }
 
 # BIOMASS/ABUNDANCE ----
  if("bioabund" %in% output){
  # Calculate mature/immature biomass/abundance, combining model uncertainty with design-based expansion uncertainty ----
  nsim <- ncol(pmat.sim)
  bioabund.df <- data.frame()
  
  for(ii in 1:nsim){
   message(paste0("Calculating bioabund sim ", ii))
    
    fit.sim <- pmat.sim[,ii]
    
    # replace PROP_MATURE with each model simulation draw
    crab_data$specimen <- cbind(sub1, fit.sim) %>%
      rename(PROP_MATURE = fit.sim) %>%
      mutate(SAMPLING_FACTOR_MATURE = SAMPLING_FACTOR * PROP_MATURE,
             SAMPLING_FACTOR_IMMATURE = SAMPLING_FACTOR - SAMPLING_FACTOR_MATURE)
    
    # calculate bioabund for each simulation
    crab_data$specimen <- crab_data$specimen %>%
      mutate(SAMPLING_FACTOR = SAMPLING_FACTOR_MATURE) # specifying sampling factor as mat_sf so crabpack recognizes
    # Mature
    bioabund_sim_mature <-  crabpack::calc_bioabund(crab_data = crab_data, species = species,
                                                    size_min = NULL, size_max = NULL,  sex = "male",
                                                    shell_condition = "new_hardshell") %>%
      mutate(CATEGORY = "Mature male")
    
    # Immature
    crab_data$specimen <- crab_data$specimen %>%
      mutate(SAMPLING_FACTOR = SAMPLING_FACTOR_IMMATURE) # specifying sampling factor as mat_sf so crabpack recognizes
    
    bioabund_sim_immature <-  crabpack::calc_bioabund(crab_data = crab_data, species = species,
                                                      size_min = NULL, size_max = NULL,  sex = "male",
                                                      shell_condition = "new_hardshell") %>%
      mutate(CATEGORY = "Immature male")
    
    # Bind
    bioabund_sim <- rbind(bioabund_sim_mature, bioabund_sim_immature)
    
    bioabund.df <- rbind(bioabund.df,  bioabund_sim %>% mutate(sim = ii))
  }
  
  # Now propagate both model and survey uncertainty into biomass/abundance
  bioabund.df2 <-
    bioabund.df %>%
    group_by(YEAR, SPECIES, DISTRICT, CATEGORY) %>%
    summarise(
      NSIM = n(),
      
      ## Means across maturity draws
      ABUNDANCE_MEAN    = mean(ABUNDANCE,   na.rm = TRUE),
      BIOMASS_MT_MEAN   = mean(BIOMASS_MT,  na.rm = TRUE),
      BIOMASS_LBS_MEAN  = mean(BIOMASS_LBS, na.rm = TRUE),
      
      ## Between-simulation variance (model uncertainty)
      VAR_ABUNDANCE_between   = var(ABUNDANCE,   na.rm = TRUE),
      VAR_BIOMASS_MT_between  = var(BIOMASS_MT,  na.rm = TRUE),
      VAR_BIOMASS_LBS_between = var(BIOMASS_LBS, na.rm = TRUE),
      
      ## Within-simulation variance (survey design, from CIs)
      VAR_ABUNDANCE_within   = mean( (ABUNDANCE_CI   / 1.96)^2, na.rm = TRUE ),
      VAR_BIOMASS_MT_within  = mean( (BIOMASS_MT_CI  / 1.96)^2, na.rm = TRUE ),
      VAR_BIOMASS_LBS_within = mean( (BIOMASS_LBS_CI / 1.96)^2, na.rm = TRUE ),
      
      ## Total variance = within + between
      ABUNDANCE_VAR    = VAR_ABUNDANCE_within   + VAR_ABUNDANCE_between,
      BIOMASS_MT_VAR   = VAR_BIOMASS_MT_within  + VAR_BIOMASS_MT_between,
      BIOMASS_LBS_VAR  = VAR_BIOMASS_LBS_within + VAR_BIOMASS_LBS_between,
      
      ## SDs and 95% CIs
      ABUNDANCE_SD     = sqrt(ABUNDANCE_VAR),
      BIOMASS_MT_SD    = sqrt(BIOMASS_MT_VAR),
      BIOMASS_LBS_SD   = sqrt(BIOMASS_LBS_VAR),
      
      ABUNDANCE_CI   = 1.96 * ABUNDANCE_SD,
      BIOMASS_MT_CI  = 1.96 * BIOMASS_MT_SD,
      BIOMASS_LBS_CI = 1.96 * BIOMASS_LBS_SD,

    )
  
  mature_bioabund <-
    bioabund.df2 %>%
    dplyr::select(!c(VAR_ABUNDANCE_between, VAR_ABUNDANCE_within,
                     VAR_BIOMASS_MT_between, VAR_BIOMASS_MT_within,
                     VAR_BIOMASS_LBS_between, VAR_BIOMASS_LBS_within,
                     ABUNDANCE_VAR, BIOMASS_MT_VAR, BIOMASS_LBS_VAR,
                     ABUNDANCE_SD, BIOMASS_MT_SD, BIOMASS_LBS_SD)) %>%
    rename(ABUNDANCE      = ABUNDANCE_MEAN,
           BIOMASS_MT     = BIOMASS_MT_MEAN,
           BIOMASS_LBS    = BIOMASS_LBS_MEAN) %>%
    right_join(expand.grid(YEAR     = seq(min(.$YEAR), max(.$YEAR)), CATEGORY = c("Mature male", "Immature male"))) %>%
    mutate(Estimator     = "sdmTMB",
           ABUNDANCE     = ABUNDANCE    / 1e6,
           ABUNDANCE_CI  = ABUNDANCE_CI / 1e6) %>%
    mutate(
      # Assign NAs for years with no snow chela data
      ABUNDANCE = case_when(YEAR %in% c(2008, 2012, 2014, 2016, 2020) & SPECIES == "SNOW" ~ NA,
                            TRUE ~ ABUNDANCE),
      ABUNDANCE_CI = case_when(YEAR %in% c(2008, 2012, 2014, 2016, 2020)  & SPECIES == "SNOW" ~ NA,
                               TRUE ~ ABUNDANCE_CI),
      BIOMASS_MT = case_when(YEAR %in% c(2008, 2012, 2014, 2016, 2020) & SPECIES == "SNOW" ~ NA,
                             TRUE ~ BIOMASS_MT),
      BIOMASS_MT_CI = case_when(YEAR %in% c(2008, 2012, 2014, 2016, 2020)  & SPECIES == "SNOW" ~ NA,
                                TRUE ~ BIOMASS_MT_CI),
      BIOMASS_LBS = case_when(YEAR %in% c(2008, 2012, 2014, 2016, 2020)  & SPECIES == "SNOW" ~ NA,
                              TRUE ~ BIOMASS_LBS),
      BIOMASS_LBS_CI = case_when(YEAR %in% c(2008, 2012, 2014, 2016, 2020) & SPECIES == "SNOW" ~ NA,
                                 TRUE ~ BIOMASS_LBS_CI),
      # Assign NAs for years with no Tanner data by district
      ABUNDANCE = case_when(YEAR %in% c(2011, 2013, 2015, 2020) & SPECIES == "TANNER" & DISTRICT == "E166" ~ NA,
                            YEAR %in% c(2013, 2015, 2020) & SPECIES == "TANNER" & DISTRICT == "W166" ~ NA,
                            TRUE ~ ABUNDANCE),
      ABUNDANCE_CI = case_when(YEAR %in% c(2011, 2013, 2015, 2020) & SPECIES == "TANNER" & DISTRICT == "E166" ~ NA,
                               YEAR %in% c(2013, 2015, 2020) & SPECIES == "TANNER" & DISTRICT == "W166" ~ NA,
                               TRUE ~ ABUNDANCE_CI),
      BIOMASS_MT = case_when(YEAR %in% c(2011, 2013, 2015, 2020) & SPECIES == "TANNER" & DISTRICT == "E166" ~ NA,
                             YEAR %in% c(2013, 2015, 2020) & SPECIES == "TANNER" & DISTRICT == "W166" ~ NA,
                             TRUE ~ BIOMASS_MT),
      BIOMASS_MT_CI = case_when(YEAR %in% c(2011, 2013, 2015, 2020) & SPECIES == "TANNER" & DISTRICT == "E166" ~ NA,
                                YEAR %in% c(2013, 2015, 2020) & SPECIES == "TANNER" & DISTRICT == "W166" ~ NA,
                                TRUE ~ BIOMASS_MT_CI),
      BIOMASS_LBS = case_when(YEAR %in% c(2011, 2013, 2015, 2020) & SPECIES == "TANNER" & DISTRICT == "E166" ~ NA,
                              YEAR %in% c(2013, 2015, 2020) & SPECIES == "TANNER" & DISTRICT == "W166" ~ NA,
                              TRUE ~ BIOMASS_LBS),
      BIOMASS_LBS_CI = case_when(YEAR %in% c(2011, 2013, 2015, 2020) & SPECIES == "TANNER" & DISTRICT == "E166" ~ NA,
                                 YEAR %in% c(2013, 2015, 2020) & SPECIES == "TANNER" & DISTRICT == "W166" ~ NA,
                                 TRUE ~ BIOMASS_LBS_CI))
  
  outputs$mature_bioabund <- mature_bioabund
  #write.csv(mature_bioabund, paste0("./Maturity data processing/Output/", species, "_maturebioabund.csv"))
  }
  
 # RETURN OUTPUTS ----
  return(outputs)
  
}



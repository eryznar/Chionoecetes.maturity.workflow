# PURPOSE: to generate ogives, SAM, and mature bio/abund for sdmTMB maturity models and compare to legacy workflow

# Author: Emily Ryznar

# NOTES:


# LOAD LIBS/PARAMS ---------------------------------------------------------------------------------------
source("./Maturity data processing/Scripts/load_libs_params.R")

# CREATE FUNCTION ----
# Create function to simulate from sdmTMB model, and calculate ogives, SAM, and mature bioabund with uncertainty
calc_maturepop_estimates <- function(model, crab_data, years, species){
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
  
  
  # 1) Simulate sdmTMB model predictive uncertainty in pmat from joint-precision matrix ----
  set.seed(1)
  print("Simulating from sdmTMB model")
  pmat.sim <- predict(model, sub1, type = "response", nsim = 500)
  
  
  # 2) Calculate weighted mean pmat by simulation, size bin, district, year ----
  wmean_pmat <- cbind(sub1, pmat.sim) %>%
    group_by(YEAR, SPECIES, DISTRICT, SIZE_5MM) %>%
    summarise(
      across(
        matches("^[0-9]+$"),
        ~ sum(.x * SAMPLING_FACTOR) / sum(SAMPLING_FACTOR),
        .names = "wmean_{.col}"
      ),
      .groups = "drop"
    ) %>%
    pivot_longer(
      starts_with("wmean_"),
      names_to  = "sim",
      values_to = "wmean"
    ) %>%
    mutate(
      sim = as.integer(gsub("wmean_", "", sim))
    )
  
  # 3) Calculate yearly ogives with uncertainty ----
  print("Calculating ogives")
  ogives <- wmean_pmat %>%
    group_by(YEAR, SPECIES, DISTRICT, SIZE_5MM) %>%
    summarise(
      PROP_MATURE_mean = mean(wmean),
      PROP_MATURE_lo   = quantile(wmean, 0.025),
      PROP_MATURE_hi   = quantile(wmean, 0.975),
      .groups = "drop"
    ) 
  
  write.csv(ogives, paste0("./Maturity data processing/Output/", species, "_ogives.csv"))
  
  
  # 4) Calculate yearly SAM with uncertainty ----
  # Specify SAM function
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
  
  # Calculate SAM
  print("Calculating SAM")
    SAM <- wmean_pmat %>%
      group_by(YEAR, SPECIES, DISTRICT, sim) %>% # first get SAM by simulation
      summarise(
        SAM = get_sam(SIZE_5MM, wmean),
        .groups = "drop"
      ) %>%
      group_by(YEAR, SPECIES, DISTRICT) %>% # then summarise across simulations
      reframe(
        SAM_mean = mean(SAM, na.rm = TRUE),
        SAM_sd   = sd(SAM,   na.rm = TRUE),
        SAM_lo   = quantile(SAM, 0.025, na.rm = TRUE),
        SAM_hi   = quantile(SAM, 0.975, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      right_join(., data.frame(YEAR = seq(min(.$YEAR), max(.$YEAR))))
    
    write.csv(SAM, paste0("./Maturity data processing/Output/", species, "_SAM.csv"))
    
  # 5) Calculate mature/immature biomass/abundance, combining model uncertainty with design-based expansion uncertainty ----
    nsim <- ncol(pmat.sim)
    bioabund.df <- data.frame()

    for(ii in 1:nsim){
      print(paste0("Calculating bioabund sim ", ii))

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
        
        .groups = "drop"
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

    write.csv(mature_bioabund, paste0("./Maturity data processing/Output/", species, "_maturebioabund.csv"))
    
 return(list(ogives = ogives, SAM = SAM, mature_bioabund = mature_bioabund))
    #return(list(ogives = ogives, SAM = SAM))
    
}

# SNOW CRAB ----
  # Specify function parameters ----
  snow_mod <- readRDS("./Maturity data processing/Doc/Snow models/sdmTMB_spVAR_SIZE_k300.rda")
  snow_dat <- readRDS("./Maturity data processing/Data/snow_survey_specimenEBS.rda")
  snow_yrs <- c(1989:2019, 2021:2025)
  species <- "SNOW"
  
  # Run function ----
  calc_maturepop_estimates(snow_mod, snow_dat, snow_yrs, species) -> snow.out
  
  # OGIVES ----
  # sdmTMB
  snow_ogive_sdmTMB <- snow.out$ogives %>%
    mutate(Estimator = "sdmTMB")
  
  # legacy using gam prop_mature interpolation of 10mm down to 5mm bins
  snow_legacy_ogives <- read.csv("./Maturity data processing/Doc/snow_pmat_5mminterp.csv") %>%
    mutate(SPECIES = "SNOW", DISTRICT = "ALL") %>%
    dplyr::select(!X) %>%
    rename(SIZE = SIZE_BIN) %>%
    mutate(Estimator = "Legacy")
  
  # Combine and plot
  ogive_plot_dat <- rbind(snow_ogive_sdmTMB %>% rename(PROP_MATURE = PROP_MATURE_mean, SIZE = SIZE_5MM), snow_legacy_ogives %>% mutate(PROP_MATURE_lo = NA, PROP_MATURE_hi = NA))
  
  ggplot(ogive_plot_dat,  aes(SIZE, PROP_MATURE, color = Estimator))+
    geom_ribbon(aes(SIZE, ymin = PROP_MATURE_lo, ymax = PROP_MATURE_hi, fill = Estimator), color = NA, alpha = 0.35)+
    geom_line(linewidth = 1)+
    scale_color_manual(values = c(
      "Legacy"       = "darkgoldenrod",
      "sdmTMB"       = "cadetblue",name = ""))+
    scale_fill_manual(values = c(
      "Legacy"       = "darkgoldenrod",
      "sdmTMB"       = "cadetblue"), name = "", guide= "none")+
    facet_wrap(~YEAR)+
    geom_rug(ogive_plot_dat %>% filter(Estimator == "sdmTMB"), mapping =aes(x = SIZE), sides = "b", inherit.aes = TRUE) +
    geom_hline(yintercept = 0.5, linetype = "dashed")+
    theme_bw()+
    xlim(0, 150)+
    ggtitle("Snow")+
    #scale_x_continuous(breaks = seq(0, 175, by = 10))+
    ylab("Proportion mature")+
    xlab("Carapace width (mm)")+
    theme(legend.position = "bottom", legend.direction = "horizontal",
          legend.text = element_text(size = 12),
          legend.title = element_blank(),
          axis.title = element_text(size = 12),
          strip.text = element_text(size = 12),
          axis.text = element_text(size = 10))
  
  ggsave("./Maturity data processing/Doc/snow_ogive_comparison_legacy.sdmTMBwitherror.png", width  = 10, height = 10, units = "in")
  
  ggplot(snow_legacy_ogives %>% filter(!YEAR %in% c(2008, 2012, 2014, 2016)), aes(SIZE, PROP_MATURE, color = Estimator, group = YEAR))+
    geom_line(linewidth = 1)+
    scale_color_manual(values = c(
      "Legacy"       = "darkgoldenrod"))+
    #facet_wrap(~YEAR)+
    geom_hline(yintercept = 0.5, linetype = "dashed")+
    theme_bw()+
    ggtitle("Snow")+
    #scale_x_continuous(breaks = seq(0, 175, by = 10))+
    ylab("Probability mature")+
    xlab("Carapace width (mm)")+
    theme(legend.position = "bottom", legend.direction = "horizontal",
          legend.text = element_text(size = 12),
          legend.title = element_blank(),
          axis.title = element_text(size = 12),
          strip.text = element_text(size = 12),
          axis.text = element_text(size = 10)) -> l.1
  # SAM ----
  # sdmTMB SAM
  snow_SAM <- snow.out$SAM
  
  # Legacy SAM
  legacy_SAM <- read.csv("./Maturity data processing/Data/snow_legacy_ogives.csv") %>% 
    dplyr::select(YEAR, B_EST, B_SE) %>% 
    distinct() %>%
    rename(SAM = B_EST) %>%
    right_join(., data.frame(YEAR = seq(min(.$YEAR), max(.$YEAR)))) 
  
  # Plot
  ggplot()+
    geom_line(legacy_SAM, mapping = aes(YEAR, SAM, color = as.factor(1)), linewidth = 1)+
    geom_point(legacy_SAM, mapping = aes(YEAR, SAM, color = as.factor(1)), size = 2)+
    geom_errorbar(legacy_SAM, mapping = aes(YEAR, ymin = SAM-B_SE*1.96, ymax = SAM + B_SE*1.96, color = as.factor(1)), linewidth = 1)+
    geom_line(snow_SAM, mapping = aes(YEAR, SAM_mean, color = as.factor(2)), linewidth = 1)+
    geom_point(snow_SAM, mapping = aes(YEAR, SAM_mean, color = as.factor(2)), size = 2)+
    geom_errorbar(snow_SAM, mapping = aes(YEAR, ymin = SAM_lo, ymax = SAM_hi, color = as.factor(2)), linewidth = 1)+
    scale_color_manual(values = c("darkgoldenrod", "cadetblue"), labels = c("Legacy", "sdmTMB"), name = "")+
    scale_fill_manual(values = c("darkgoldenrod", "cadetblue"), labels = c("Legacy", "sdmTMB"), name = "")+
    theme_bw()+
    ggtitle("Snow")+
    xlab("Year")+
    ylab("Size at 50% maturity (mm)")+
    theme(legend.position = "bottom", legend.direction = "horizontal",
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          legend.text = element_text(size = 12)) -> ls.1
  
  ggsave("./Maturity data processing/Doc/snow_SAM_comparison.png", width = 8, height = 6, units = "in")
  
  # Plot
  ggplot()+
    geom_line(legacy_SAM, mapping = aes(YEAR, SAM, color = as.factor(1)), linewidth = 1)+
    geom_point(legacy_SAM, mapping = aes(YEAR, SAM, color = as.factor(1)), size = 2)+
    geom_errorbar(legacy_SAM, mapping = aes(YEAR, ymin = SAM-B_SE*1.96, ymax = SAM + B_SE*1.96, color = as.factor(1)), linewidth = 1)+
    scale_color_manual(values = c("darkgoldenrod"), labels = c("Legacy"), name = "")+
    scale_fill_manual(values = c("darkgoldenrod"), labels = c("Legacy"), name = "")+
    theme_bw()+
    ggtitle("Snow")+
    xlab("Year")+
    ylab("Size at 50% maturity (mm)")+
    theme(legend.position = "bottom", legend.direction = "horizontal",
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          legend.text = element_text(size = 12)) -> ls.1
  
  
  # MATURE BIOMASS/ABUNDANCE ----
  # sdmTMB
  snow_sdmTMB_bioabund <- snow.out$mature_bioabund %>% 
                    dplyr::select(!NSIM)
  
  ggplot(snow_sdmTMB_bioabund %>% filter(YEAR >= 1990), aes(YEAR, ABUNDANCE))+
    geom_point()+
    geom_line()+
    facet_wrap(~CATEGORY)
  
  # Legacy
  snow_legacy_spec <- snow_dat
  snow_legacy_spec$specimen <- snow_dat$specimen %>%
    filter(YEAR %in% snow_yrs, SPECIES == "SNOW", SHELL_CONDITION == 2, SEX == 1) %>%
    mutate(BIN_5MM = cut_width(SIZE_1MM, width = 5, center = 2.5, closed = "left", dig.lab = 4),
           BIN2 = BIN_5MM) %>%
    separate(BIN2, sep = ",", into = c("LOWER", "UPPER")) %>%
    mutate(LOWER = as.numeric(sub('.', '', LOWER)),
           UPPER = as.numeric(gsub('.$', '', UPPER)),
           SIZE_5MM = (UPPER + LOWER)/2) %>%
    dplyr::select(!c(BIN_5MM, LOWER, UPPER)) 
  
  snow_pars <- read.csv("./Maturity data processing/Data/snow_legacy_modelpars.csv") 
  
  snow_legacy_spec$specimen <-  snow_legacy_spec$specimen %>% 
    right_join(., snow_pars) %>%
    mutate(PROP_MATURE = (1/(1 + exp(-A_EST * (SIZE_5MM - B_EST))))) %>%
    mutate(SAMPLING_FACTOR = SAMPLING_FACTOR * PROP_MATURE) %>%
    dplyr::select(!X)
  
  # Use crabpack to calculate bioabund
  snow_legacy_bioabund <-  crabpack::calc_bioabund(crab_data = snow_legacy_spec, species = "SNOW", 
                                                   size_min = NULL, size_max = NULL,  sex = "male", 
                                                   shell_condition = "new_hardshell") %>%
    right_join(., data.frame(YEAR = seq(min(.$YEAR), max(.$YEAR)))) %>%
    mutate(Estimator = "Legacy",
           ABUNDANCE = ABUNDANCE/1e6,
           ABUNDANCE_CI = ABUNDANCE_CI/1e6)
 
  # Bind and specify missing years
  snow_bioabund_dat <- rbind(snow_legacy_bioabund%>%
                               mutate(CATEGORY = "Mature male") %>%
                               dplyr::select(names(.)[names(.) %in% colnames(snow_sdmTMB_bioabund)]), 
                             snow_sdmTMB_bioabund %>% filter(CATEGORY == "Mature male"))
  
  # Plot
  p1 <- ggplot()+
    geom_point(snow_bioabund_dat %>% filter(YEAR %in% c(2013, 2015)), mapping = aes(YEAR, ABUNDANCE, color = Estimator))+
    geom_errorbar(snow_bioabund_dat %>% filter(YEAR %in% c(2013, 2015)), mapping = aes(x = YEAR, ymin = ABUNDANCE - ABUNDANCE_CI, ymax = ABUNDANCE + ABUNDANCE_CI, 
                                                                                       color= Estimator), linewidth = 1)+
    geom_ribbon(snow_bioabund_dat, mapping = aes(x = YEAR, ymin = ABUNDANCE - ABUNDANCE_CI, ymax = ABUNDANCE + ABUNDANCE_CI, 
                                                 fill= Estimator), alpha = 0.15)+
    geom_line(snow_bioabund_dat, mapping = aes(YEAR, ABUNDANCE, color = Estimator), linewidth = 0.75)+
    scale_color_manual(values = c("darkgoldenrod", "cadetblue"), name = "")+
    scale_fill_manual(values = c("darkgoldenrod", "cadetblue"), name = "")+
    theme_bw()+
    ylab("Abundance (millions)")+
    ggtitle("Snow")+
    xlab("Year")+
    theme(
          axis.text.y = element_text(size = 12), 
          axis.title.y = element_text(size = 12),
          legend.position = "bottom",
          legend.direction = "horizontal")
  
  p2 <- ggplot()+
    geom_point(snow_bioabund_dat %>% filter(YEAR %in% c(2013, 2015)), mapping = aes(YEAR, BIOMASS_MT, color = Estimator))+
    geom_errorbar(snow_bioabund_dat %>% filter(YEAR %in% c(2013, 2015)), mapping = aes(x = YEAR, ymin = BIOMASS_MT - BIOMASS_MT_CI, ymax = BIOMASS_MT + BIOMASS_MT_CI, 
                                                                                       color= Estimator) , linewidth = 1)+
    geom_ribbon(snow_bioabund_dat, mapping = aes(x = YEAR, ymin = BIOMASS_MT - BIOMASS_MT_CI, ymax = BIOMASS_MT + BIOMASS_MT_CI, 
                                                 fill= Estimator), alpha = 0.15)+
    geom_line(snow_bioabund_dat, mapping = aes(YEAR, BIOMASS_MT, color = Estimator), linewidth = 0.75)+
    scale_color_manual(values = c("darkgoldenrod", "cadetblue"), name = "")+
    scale_fill_manual(values = c("darkgoldenrod", "cadetblue"), name = "")+
    theme_bw()+
    ylab("Biomass (tons)")+
    #ggtitle("Morphometric mature male biomass (SH2)")+
    xlab("Year")+
    theme(axis.text = element_text(size = 12),
          legend.position = "bottom",
          legend.direction = "horizontal",
          axis.title = element_text(size = 12))
  
  p1/p2 + plot_layout(guides = "collect") & 
    theme(legend.position = "bottom", legend.text = element_text(size = 12))
  
  ggsave("./Maturity data processing/Doc/snow_bioabund_comparison_legacy.sdmTMB.png", width = 7, height = 6, units = "in")
  
# TANNER CRAB ----
  # Specify function parameters ----
  tanner_mod <- readRDS("./Maturity data processing/Doc/Tanner models/sdmTMB_spVAR_SIZE_k200.rda")
  tanner_dat <- readRDS("./Maturity data processing/Data/tanner_survey_specimenEBS.rda")
  tanner_yrs <- c(1990:2019, 2021:2025)
  species <- "TANNER"
  
  # Run function ----
  calc_maturepop_estimates(tanner_mod, tanner_dat, tanner_yrs, species) -> tanner.out
  
  # OGIVES ----
  # sdmTMB
  tanner_ogive_sdmTMB <- tanner.out$ogives %>%
    mutate(Estimator = "sdmTMB")
  
  # legacy using gam prop_mature interpolation of 10mm down to 5mm bins
  tanner_legacy_ogives <- rbind(read.csv("./Maturity data processing/Doc/tanE_pmat_5mminterp.csv") %>% mutate(DISTRICT = "E166"),
                                read.csv("./Maturity data processing/Doc/tanW_pmat_5mminterp.csv") %>% mutate(DISTRICT = "W166"))%>%
    mutate(SPECIES = "TANNER") %>%
    dplyr::select(!X) %>%
    rename(SIZE = SIZE_BIN) %>%
    mutate(Estimator = "Legacy")
  
  # Combine 
  ogive_plot_dat <- rbind(tanner_ogive_sdmTMB %>% rename(PROP_MATURE = PROP_MATURE_mean, SIZE = SIZE_5MM), tanner_legacy_ogives %>% mutate(PROP_MATURE_lo = NA, PROP_MATURE_hi = NA))
  
  
  # Plot
  ggplot(ogive_plot_dat %>% filter(DISTRICT == "W166" & !YEAR %in% c(2013, 2015)),  aes(SIZE, PROP_MATURE, color = Estimator))+
    geom_ribbon(aes(SIZE, ymin = PROP_MATURE_lo, ymax = PROP_MATURE_hi, fill = Estimator), color = NA, alpha = 0.35)+
    geom_line(linewidth = 1)+
    scale_color_manual(values = c(
      "Legacy"       = "darkgoldenrod",
      "sdmTMB"       = "cadetblue",name = ""))+
    scale_fill_manual(values = c(
      "Legacy"       = "darkgoldenrod",
      "sdmTMB"       = "cadetblue"), name = "", guide= "none")+
    facet_wrap(~YEAR)+
    geom_rug(ogive_plot_dat %>% filter(Estimator == "sdmTMB", DISTRICT == "W166" & !YEAR %in% c(2013, 2015)), mapping =aes(x = SIZE), sides = "b", inherit.aes = TRUE) +
    geom_hline(yintercept = 0.5, linetype = "dashed")+
    theme_bw()+
    ggtitle("Tanner West")+
    ylab("Proportion mature")+
    xlab("Carapace width (mm)")+
    theme(legend.position = "bottom", legend.direction = "horizontal",
          legend.text = element_text(size = 12),
          legend.title = element_blank(),
          axis.title = element_text(size = 12),
          strip.text = element_text(size = 12),
          axis.text = element_text(size = 10))
  
  ggsave("./Maturity data processing/Doc/tannerW166_ogive_comparison_legacy.sdmTMBwitherror.png", width  = 10, height = 10, units = "in")
  
  
  ggplot(ogive_plot_dat %>% filter(DISTRICT == "E166" & !YEAR %in% c(2011, 2013, 2015)),  aes(SIZE, PROP_MATURE, color = Estimator))+
           geom_ribbon(aes(SIZE, ymin = PROP_MATURE_lo, ymax = PROP_MATURE_hi, fill = Estimator), color = NA, alpha = 0.35)+
           geom_line(linewidth = 1)+
           scale_color_manual(values = c(
             "Legacy"       = "darkgoldenrod",
             "sdmTMB"       = "cadetblue",name = ""))+
           scale_fill_manual(values = c(
             "Legacy"       = "darkgoldenrod",
             "sdmTMB"       = "cadetblue"), name = "", guide= "none")+
           facet_wrap(~YEAR)+
           geom_rug(ogive_plot_dat %>% filter(Estimator == "sdmTMB", DISTRICT == "E166" & !YEAR %in% c(2011, 2013, 2015)), mapping =aes(x = SIZE), sides = "b", inherit.aes = TRUE) +
           geom_hline(yintercept = 0.5, linetype = "dashed")+
           theme_bw()+
           ggtitle("Tanner East")+
         ylab("Proportion mature")+
           xlab("Carapace width (mm)")+
           theme(legend.position = "bottom", legend.direction = "horizontal",
                 legend.text = element_text(size = 12),
                 legend.title = element_blank(),
                 axis.title = element_text(size = 12),
                 strip.text = element_text(size = 12),
                 axis.text = element_text(size = 10))
         
  ggsave("./Maturity data processing/Doc/tannerE166_ogive_comparison_legacy.sdmTMBwitherror.png", width  = 10, height = 10, units = "in")
         
  ggplot(tanner_legacy_ogives %>% filter(DISTRICT == "E166"), aes(SIZE, PROP_MATURE, color = Estimator, group = YEAR))+
    geom_line(linewidth = 1)+
    scale_color_manual(values = c(
      "Legacy"       = "darkgoldenrod"))+
    #facet_wrap(~YEAR)+
    geom_hline(yintercept = 0.5, linetype = "dashed")+
    theme_bw()+
    ggtitle("Tanner East")+
    #scale_x_continuous(breaks = seq(0, 175, by = 10))+
    ylab("Probability mature")+
    xlab("Carapace width (mm)")+
    theme(legend.position = "bottom", legend.direction = "horizontal",
          legend.text = element_text(size = 12),
          legend.title = element_blank(),
          axis.title = element_text(size = 12),
          strip.text = element_text(size = 12),
          axis.text = element_text(size = 10)) -> l.2
  
  ggplot(tanner_legacy_ogives %>% filter(DISTRICT == "W166"), aes(SIZE, PROP_MATURE, color = Estimator, group = YEAR))+
    geom_line(linewidth = 1)+
    scale_color_manual(values = c(
      "Legacy"       = "darkgoldenrod"))+
    #facet_wrap(~YEAR)+
    geom_hline(yintercept = 0.5, linetype = "dashed")+
    theme_bw()+
    ggtitle("Tanner West")+
    #scale_x_continuous(breaks = seq(0, 175, by = 10))+
    ylab("Probability mature")+
    xlab("Carapace width (mm)")+
    theme(legend.position = "bottom", legend.direction = "horizontal",
          legend.text = element_text(size = 12),
          legend.title = element_blank(),
          axis.title = element_text(size = 12),
          strip.text = element_text(size = 12),
          axis.text = element_text(size = 10)) -> l.3
  
  l.1/l.2/l.3 +
    plot_layout(guides = "collect", axis_titles = "collect") & 
    theme(
      legend.position = "bottom",
      axis.title.x = element_text(),   # keep x title only on bottom row
      axis.title.y = element_text()    # shared y title
    )
  
  ggsave("./Maturity data processing/Figures/legacy_ogives_CPT.png", width = 7, height = 8)
  # SAM ----
  # sdmTMB SAM
  tanner_SAM <- tanner.out$SAM %>%
                na.omit() %>%
             right_join(., expand.grid(YEAR = seq(min(.$YEAR), max(.$YEAR)), DISTRICT = c("E166", "W166"),
                              Estimator = "sdmTMB"))
  
  # Legacy SAM
  legacy_SAM <- read.csv("./Maturity data processing/Data/tanner_legacy_ogives.csv") %>% 
            filter(DISTRICT != "ALL") %>%
            dplyr::select(YEAR, DISTRICT, B_EST, B_SE) %>% 
            distinct() %>%
            rename(SAM = B_EST) %>%
            right_join(., expand.grid(YEAR = seq(min(.$YEAR), max(.$YEAR)), DISTRICT = c("E166", "W166"),
                                      Estimator = "Legacy")) %>%
            distinct() 
  
  
  # Plot
  ggplot()+
    geom_line(legacy_SAM, mapping = aes(YEAR, SAM, color = as.factor(1)), linewidth = 1)+
    geom_point(legacy_SAM, mapping = aes(YEAR, SAM, color = as.factor(1)), size = 2)+
    geom_errorbar(legacy_SAM, mapping = aes(YEAR, ymin = SAM-B_SE*1.96, ymax = SAM + B_SE*1.96, color = as.factor(1)), linewidth = 1)+
    geom_line(tanner_SAM, mapping = aes(YEAR, SAM_mean, color = as.factor(2)), linewidth = 1)+
    geom_point(tanner_SAM, mapping = aes(YEAR, SAM_mean, color = as.factor(2)), size = 2)+
    geom_errorbar(tanner_SAM, mapping = aes(YEAR, ymin = SAM_lo, ymax = SAM_hi, color = as.factor(2)), linewidth = 1)+
    scale_color_manual(values = c("darkgoldenrod", "cadetblue"), labels = c("Legacy", "sdmTMB"), name = "")+
    scale_fill_manual(values = c("darkgoldenrod", "cadetblue"), labels = c("Legacy", "sdmTMB"), name = "")+
    theme_bw()+
    facet_wrap(~DISTRICT, nrow = 2, scales = "free_y")+
    xlab("Year")+
    ylab("Size at 50% maturity (mm)")+
    theme(legend.position = "bottom", legend.direction = "horizontal",
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          legend.text = element_text(size = 12))
  
  ggsave("./Maturity data processing/Doc/tanner_SAM_comparison.png", width = 8, height = 6, units = "in")
  
  # Plot
  eSAM <- legacy_SAM %>% filter(DISTRICT == "E166")
  eSAM2 <- tanner_SAM %>% filter(DISTRICT == "E166")
  ggplot()+
    geom_line(eSAM, mapping = aes(YEAR, SAM, color = as.factor(1)), linewidth = 1)+
    geom_point(eSAM, mapping = aes(YEAR, SAM, color = as.factor(1)), size = 2)+
    geom_errorbar(eSAM, mapping = aes(YEAR, ymin = SAM-B_SE*1.96, ymax = SAM + B_SE*1.96, color = as.factor(1)), linewidth = 1)+
    geom_line(eSAM2, mapping = aes(YEAR, SAM_mean, color = as.factor(2)), linewidth = 1)+
    geom_point(eSAM2, mapping = aes(YEAR, SAM_mean, color = as.factor(2)), size = 2)+
    geom_errorbar(eSAM2, mapping = aes(YEAR, ymin = SAM_lo, ymax = SAM_hi, color = as.factor(2)), linewidth = 1)+
    scale_color_manual(values = c("darkgoldenrod", "cadetblue"), labels = c("Legacy", "sdmTMB"), name = "")+
    scale_fill_manual(values = c("darkgoldenrod", "cadetblue"), labels = c("Legacy", "sdmTMB"), name = "")+
    theme_bw()+
    ggtitle("Tanner East")+
    xlab("Year")+
    ylab("Size at 50% maturity (mm)")+
    theme(legend.position = "bottom", legend.direction = "horizontal",
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          legend.text = element_text(size = 12)) -> ls.2
  
  wSAM <- legacy_SAM %>% filter(DISTRICT == "W166")
  wSAM2 <- tanner_SAM %>% filter(DISTRICT == "W166")
  ggplot()+
    geom_line(wSAM, mapping = aes(YEAR, SAM, color = as.factor(1)), linewidth = 1)+
    geom_point(wSAM, mapping = aes(YEAR, SAM, color = as.factor(1)), size = 2)+
    geom_errorbar(wSAM, mapping = aes(YEAR, ymin = SAM-B_SE*1.96, ymax = SAM + B_SE*1.96, color = as.factor(1)), linewidth = 1)+
    geom_line(wSAM2, mapping = aes(YEAR, SAM_mean, color = as.factor(2)), linewidth = 1)+
    geom_point(wSAM2, mapping = aes(YEAR, SAM_mean, color = as.factor(2)), size = 2)+
    geom_errorbar(wSAM2, mapping = aes(YEAR, ymin = SAM_lo, ymax = SAM_hi, color = as.factor(2)), linewidth = 1)+
    scale_color_manual(values = c("darkgoldenrod", "cadetblue"), labels = c("Legacy", "sdmTMB"), name = "")+
    scale_fill_manual(values = c("darkgoldenrod", "cadetblue"), labels = c("Legacy", "sdmTMB"), name = "")+
    theme_bw()+
    ggtitle("Tanner West")+
    xlab("Year")+
    ylab("Size at 50% maturity (mm)")+
    theme(legend.position = "bottom", legend.direction = "horizontal",
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          legend.text = element_text(size = 12)) -> ls.3
  
  
  ls.1/ls.2/ls.3 + plot_layout(guides = "collect", axes = "collect") & theme(legend.position = "bottom", legend.direction = "horizontal")
  
  ggsave("./Maturity data processing/Figures/legacysdmTMB_SAM_CPT.png", width = 7, height = 8)
  
  
  eSAM <- legacy_SAM %>% filter(DISTRICT == "E166")
  ggplot()+
    geom_line(eSAM, mapping = aes(YEAR, SAM, color = as.factor(1)), linewidth = 1)+
    geom_point(eSAM, mapping = aes(YEAR, SAM, color = as.factor(1)), size = 2)+
    geom_errorbar(eSAM, mapping = aes(YEAR, ymin = SAM-B_SE*1.96, ymax = SAM + B_SE*1.96, color = as.factor(1)), linewidth = 1)+
    scale_color_manual(values = c("darkgoldenrod"), labels = c("Legacy"), name = "")+
    scale_fill_manual(values = c("darkgoldenrod"), labels = c("Legacy"), name = "")+
    theme_bw()+
    ggtitle("Tanner East")+
    xlab("Year")+
    ylab("Size at 50% maturity (mm)")+
    theme(legend.position = "bottom", legend.direction = "horizontal",
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          legend.text = element_text(size = 12)) -> ls.2
  
  wSAM <- legacy_SAM %>% filter(DISTRICT == "W166")
  ggplot()+
    geom_line(wSAM, mapping = aes(YEAR, SAM, color = as.factor(1)), linewidth = 1)+
    geom_point(wSAM, mapping = aes(YEAR, SAM, color = as.factor(1)), size = 2)+
    geom_errorbar(wSAM, mapping = aes(YEAR, ymin = SAM-B_SE*1.96, ymax = SAM + B_SE*1.96, color = as.factor(1)), linewidth = 1)+
    scale_color_manual(values = c("darkgoldenrod"), labels = c("Legacy"), name = "")+
    scale_fill_manual(values = c("darkgoldenrod"), labels = c("Legacy"), name = "")+
    theme_bw()+
    ggtitle("Tanner West")+
    xlab("Year")+
    ylab("Size at 50% maturity (mm)")+
    theme(legend.position = "bottom", legend.direction = "horizontal",
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          legend.text = element_text(size = 12)) -> ls.3
  
  ls.1/ls.2/ls.3 + plot_layout(guides = "collect", axes = "collect") & theme(legend.position = "bottom")
  ggsave("./Maturity data processing/Figures/legacy_SAM_CPT.png", width = 7, height = 8)
  
  
  # MATURE BIOMASS/ABUNDANCE ----
  # sdmTMB
  tanner_sdmTMB_bioabund <- tanner.out$mature_bioabund
  
  # Legacy
  tanner_legacy_spec <- tanner_dat
  tanner_legacy_spec$specimen <- tanner_dat$specimen %>%
    filter(YEAR %in% tanner_yrs, SPECIES == species, SHELL_CONDITION == 2, SEX == 1) %>%
    mutate(BIN_5MM = cut_width(SIZE_1MM, width = 5, center = 2.5, closed = "left", dig.lab = 4),
           BIN2 = BIN_5MM) %>%
    separate(BIN2, sep = ",", into = c("LOWER", "UPPER")) %>%
    mutate(LOWER = as.numeric(sub('.', '', LOWER)),
           UPPER = as.numeric(gsub('.$', '', UPPER)),
           SIZE_5MM = (UPPER + LOWER)/2) %>%
    dplyr::select(!c(BIN_5MM, LOWER, UPPER)) 
  
  tanner_pars <- read.csv("./Maturity data processing/Data/tanner_legacy_modelpars.csv") 
  
  tanner_legacy_spec$specimen <-  tanner_legacy_spec$specimen %>% 
    right_join(., tanner_pars) %>%
    mutate(PROP_MATURE = (1/(1 + exp(-A_EST * (SIZE_5MM - B_EST))))) %>%
    mutate(SAMPLING_FACTOR = SAMPLING_FACTOR * PROP_MATURE) %>%
    dplyr::select(!X)
  
  # Use crabpack to calculate bioabund
  tanner_legacy_bioabund <-  crabpack::calc_bioabund(crab_data = tanner_legacy_spec, species = "TANNER", 
                                                   size_min = NULL, size_max = NULL,  sex = "male", 
                                                   shell_condition = "new_hardshell") %>%
    right_join(., data.frame(YEAR = seq(min(.$YEAR), max(.$YEAR)))) %>%
    mutate(Estimator = "Legacy",
           ABUNDANCE = ABUNDANCE/1e6,
           ABUNDANCE_CI = ABUNDANCE_CI/1e6)
  
  # Bind and specify missing years
  tanner_bioabund_dat <- rbind(tanner_legacy_bioabund%>%
                            dplyr::select(names(.)[names(.) %in% colnames(tanner_sdmTMB_bioabund)]), tanner_sdmTMB_bioabund)
  
  # Plot
  e1 <- ggplot(tanner_bioabund_dat %>% filter(DISTRICT == "E166"), aes(YEAR, ABUNDANCE))+
    geom_point(tanner_bioabund_dat %>% filter(DISTRICT == "E166", YEAR %in% c(2012, 2014)),
               mapping = aes(color = Estimator))+
    geom_errorbar(tanner_bioabund_dat %>% filter(DISTRICT == "E166", YEAR %in% c(2012, 2014)),
                  mapping = aes(x = YEAR, ymin = ABUNDANCE - ABUNDANCE_CI, ymax = ABUNDANCE + ABUNDANCE_CI, 
                                color= Estimator), linewidth = 1)+
    geom_ribbon(mapping = aes(x = YEAR, ymin = ABUNDANCE - ABUNDANCE_CI, ymax = ABUNDANCE + ABUNDANCE_CI, 
                              fill= Estimator), alpha = 0.15)+
    geom_line(mapping = aes(color = Estimator), linewidth = 0.75)+
    scale_color_manual(values = c("darkgoldenrod", "cadetblue"), name = "")+
    scale_fill_manual(values = c("darkgoldenrod", "cadetblue"), name = "")+
    theme_bw()+
    ylab("Abundance (millions)")+
    ggtitle("Tanner East")+
    xlab("Year")+
    theme(
          axis.text.y = element_text(size = 12), 
          axis.title.y = element_text(size = 12),
          legend.position = "bottom",
          legend.direction = "horizontal")
  
  e2 <- ggplot(tanner_bioabund_dat %>% filter(DISTRICT == "E166"), aes(YEAR, BIOMASS_MT))+
    geom_point(tanner_bioabund_dat %>% filter(DISTRICT == "E166", YEAR %in% c(2012, 2014)),
               mapping = aes(color = Estimator))+
    geom_errorbar(tanner_bioabund_dat %>% filter(DISTRICT == "E166", YEAR %in% c(2012, 2014)),
                  mapping = aes(x = YEAR, ymin = BIOMASS_MT - BIOMASS_MT_CI, ymax = BIOMASS_MT + BIOMASS_MT_CI, 
                                color= Estimator), linewidth = 1)+
    geom_ribbon(mapping = aes(x = YEAR, ymin = BIOMASS_MT - BIOMASS_MT_CI, ymax = BIOMASS_MT + BIOMASS_MT_CI, 
                              fill= Estimator), alpha = 0.15)+
    geom_line(mapping = aes(color = Estimator), linewidth = 0.75)+
    scale_color_manual(values = c("darkgoldenrod", "cadetblue"), name = "")+
    scale_fill_manual(values = c("darkgoldenrod", "cadetblue"), name = "")+
    theme_bw()+
    ylab("Biomass (tons)")+
    #ggtitle("Morphometric mature male biomass (SH2)")+
    xlab("Year")+
    theme(axis.text = element_text(size = 12),
          legend.position = "bottom",
          legend.direction = "horizontal",
          axis.title = element_text(size = 12))
  
  e1/e2 + plot_layout(guides = "collect") & 
    theme(legend.position = "bottom", legend.text = element_text(size = 12))
  
  ggsave("./Maturity data processing/Doc/tannerE_bioabund_comparison_legacy.sdmTMB.png", width = 7, height = 6, units = "in")
  
  w1 <- ggplot(tanner_bioabund_dat %>% filter(DISTRICT == "W166"), aes(YEAR, ABUNDANCE))+
    geom_point(tanner_bioabund_dat %>% filter(DISTRICT == "W166", YEAR %in% c(2014)),
               mapping = aes(color = Estimator))+
    geom_errorbar(tanner_bioabund_dat %>% filter(DISTRICT == "W166", YEAR %in% c(2014)),
                  mapping = aes(x = YEAR, ymin = ABUNDANCE - ABUNDANCE_CI, ymax = ABUNDANCE + ABUNDANCE_CI, 
                                color= Estimator), linewidth = 1)+
    geom_ribbon(mapping = aes(x = YEAR, ymin = ABUNDANCE - ABUNDANCE_CI, ymax = ABUNDANCE + ABUNDANCE_CI, 
                              fill= Estimator), alpha = 0.15)+
    geom_line(mapping = aes(color = Estimator), linewidth = 0.75)+
    scale_color_manual(values = c("darkgoldenrod", "cadetblue"), name = "")+
    scale_fill_manual(values = c("darkgoldenrod", "cadetblue"), name = "")+
    theme_bw()+
    ylab("Abundance (millions)")+
    ggtitle("Tanner West")+
    xlab("Year")+
    theme(axis.text.y = element_text(size = 12), 
          axis.title.y = element_text(size = 12),
          legend.position = "bottom",
          legend.direction = "horizontal")
  
  w2 <- ggplot(tanner_bioabund_dat %>% filter(DISTRICT == "W166"), aes(YEAR, BIOMASS_MT))+
    geom_point(tanner_bioabund_dat %>% filter(DISTRICT == "W166", YEAR %in% c(2014)),
               mapping = aes(color = Estimator))+
    geom_errorbar(tanner_bioabund_dat %>% filter(DISTRICT == "W166", YEAR %in% c(2014)),
                  mapping = aes(x = YEAR, ymin = BIOMASS_MT - BIOMASS_MT_CI, ymax = BIOMASS_MT + BIOMASS_MT_CI, 
                                color= Estimator), linewidth = 1)+
    geom_ribbon(mapping = aes(x = YEAR, ymin = BIOMASS_MT - BIOMASS_MT_CI, ymax = BIOMASS_MT + BIOMASS_MT_CI, 
                              fill= Estimator), alpha = 0.15)+
    geom_line(mapping = aes(color = Estimator), linewidth = 0.75)+
    scale_color_manual(values = c("darkgoldenrod", "cadetblue"), name = "")+
    scale_fill_manual(values = c("darkgoldenrod", "cadetblue"), name = "")+
    theme_bw()+
    ylab("Biomass (tons)")+
    #ggtitle("Morphometric mature male biomass (SH2)")+
    xlab("Year")+
    theme(axis.text = element_text(size = 12),
          legend.position = "bottom",
          legend.direction = "horizontal",
          axis.title = element_text(size = 12))
  
  w1/w2 + plot_layout(guides = "collect") & 
    theme(legend.position = "bottom", legend.text = element_text(size = 12))
  
  ggsave("./Maturity data processing/Doc/tannerW_bioabund_comparison_legacy.sdmTMB.png", width = 7, height = 6, units = "in")
  
  
  p1/e1/w1 + plot_layout(guides = "collect", axes = "collect") & theme(legend.position = "bottom", legend.text=element_text(size = 12))
  ggsave("./Maturity data processing/Figures/legacysdmTMB_bioabund_CPT.png", width = 7, height = 8)
  
  
  p1 <- ggplot()+
    geom_point(snow_bioabund_dat %>% filter(YEAR %in% c(2013, 2015), Estimator == "Legacy"), mapping = aes(YEAR, ABUNDANCE, color = Estimator))+
    geom_errorbar(snow_bioabund_dat %>% filter(YEAR %in% c(2013, 2015), Estimator == "Legacy"), mapping = aes(x = YEAR, ymin = ABUNDANCE - ABUNDANCE_CI, ymax = ABUNDANCE + ABUNDANCE_CI, 
                                                                                       color= Estimator), linewidth = 1)+
    geom_ribbon(snow_bioabund_dat %>% filter(Estimator == "Legacy"), mapping = aes(x = YEAR, ymin = ABUNDANCE - ABUNDANCE_CI, ymax = ABUNDANCE + ABUNDANCE_CI, 
                                                 fill= Estimator), alpha = 0.15)+
    geom_line(snow_bioabund_dat %>% filter(Estimator == "Legacy"), mapping = aes(YEAR, ABUNDANCE, color = Estimator), linewidth = 0.75)+
    scale_color_manual(values = c("darkgoldenrod", "cadetblue"), name = "")+
    scale_fill_manual(values = c("darkgoldenrod", "cadetblue"), name = "")+
    theme_bw()+
    ylab("Abundance (millions)")+
    ggtitle("Snow")+
    xlab("Year")+
    theme(axis.text.y = element_text(size = 12), 
          axis.title.y = element_text(size = 12),
          legend.position = "bottom",
          legend.direction = "horizontal")
  
  
  e1 <- ggplot(tanner_bioabund_dat %>% filter(DISTRICT == "E166", Estimator =="Legacy"), aes(YEAR, ABUNDANCE))+
    geom_point(tanner_bioabund_dat %>% filter(DISTRICT == "E166", Estimator =="Legacy", YEAR %in% c(2012, 2014)),
               mapping = aes(color = Estimator))+
    geom_errorbar(tanner_bioabund_dat %>% filter(DISTRICT == "E166", Estimator =="Legacy", YEAR %in% c(2012, 2014)),
                  mapping = aes(x = YEAR, ymin = ABUNDANCE - ABUNDANCE_CI, ymax = ABUNDANCE + ABUNDANCE_CI, 
                                color= Estimator), linewidth = 1)+
    geom_ribbon(mapping = aes(x = YEAR, ymin = ABUNDANCE - ABUNDANCE_CI, ymax = ABUNDANCE + ABUNDANCE_CI, 
                              fill= Estimator), alpha = 0.15)+
    geom_line(mapping = aes(color = Estimator), linewidth = 0.75)+
    scale_color_manual(values = c("darkgoldenrod", "cadetblue"), name = "")+
    scale_fill_manual(values = c("darkgoldenrod", "cadetblue"), name = "")+
    theme_bw()+
    ylab("Abundance (millions)")+
    ggtitle("Tanner East")+
    xlab("Year")+
    theme(axis.text.y = element_text(size = 12), 
          axis.title.y = element_text(size = 12),
          legend.position = "bottom",
          legend.direction = "horizontal")
  
  w1 <- ggplot(tanner_bioabund_dat %>% filter(DISTRICT == "W166", Estimator == "Legacy"), aes(YEAR, ABUNDANCE))+
    geom_point(tanner_bioabund_dat %>% filter(DISTRICT == "W166",Estimator == "Legacy", YEAR %in% c(2014)),
               mapping = aes(color = Estimator))+
    geom_errorbar(tanner_bioabund_dat %>% filter(DISTRICT == "W166", Estimator == "Legacy", YEAR %in% c(2014)),
                  mapping = aes(x = YEAR, ymin = ABUNDANCE - ABUNDANCE_CI, ymax = ABUNDANCE + ABUNDANCE_CI, 
                                color= Estimator), linewidth = 1)+
    geom_ribbon(mapping = aes(x = YEAR, ymin = ABUNDANCE - ABUNDANCE_CI, ymax = ABUNDANCE + ABUNDANCE_CI, 
                              fill= Estimator), alpha = 0.15)+
    geom_line(mapping = aes(color = Estimator), linewidth = 0.75)+
    scale_color_manual(values = c("darkgoldenrod", "cadetblue"), name = "")+
    scale_fill_manual(values = c("darkgoldenrod", "cadetblue"), name = "")+
    theme_bw()+
    ylab("Abundance (millions)")+
    ggtitle("Tanner West")+
    xlab("Year")+
    theme( axis.text.y = element_text(size = 12), 
          axis.title.y = element_text(size = 12),
          legend.position = "bottom",
          legend.direction = "horizontal")
  
  p1/e1/w1 + plot_layout(guides = "collect", axes = "collect") & theme(legend.position = "bottom", 
                                                                       axis.text.x = element_text(size = 12),
                                                                       axis.title.x = element_text(size = 12),
                                                                       legend.text = element_text(size = 12))
  
  
  ggsave("./Maturity data processing/Figures/legacy_bioabund_CPT.png", width = 7, height = 8)
  
  
  
  
  
  
  
  
  
  
  
  
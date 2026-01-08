# PURPOSE: to create function to calculate mature cpue, biomass, and abundance, using different estimators for expanding
# chela subsample to specimen subsample

# Author: Emily Ryznar

# NOTES:

# SPECIFY FUNCTIONS --------------------------------------------------------------------------------------
calc_mature_cpue <- function(crab_data, species, years, size_min, size_max, sex, shell_condition, estimator){
  
  # Add in 5mm bins to crab_data
  crab_data$specimen <-  crab_data$specimen%>%
    mutate(BIN_5MM = cut_width(SIZE_1MM, width = 5, center = 2.5, closed = "left", dig.lab = 4),
           BIN2 = BIN_5MM) %>%
    separate(BIN2, sep = ",", into = c("LOWER", "UPPER")) %>%
    mutate(LOWER = as.numeric(sub('.', '', LOWER)),
           UPPER = as.numeric(gsub('.$', '', UPPER)),
           SIZE_5MM = (UPPER + LOWER)/2) %>%
    dplyr::select(!c(BIN_5MM, LOWER, UPPER)) %>%
    mutate(BIN_10MM = cut_width(SIZE_1MM, width = 10, center = 5, closed = "left", dig.lab = 4),
           BIN2 = BIN_10MM) %>%
    separate(BIN2, sep = ",", into = c("LOWER", "UPPER")) %>%
    mutate(LOWER = as.numeric(sub('.', '', LOWER)),
           UPPER = as.numeric(gsub('.$', '', UPPER)),
           SIZE_10MM = (UPPER + LOWER)/2) %>%
    dplyr::select(!c(BIN_10MM, LOWER, UPPER))
  
  # Specify corresponding sex code
  if(sex == "male"){
    sx_num = 1
  } else if(sex == "female"){
    sx_num = 2
  }
  
  # Specify corresponding SC codes
  if(shell_condition == "soft molting"){
    sc_num <- 1
  } else if(shell_condition == "new hardshell"){
    sc_num <- 2
  } else if(shell_condition == "oldshell"){
    sc_num <- 3
  } else if(shell_condition == "very oldshell"){
    sc_num <- 4
  }
  
  # Specify estimators
  if(estimator == "MWK"){
    mwk <- read.csv("./Maturity research/Output/MWK.csv")
    
    # Filter specimen data by params (not size yet for full join), join with mwk
    crab_data$specimen <- crab_data$specimen %>% 
      filter(YEAR %in% years, SPECIES %in% species, SHELL_CONDITION == sc_num, SEX == sx_num) %>%
      right_join(., mwk, by = c("YEAR", "STATION_ID", "LATITUDE", "LONGITUDE", "SIZE_1MM", "SIZE_5MM")) %>%
      filter(!YEAR %in% 2012) %>%
      mutate(SAMPLING_FACTOR = SAMPLING_FACTOR * PROP_MATURE)
    
  } else if(estimator == "GAM"){
    gam_mods <- readRDS("./Maturity research/Output/gam_mods.rda")
    
    gam_preds <- data.frame()
    
    for(ii in 1:length(years)){
      print(paste("Predicting GAM", years[ii]))
      
      # select year model
      mod <- gam_mods[[as.character(years[ii])]]
      
      # filter specimen data by year
      crab_data$specimen %>% filter(YEAR == years[ii]) -> sub1
      
      # predict model
      pp <- cbind(sub1, predict(mod, sub1, , type = "response", se = TRUE))
      
      # bind predictions across year
      gam_preds <- rbind(gam_preds, pp) 
      
    }
    
    # Filter predicted specimen data by params (not size yet for full join)
    crab_data$specimen <- gam_preds %>% 
      rename(PROP_MATURE = fit) %>%
      filter(YEAR %in% years, SPECIES %in% species, SHELL_CONDITION == sc_num, SEX == sx_num) %>%
      mutate(SAMPLING_FACTOR = SAMPLING_FACTOR * PROP_MATURE)
    
  } else if(estimator == "spatial GAM"){
    gam_mods <- readRDS("./Maturity research/Output/gam_mods_spatial.rda")
    
    gam_preds <- data.frame()
    
    for(ii in 1:length(years)){
      print(paste("Predicting spatial GAM", years[ii]))
      
      # select year model
      mod <- gam_mods[[as.character(years[ii])]]
      
      # filter specimen data by year
      crab_data$specimen %>% filter(YEAR == years[ii]) -> sub1
      
      # predict model
      pp <- cbind(sub1, predict(mod, sub1, , type = "response", se = TRUE))
      
      # bind predictions across year
      gam_preds <- rbind(gam_preds, pp) 
      
    }
    
    # Filter predicted specimen data by params (not size yet for full join)
    crab_data$specimen <- gam_preds %>% 
      rename(PROP_MATURE = fit) %>%
      filter(YEAR %in% years, SPECIES %in% species, SHELL_CONDITION == sc_num, SEX == sx_num) %>%
      mutate(SAMPLING_FACTOR = SAMPLING_FACTOR * PROP_MATURE)
    
  } else if(estimator == "pooled GAM"){
    print("Predicting pooled GAM")
    
    mod <- readRDS("./Maturity research/Output/gam_pooled.rda")
    
    # filter specimen data by year
    crab_data$specimen -> sub1
    
    # predict model
    gam_preds <- cbind(sub1, predict(mod, sub1, , type = "response", se = TRUE))
    
    
    # Filter predicted specimen data by params (not size yet for full join)
    crab_data$specimen <- gam_preds %>% 
      rename(PROP_MATURE = fit) %>%
      filter(YEAR %in% years, SPECIES %in% species, SHELL_CONDITION == sc_num, SEX == sx_num) %>%
      mutate(SAMPLING_FACTOR = SAMPLING_FACTOR * PROP_MATURE)
  } else if(estimator == "tinyVAST"){
    print("Predicting tinyVAST")
    
    
    mod <- readRDS("./Maturity research/Output/tinyVAST/tinyVAST_fullsize_correlation.rda")
    
    # filter specimen data by year
    crab_data$specimen -> sub1
    
    # predict model
    preds <- cbind(sub1, fit = predict(mod, sub1, , type = "response"))
    
    
    # Filter predicted specimen data by params (not size yet for full join)
    crab_data$specimen <-preds %>% 
      rename(PROP_MATURE = fit) %>%
      filter(YEAR %in% years, SPECIES %in% species, SHELL_CONDITION == sc_num, SEX == sx_num) %>%
      mutate(SAMPLING_FACTOR = SAMPLING_FACTOR * PROP_MATURE)
  } else if(estimator == "sdmTMB"){
    print("Predicting sdmTMB")
    
    mod <- readRDS("./Maturity research/Output/sdmTMB/s(SIZE)_spvar_iid_90_sdmTMB.rda")
    
    # filter specimen data by year and transform to sdmTMB coordinates
    crab_data$specimen %>%
      st_as_sf(., coords = c("LONGITUDE", "LATITUDE"), crs = "+proj=longlat +datum=WGS84") %>%
      st_transform(., crs = "+proj=utm +zone=2") %>%
      cbind(st_coordinates(.)) %>%
      as.data.frame(.) %>%
      mutate(LATITUDE = Y/1000, # scale to km so values don't get too large
             LONGITUDE = X/1000,
             YEAR_F = as.factor(YEAR),
             YEAR_SCALED = as.numeric(scale(mod.dat$YEAR))) -> sub1
    
    # predict model
    pmat.sim <- predict(mod, sub1, type = "response", nsim = 500)
    
    fit <- apply(pmat.sim, 1, mean)
    se_fit   <- apply(pmat.sim, 1, sd)
    
    preds <- cbind(sub1, fit, se_fit)
    
    
    # Filter predicted specimen data by params (not size yet for full join)
    crab_data$specimen <-preds %>% 
      rename(PROP_MATURE = fit) %>%
      filter(YEAR %in% years, SPECIES %in% species, SHELL_CONDITION == sc_num, SEX == sx_num) %>%
      mutate(SAMPLING_FACTOR = SAMPLING_FACTOR * PROP_MATURE)
  }
  
  # Calculate mature cpue by params
  mat_cpue <- crabpack::calc_cpue(crab_data = crab_data, years = years, size_min = size_min, size_max = size_max, 
                                  species = species, sex = sex, shell_condition = shell_condition)
  
  # Return
  return(mat_cpue)
}

calc_mature_bioabund <- function(crab_data, species, years, size_min, size_max, sex, shell_condition, estimator){
  
  # Add in 5mm bins to crab_data
  crab_data$specimen <-  crab_data$specimen%>%
    mutate(BIN_5MM = cut_width(SIZE_1MM, width = 5, center = 2.5, closed = "left", dig.lab = 4),
           BIN2 = BIN_5MM) %>%
    separate(BIN2, sep = ",", into = c("LOWER", "UPPER")) %>%
    mutate(LOWER = as.numeric(sub('.', '', LOWER)),
           UPPER = as.numeric(gsub('.$', '', UPPER)),
           SIZE_5MM = (UPPER + LOWER)/2) %>%
    dplyr::select(!c(BIN_5MM, LOWER, UPPER)) %>%
    mutate(BIN_10MM = cut_width(SIZE_1MM, width = 10, center = 5, closed = "left", dig.lab = 4),
           BIN2 = BIN_10MM) %>%
    separate(BIN2, sep = ",", into = c("LOWER", "UPPER")) %>%
    mutate(LOWER = as.numeric(sub('.', '', LOWER)),
           UPPER = as.numeric(gsub('.$', '', UPPER)),
           SIZE_10MM = (UPPER + LOWER)/2) %>%
    dplyr::select(!c(BIN_10MM, LOWER, UPPER))
  
  # Specify corresponding sex code
  if(sex == "male"){
    sx_num = 1
  } else if(sex == "female"){
    sx_num = 2
  }
  
  # Specify corresponding SC codes
  if(shell_condition == "soft molting"){
    sc_num <- 1
  } else if(shell_condition == "new hardshell"){
    sc_num <- 2
  } else if(shell_condition == "oldshell"){
    sc_num <- 3
  } else if(shell_condition == "very oldshell"){
    sc_num <- 4
  }
  
  # Specify estimators
  if(estimator == "MWK"){
    mwk <- read.csv("./Maturity research/Output/MWK.csv")
    
    # Filter specimen data by params (not size yet for full join), join with mwk
    crab_data$specimen <- crab_data$specimen %>% 
      filter(YEAR %in% years, SPECIES %in% species, SHELL_CONDITION == sc_num, SEX == sx_num) %>%
      right_join(., mwk, by = c("YEAR", "STATION_ID", "LATITUDE", "LONGITUDE", "SIZE_1MM", "SIZE_5MM")) %>%
      filter(!YEAR %in% 2012) %>%
      mutate(SAMPLING_FACTOR = SAMPLING_FACTOR * PROP_MATURE)
    
  } else if(estimator == "GAM"){
    gam_mods <- readRDS("./Maturity research/Output/gam_mods.rda")
    
    gam_preds <- data.frame()
    
    for(ii in 1:length(years)){
      print(paste("Predicting GAM", years[ii]))
      
      # select year model
      mod <- gam_mods[[as.character(years[ii])]]
      
      # filter specimen data by year
      crab_data$specimen %>% filter(YEAR == years[ii]) -> sub1
      
      # predict model
      pp <- cbind(sub1, predict(mod, sub1, , type = "response", se = TRUE))
      
      # bind predictions across year
      gam_preds <- rbind(gam_preds, pp) 
      
    }
    
    # Filter predicted specimen data by params (not size yet for full join)
    crab_data$specimen <- gam_preds %>% 
      rename(PROP_MATURE = fit) %>%
      filter(YEAR %in% years, SPECIES %in% species, SHELL_CONDITION == sc_num, SEX == sx_num) %>%
      mutate(SAMPLING_FACTOR = SAMPLING_FACTOR * PROP_MATURE)
    
  } else if(estimator == "spatial GAM"){
    gam_mods <- readRDS("./Maturity research/Output/gam_mods_spatial.rda")
    
    gam_preds <- data.frame()
    
    for(ii in 1:length(years)){
      print(paste("Predicting spatial GAM", years[ii]))
      
      # select year model
      mod <- gam_mods[[as.character(years[ii])]]
      
      # filter specimen data by year
      crab_data$specimen %>% filter(YEAR == years[ii]) -> sub1
      
      # predict model
      pp <- cbind(sub1, predict(mod, sub1, , type = "response", se = TRUE))
      
      # bind predictions across year
      gam_preds <- rbind(gam_preds, pp) 
      
    }
    
    # Filter predicted specimen data by params (not size yet for full join)
    crab_data$specimen <- gam_preds %>% 
      rename(PROP_MATURE = fit) %>%
      filter(YEAR %in% years, SPECIES %in% species, SHELL_CONDITION == sc_num, SEX == sx_num) %>%
      mutate(SAMPLING_FACTOR = SAMPLING_FACTOR * PROP_MATURE)
  } else if(estimator == "pooled GAM"){
    print("Predicting pooled GAM")
    
    mod <- readRDS("./Maturity research/Output/gam_pooled.rda")
    
    # filter specimen data by year
    crab_data$specimen -> sub1
    
    # predict model
    gam_preds <- cbind(sub1, predict(mod, sub1, , type = "response", se = TRUE))
    
    
    # Filter predicted specimen data by params (not size yet for full join)
    crab_data$specimen <- gam_preds %>% 
      rename(PROP_MATURE = fit) %>%
      filter(YEAR %in% years, SPECIES %in% species, SHELL_CONDITION == sc_num, SEX == sx_num) %>%
      mutate(SAMPLING_FACTOR = SAMPLING_FACTOR * PROP_MATURE)
    
  } else if(estimator == "tinyVAST"){
    print("Predicting tinyVAST")
    
    
    mod <- readRDS("./Maturity research/Output/tinyVAST/tinyVAST_fullsize_correlation.rda")
    
    # filter specimen data by year
    crab_data$specimen -> sub1
    
    # predict model
    preds <- cbind(sub1, fit = predict(mod, sub1, , type = "response"))
    
    
    # Filter predicted specimen data by params (not size yet for full join)
    crab_data$specimen <-preds %>% 
      rename(PROP_MATURE = fit) %>%
      filter(YEAR %in% years, SPECIES %in% species, SHELL_CONDITION == sc_num, SEX == sx_num) %>%
      mutate(SAMPLING_FACTOR = SAMPLING_FACTOR * PROP_MATURE)
    
  } else if(estimator == "sdmTMB"){
    print("Predicting sdmTMB")
    
    mod <- readRDS("./Maturity research/Output/sdmTMB/s(SIZE)_spvar_iid_90_sdmTMB.rda")

    # filter specimen data by year and transform to sdmTMB coordinates
    crab_data$specimen %>%
      st_as_sf(., coords = c("LONGITUDE", "LATITUDE"), crs = "+proj=longlat +datum=WGS84") %>%
      st_transform(., crs = "+proj=utm +zone=2") %>%
      cbind(st_coordinates(.)) %>%
      as.data.frame(.) %>%
      mutate(LATITUDE = Y/1000, # scale to km so values don't get too large
             LONGITUDE = X/1000,
             YEAR_F = as.factor(YEAR),
             YEAR_SCALED = as.numeric(scale(mod.dat$YEAR))) -> sub1
    
    # predict model
    pmat.sim <- predict(mod, sub1, type = "response", nsim = 500)
    
    fit <- apply(pmat.sim, 1, mean)
    se_fit   <- apply(pmat.sim, 1, sd)
    
    preds <- cbind(sub1, fit, se_fit)
    
    
    # Filter predicted specimen data by params (not size yet for full join)
    crab_data$specimen <-preds %>% 
      rename(PROP_MATURE = fit) %>%
      filter(YEAR %in% years, SPECIES %in% species, SHELL_CONDITION == sc_num, SEX == sx_num) %>%
      mutate(SAMPLING_FACTOR = SAMPLING_FACTOR * PROP_MATURE)
  }
  
  
  
  # Calculate mature bioabund by params
  mat_bioabund <- crabpack::calc_bioabund(crab_data = crab_data, years = years, size_min = size_min, size_max = size_max,
                                          species = "SNOW", sex = sex, shell_condition = shell_condition)
  
  # Return
  return(list(ogive_dat = crab_data$specimen, mature_bioabund = mat_bioabund))
}


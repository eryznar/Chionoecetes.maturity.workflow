# PURPOSE: 
# To get updated survey data by year for chela maturity processing

# Notes:
# 1) What is the best way to update data? Use the chela database and then specimen data each new year? Or will
# the chela database be updated in time?

# LOAD LIBS/PARAMS -----
source("./Maturity data processing/Scripts/1) load_libs_params.R")

# Set channel
channel <- "API"

# GET SNOW SURVEY SPECIMEN/CPUE/BIOABUND DATA FROM CRABPACK ----
# Pull specimen data
species <- "SNOW"
specimen_data <- crabpack::get_specimen_data(species = species,
                                             region = "EBS",
                                             #years = years, # set in libs/params script
                                             channel = channel)

saveRDS(specimen_data, "./Maturity data processing/Data/snow_survey_specimenEBS.rda")

# Calculate per-station CPUE by 1mm bins for weighting for ALL shell 2 males, not just chela msrd (already filters ! HT17)
snow_cpue <- calc_cpue(crab_data = readRDS("./Maturity data processing/Data/snow_survey_specimenEBS.rda"),
                       species = species,
                       years = years, # set in libs/params script
                       sex = "male",
                       shell_condition = "new_hardshell",
                       bin_1mm = TRUE) %>%
  dplyr::select(YEAR, SPECIES, STATION_ID, LATITUDE, LONGITUDE, SIZE_1MM, CPUE, CPUE_MT)

saveRDS(snow_cpue, "./Maturity data processing/Data/snow_survey_CPUE_male1mm.rda")

# Calculate biomass and abundance for small/large size bins 
snow_bioabund_small <- calc_bioabund(crab_data = readRDS("./Maturity data processing/Data/snow_survey_specimenEBS.rda"),
                                     species = species,
                                     years = years, # set in libs/params script
                                     sex = "male",
                                     shell_condition = "new_hardshell",
                                     size_min = 55,
                                     size_max = 65)

write.csv(snow_bioabund_small, "./Maturity research/Data/snow_bioabund_smallbin.csv")

snow_bioabund_large <- calc_bioabund(crab_data = readRDS("./Maturity data processing/Data/snow_survey_specimenEBS.rda"),
                                     species = species,
                                     years = years, # set in libs/params script
                                     sex = "male",
                                     shell_condition = "new_hardshell",
                                     size_min = 95,
                                     size_max = 105)

write.csv(snow_bioabund_large, "./Maturity research/Data/snow_bioabund_largebin.csv")


# GET TANNER SURVEY SPECIMEN/CPUE/BIOABUND DATA FROM CRABPACK ----
# Pull specimen data
species <- "TANNER"
specimen_data <- crabpack::get_specimen_data(species = species,
                                             region = "EBS",
                                             #years = years, # set in libs/params script
                                             channel = channel)

saveRDS(specimen_data, "./Maturity data processing/Data/tanner_survey_specimenEBS.rda")

# Calculate per-station CPUE by 1mm bins for weighting for ALL shell 2 males, not just chela msrd (already filters ! HT17)
tanner_cpue <- calc_cpue(crab_data = readRDS("./Maturity data processing/Data/tanner_survey_specimenEBS.rda"),
                       species = species,
                       years = years, # set in libs/params script
                       sex = "male",
                       shell_condition = "new_hardshell",
                       bin_1mm = TRUE) %>%
  dplyr::select(YEAR, SPECIES, STATION_ID, LATITUDE, LONGITUDE, SIZE_1MM, CPUE, CPUE_MT)

saveRDS(tanner_cpue, "./Maturity data processing/Data/tanner_survey_CPUE_male1mm.rda")


# GET CHELA-DATA FOR SNOW AND TANNER THAT INCLUDES CHELA DB DATA ----
chela_db <- rbind(read.csv(paste0(data_dir, "specimen_chela.csv")), # already != HT 17, only shell 2, no special projects
  read.csv(paste0(data_dir, "specimen_chela_", current.year, ".csv")) %>%
  filter(YEAR == 2025) %>% dplyr::select(!SAMPLING_FACTOR)) %>%
  filter(HAUL_TYPE !=17, SEX == 1, SHELL_CONDITION == 2, is.na(CHELA_HEIGHT) == FALSE) %>% # filter for males, sh2, only chela msrd, not HT17
  mutate(RATIO = SIZE/CHELA_HEIGHT) %>%
  filter(RATIO > 2 & RATIO < 35) %>% # filter extreme measurements
  dplyr::select(!c(RATIO)) %>%
  mutate(LN_CH = log(CHELA_HEIGHT),
         LN_CW = log(SIZE))

write.csv(chela_db, "./Maturity data processing/Data/snow_tanner_cheladatabase.csv")


# Get maturity data----
species <- "SNOW"
param_data <- crabpack::get_male_maturity(species = species,
                                             region = "EBS",
                                             channel = channel)

write.csv(param_data$model_parameters, "./Maturity data processing/Data/snow_legacy_modelpars.csv")


write.csv(param_data$male_mat_ratio, "./Maturity data processing/Data/snow_legacy_matmale_ratio.csv")
legacy <- expand.grid(SPECIES = species,
            SIZE = seq(1, 200, by = 5), # check on this max limit for Tanner!
            YEAR = unique(param_data$model_parameters %>% filter(!is.na(A_EST)) %>% pull(YEAR))) %>%
          left_join(param_data$model_parameters) %>%
          mutate(PROP_MATURE = (1/(1 + exp(-A_EST * (SIZE - B_EST)))))

write.csv(legacy, "./Maturity data processing/Data/snow_legacy_ogives.csv")

species <- "TANNER"
param_data <- crabpack::get_male_maturity(species = species,
                                          region = "EBS",
                                          channel = channel)

write.csv(param_data$model_parameters, "./Maturity data processing/Data/tanner_legacy_modelpars.csv")


write.csv(param_data$male_mat_ratio, "./Maturity data processing/Data/tanner_legacy_matmale_ratio.csv")
legacy <- expand.grid(SPECIES = species,
                      DISTRICT = c("ALL", "W166", "E166"),
                      SIZE = c(1:250), # check on this max limit for Tanner!
                      YEAR = unique(param_data$model_parameters %>% filter(!is.na(A_EST)) %>% pull(YEAR))) %>%
  left_join(param_data$model_parameters) %>%
  mutate(PROP_MATURE = (1/(1 + exp(-A_EST * (SIZE - B_EST)))))

write.csv(legacy, "./Maturity data processing/Data/tanner_legacy_ogives.csv")

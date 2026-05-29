# PURPOSE: 
# To get updated survey data by year for chela maturity processing

# Notes:
# 1) What is the best way to update data? Use the chela database and then specimen data each new year? Or will
# the chela database be updated in time?

# LOAD LIBS/PARAMS -----
source("./Scripts/Sourced scripts/load_libs_params.R")

# Set channel
channel <- "API"

# Set recent year
current_yr <- 2025

# GET SNOW SURVEY SPECIMEN DATA FROM CRABPACK ----
# Pull specimen data
species <- "SNOW"
specimen_data <- crabpack::get_specimen_data(species = species,
                                             region = "EBS",
                                             years = 1989:current_yr,
                                             channel = channel)

saveRDS(specimen_data, "./Data/snow_survey_specimenEBS.rda")


# GET TANNER SURVEY SPECIMEN DATA FROM CRABPACK ----
# Pull specimen data
species <- "TANNER"
specimen_data <- crabpack::get_specimen_data(species = species,
                                             region = "EBS",
                                             years = 1990:current_yr,
                                             channel = channel)

saveRDS(specimen_data, "./Data/tanner_survey_specimenEBS.rda")


# GET CHELA-DATA FOR SNOW AND TANNER THAT INCLUDES CHELA DB DATA ----
chela_db <- rbind(read.csv(paste0(data_dir, "specimen_chela.csv")), # includes special project chela but not >2024; already != HT 17, only shell 2, 
                  read.csv(paste0(data_dir, "specimen_chela_", current_yr, ".csv")) %>% # does not include special project chela but does include >2024
                    filter(YEAR == 2025) %>% dplyr::select(!SAMPLING_FACTOR)) %>% # filtering this database by years >2024 to add to the previous data
  filter(HAUL_TYPE !=17, SEX == 1, SHELL_CONDITION == 2, is.na(CHELA_HEIGHT) == FALSE) %>% # filter for males, sh2, only chela msrd, not HT17
  mutate(RATIO = SIZE/CHELA_HEIGHT) %>%
  filter(RATIO > 2 & RATIO < 35) %>% # filter extreme measurements
  dplyr::select(!c(RATIO)) %>%
  mutate(LN_CH = log(CHELA_HEIGHT),
         LN_CW = log(SIZE))

write.csv(chela_db, "./Data/snow_tanner_cheladatabase.csv")



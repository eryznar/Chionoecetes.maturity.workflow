# PURPOSE: to generate mature bio/abund for sdmTMB maturity models for ESPs

# Author: Emily Ryznar

# NOTES:


# LOAD LIBS/PARAMS ---------------------------------------------------------------------------------------
source("./Maturity data processing/Scripts/load_libs_params.R")

# Source function to simulate from sdmTMB model, and calculate ogives, SAM, and mature bioabund with uncertainty
source("./Maturity data processing/Scripts/calc_maturepop_estimates_function.R")

# Specify function parameters ----
model <- readRDS("./Maturity data processing/Doc/Snow models/sdmTMB_spVAR_noBIN_k300.rda")
crab_data <- readRDS("./Maturity data processing/Data/snow_survey_specimenEBS.rda")
years <- c(1988:2019, 2021:2025)
species <- "SNOW"
output = c("SAM", "bioabund")
size_min = 65
size_max = 80
district = "ALL"
size_1mm = NULL
region = "EBS"

# Run function ----
calc_maturepop_estimates(model, crab_data, years, species, 
                         region, district, size_1mm, size_min, size_max, output) -> snow.out6580


# Calc average SAM
mean_SAM <- mean(na.omit(snow.out$SAM$SAM_mean))
mean_SAM <- mean(na.omit(readRDS("./Maturity data processing/Doc/snow_matpop.rda")$SAM$SAM_mean))

# Specify no chela years
missing <- c(1988, 2008, 2012, 2014, 2016, 2020) 

# Filter spec data to years without chela data (to apply mean SAM as cutline)
crab_data$specimen %>%
   filter(REGION == region, SPECIES == species, 
         SHELL_CONDITION == 2, SEX == 1, YEAR %in% missing) %>% 
  mutate(BIN_5MM = cut_width(SIZE_1MM, width = 5, center = 2.5, closed = "left", dig.lab = 4),
         BIN2 = BIN_5MM) %>%
  separate(BIN2, sep = ",", into = c("LOWER", "UPPER")) %>%
  mutate(LOWER = as.numeric(sub('.', '', LOWER)),
         UPPER = as.numeric(gsub('.$', '', UPPER)),
         SIZE_5MM = (UPPER + LOWER)/2) %>%
  dplyr::select(!c(BIN_5MM, LOWER, UPPER)) %>%
  mutate(YEAR_F = as.factor(YEAR),
         YEAR_SCALED = scale(YEAR),
         CATEGORY = case_when(SIZE >= mean_SAM ~ "Mature male", # use mean SAM to classify as mature/immature
                              TRUE ~ "Immature male")) -> missing_years

# Calc mature and immature bioabund for missing years
missing_years %>% filter(CATEGORY == "Immature male") -> imm_male

missing_imm <- crab_data

missing_imm$specimen <- imm_male

missing_imm_bioabund <-  suppressMessages(crabpack::calc_bioabund(crab_data = missing_imm, species = "SNOW",
                                                                   size_min = 65, size_max = 80, 
                                                                   sex = "male",
                                                                   shell_condition = "new_hardshell")) %>%
  mutate(CATEGORY = "Immature male") %>%
  filter(YEAR %in% missing) %>%
  dplyr::select(YEAR, SPECIES, DISTRICT, CATEGORY, ABUNDANCE, BIOMASS_MT, BIOMASS_LBS)

snow.out6580$mature_bioabund %>% 
  filter(YEAR >=1989, CATEGORY == "Immature male", is.na(ABUNDANCE) == FALSE) %>%
  dplyr::select(YEAR, SPECIES, DISTRICT, CATEGORY, ABUNDANCE, BIOMASS_MT, BIOMASS_LBS) %>%
  mutate(ABUNDANCE = ABUNDANCE * 1e6) -> model_based

rbind(missing_imm_bioabund, model_based) -> imm_male6580

ggplot(imm_male6580, aes(YEAR, ABUNDANCE))+
  geom_line()+
  geom_point()

write.csv(imm_male6580, "./Maturity data processing/Output/snowESP_immale_6580.csv")


# Specify function parameters ----
model <- readRDS("./Maturity data processing/Doc/Snow models/sdmTMB_spVAR_noBIN_k300.rda")
crab_data <- readRDS("./Maturity data processing/Data/snow_survey_specimenEBS.rda")
years <- c(1988:2019, 2021:2025)
species <- "SNOW"
output = c("bioabund")
size_min = 40
size_max = 50
district = "ALL"
size_1mm = NULL
region = "EBS"

# Run function ----
calc_maturepop_estimates(model, crab_data, years, species, 
                         region, district, size_1mm, size_min, size_max, output) -> snow.out4050


# Calculate immature male bioabund for missing chela years based on mean SAM (see above)
missing_imm_bioabund <-  suppressMessages(crabpack::calc_bioabund(crab_data = missing_imm, species = "SNOW",
                                                                  size_min = 40, size_max = 50, 
                                                                  sex = "male",
                                                                  shell_condition = "new_hardshell")) %>%
  mutate(CATEGORY = "Immature male") %>%
  filter(YEAR %in% missing) %>%
  dplyr::select(YEAR, SPECIES, DISTRICT, CATEGORY, ABUNDANCE, BIOMASS_MT, BIOMASS_LBS)

snow.out4050$mature_bioabund %>% 
  filter(YEAR >=1989, CATEGORY == "Immature male", is.na(ABUNDANCE) == FALSE) %>%
  dplyr::select(YEAR, SPECIES, DISTRICT, CATEGORY, ABUNDANCE, BIOMASS_MT, BIOMASS_LBS) %>%
  mutate(ABUNDANCE = ABUNDANCE * 1e6) -> model_based

rbind(missing_imm_bioabund, model_based) -> imm_male4050

ggplot(imm_male4050, aes(YEAR, ABUNDANCE))+
  geom_line()+
  geom_point()

write.csv(imm_male4050, "./Maturity data processing/Output/snowESP_immale_4050.csv")

# Specify function parameters ----
model <- readRDS("./Maturity data processing/Doc/Snow models/sdmTMB_spVAR_noBIN_k300.rda")
crab_data <- readRDS("./Maturity data processing/Data/snow_survey_specimenEBS.rda")
years <- c(1988:2019, 2021:2025)
species <- "SNOW"
output = c("bioabund")
size_min = 50
size_max = 65
district = "ALL"
size_1mm = NULL
region = "EBS"

# Run function ----
calc_maturepop_estimates(model, crab_data, years, species, 
                         region, district, size_1mm, size_min, size_max, output) -> snow.out5065


# Calculate immature male bioabund for missing chela years based on mean SAM (see above)
missing_imm_bioabund <-  suppressMessages(crabpack::calc_bioabund(crab_data = missing_imm, species = "SNOW",
                                                                  size_min = 50, size_max = 65, 
                                                                  sex = "male",
                                                                  shell_condition = "new_hardshell")) %>%
  mutate(CATEGORY = "Immature male") %>%
  filter(YEAR %in% missing) %>%
  dplyr::select(YEAR, SPECIES, DISTRICT, CATEGORY, ABUNDANCE, BIOMASS_MT, BIOMASS_LBS)

snow.out5065$mature_bioabund %>% 
  filter(YEAR >=1989, CATEGORY == "Immature male", is.na(ABUNDANCE) == FALSE) %>%
  dplyr::select(YEAR, SPECIES, DISTRICT, CATEGORY, ABUNDANCE, BIOMASS_MT, BIOMASS_LBS) %>%
  mutate(ABUNDANCE = ABUNDANCE * 1e6) -> model_based

rbind(missing_imm_bioabund, model_based) -> imm_male5065

ggplot(imm_male5065, aes(YEAR, ABUNDANCE))+
  geom_line()+
  geom_point()

write.csv(imm_male5065, "./Maturity data processing/Output/snowESP_immale_5065.csv")


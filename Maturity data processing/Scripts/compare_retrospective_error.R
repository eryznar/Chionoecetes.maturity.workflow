# PURPOSE: to explore fitting sdmTMB models to estimate p(mat) for Chionoecetes

# Author: Emily Ryznar

# NOTES:


# LOAD LIBS/PARAMS ---------------------------------------------------------------------------------------
source("./Maturity data processing/Scripts/load_libs_params.R")

# Source function to simulate from sdmTMB model, and calculate ogives, SAM, and mature bioabund with uncertainty
source("./Maturity data processing/Scripts/calc_maturepop_estimates_function.R")


# LOAD DATA ----------------------------------------------------------------------------------------------
recent_years <- c(2022:2025)
# SNOW ----
# Load data ----
# Load minima data, calculate cutline params
snow_minima <- read.csv("./Maturity data processing/Output/chela_cutline_minima.csv") %>%
  filter(SPECIES == "SNOW") %>%
  mutate(BETA0 = coef(lm(MINIMUM ~ MIDPOINT))[1],
         BETA1 = coef(lm(MINIMUM ~ MIDPOINT))[2])

BETA0 <- unique(snow_minima$BETA0)
BETA1 <- unique(snow_minima$BETA1)

# Load snow specimen data (subsample 1)
snow.specimen<-  readRDS("./Maturity data processing/Data/snow_survey_specimenEBS.rda")$specimen %>%
  filter(SHELL_CONDITION == 2, SEX == 1) %>%
  mutate(SIZE_1MM = floor(SIZE),
         BIN = cut_width(SIZE, width = 5, center = 2.5, closed = "left", dig.lab = 4),
         BIN2 = BIN) %>%
  separate(BIN2, sep = ",", into = c("LOWER", "UPPER")) %>%
  mutate(LOWER = as.numeric(sub('.', '', LOWER)),
         UPPER = as.numeric(gsub('.$', '', UPPER)),
         SIZE_BINNED = (UPPER + LOWER)/2) %>%
  st_as_sf(., coords = c("LONGITUDE", "LATITUDE"), crs = "+proj=longlat +datum=WGS84") %>%
  st_transform(., crs = "+proj=utm +zone=2") %>%
  cbind(st_coordinates(.)) %>%
  as.data.frame(.) %>%
  mutate(LATITUDE = Y/1000, # scale to km so values don't get too large
         LONGITUDE = X/1000,
         SIZE_CATEGORY = as.factor(paste0("SIZE", SIZE_BINNED)),
         YEAR_F = as.factor(YEAR)) %>%
  rename(SIZE_5MM = SIZE_BINNED) %>%
  dplyr::select(YEAR, YEAR_F, STATION_ID, LATITUDE, LONGITUDE, SIZE_5MM, SIZE_CATEGORY, 
                SAMPLING_FACTOR) %>%
  filter(!YEAR %in% c(2012, recent_years))

# Load snow chela data from chela database (subsample 2)
snow.chela <-  read.csv("./Maturity data processing/Data/snow_tanner_cheladatabase.csv") %>% #already filtered appropriately
  dplyr::select(!X) %>%
  filter(SPECIES == "SNOW", SIZE >= 35 & SIZE <= 135) %>%
  mutate(CUTOFF = BETA0 + BETA1*LN_CW, # apply cutline model
         MATURE = case_when((LN_CH > CUTOFF) ~ 1,
                            TRUE ~ 0),
         SIZE_1MM = floor(SIZE),
         BIN = cut_width(SIZE, width = 5, center = 2.5, closed = "left", dig.lab = 4),
         BIN2 = BIN) %>%
  separate(BIN2, sep = ",", into = c("LOWER", "UPPER")) %>%
  mutate(LOWER = as.numeric(sub('.', '', LOWER)),
         UPPER = as.numeric(gsub('.$', '', UPPER)),
         SIZE_BINNED = (UPPER + LOWER)/2) %>%
  # MATURE = case_when((SIZE <=35) ~ 0,
  #                    (SIZE >= 135) ~ 1,
  #TRUE ~ MATURE)) %>%
  st_as_sf(., coords = c("LONGITUDE", "LATITUDE"), crs = "+proj=longlat +datum=WGS84") %>%
  st_transform(., crs = "+proj=utm +zone=2") %>%
  cbind(st_coordinates(.)) %>%
  as.data.frame(.) %>%
  #filter(N_chela > 40) %>% #necessary for consecutive size correlations to run
  mutate(LATITUDE = Y/1000, # scale to km so values don't get too large
         LONGITUDE = X/1000,
         SIZE_CATEGORY = as.factor(paste0("SIZE", SIZE_BINNED)),
         YEAR_F = as.factor(YEAR),
         YEAR_SCALED = as.numeric(scale(YEAR)),
         MATURE = case_when((SIZE <=35) ~ 0,
                            (SIZE >= 135) ~ 1,
                            TRUE ~ MATURE)) %>%
  as.data.frame(.) %>%
  rename(SIZE_5MM = SIZE_BINNED) %>%
  dplyr::select(YEAR, YEAR_F, YEAR_SCALED, STATION_ID, LATITUDE, LONGITUDE, SIZE, SIZE_5MM, SIZE_CATEGORY, MATURE) %>%
  filter(!YEAR %in% c(2012, recent_years))

# Fit models ----
mod.dat <- snow.chela

set.seed(1)

# Set params
mat.msh <- sdmTMB::make_mesh(mod.dat, c("LONGITUDE","LATITUDE"), n_knots = 300, type = "kmeans")

xtra.time <- c(2008, 2012, 2014, 2016, 2020) # missing years across all size bins

mod.2 <- sdmTMB(
  MATURE ~ s(SIZE, k = 13) + YEAR_SCALED, 
  spatial = "on",
  spatiotemporal = "iid",
  mesh = mat.msh,
  family = binomial(),
  time = "YEAR",
  spatial_varying = ~ 0 + SIZE,
  extra_time = xtra.time,
  anisotropy = TRUE,
  data = mod.dat
)

saveRDS(mod.2, "./Maturity data processing/Doc/Snow models/sdmTMB_spVAR_noBIN_k300_retrospective.rda")




# Calc pop estimates ----
model <- readRDS("./Maturity data processing/Doc/Snow models/sdmTMB_spVAR_noBIN_k300_retrospective.rda")
snow_dat <- readRDS("./Maturity data processing/Data/snow_survey_specimenEBS.rda")
snow_yrs <- c(1989:2019, 2021) # omitting recent years
crab_data = snow_dat
years = snow_yrs
species = "SNOW"
output = NULL


s.out <- calc_maturepop_estimates(model, crab_data, years, species, output)
saveRDS(s.out, "./Maturity data processing/Doc/snow_matpop_retrospective.rda")


# TANNER ----
# Load data ----
# Load minima data, calculate cutline params
tanner_minima <- read.csv("./Maturity data processing/Output/chela_cutline_minima.csv") %>%
  filter(SPECIES == "TANNER") %>%
  mutate(BETA0 = coef(lm(MINIMUM ~ MIDPOINT))[1],
         BETA1 = coef(lm(MINIMUM ~ MIDPOINT))[2])

BETA0 <- unique(tanner_minima$BETA0)
BETA1 <- unique(tanner_minima$BETA1)

# Load tanner specimen data (subsample 1)
tanner.specimen<-  readRDS("./Maturity data processing/Data/tanner_survey_specimenEBS.rda")$specimen %>%
  filter(SHELL_CONDITION == 2, SEX == 1) %>%
  mutate(SIZE_1MM = floor(SIZE),
         BIN = cut_width(SIZE, width = 5, center = 2.5, closed = "left", dig.lab = 4),
         BIN2 = BIN) %>%
  separate(BIN2, sep = ",", into = c("LOWER", "UPPER")) %>%
  mutate(LOWER = as.numeric(sub('.', '', LOWER)),
         UPPER = as.numeric(gsub('.$', '', UPPER)),
         SIZE_BINNED = (UPPER + LOWER)/2) %>%
  st_as_sf(., coords = c("LONGITUDE", "LATITUDE"), crs = "+proj=longlat +datum=WGS84") %>%
  st_transform(., crs = "+proj=utm +zone=2") %>%
  cbind(st_coordinates(.)) %>%
  as.data.frame(.) %>%
  mutate(LATITUDE = Y/1000, # scale to km so values don't get too large
         LONGITUDE = X/1000,
         SIZE_CATEGORY = as.factor(paste0("SIZE", SIZE_BINNED)),
         YEAR_F = as.factor(YEAR),
         YEAR_SCALED = scale(YEAR)) %>%
  rename(SIZE_5MM = SIZE_BINNED) %>%
  dplyr::select(YEAR, YEAR_F, YEAR_SCALED, DISTRICT, STATION_ID, LATITUDE, LONGITUDE, SIZE, SIZE_5MM, SIZE_CATEGORY, 
                SAMPLING_FACTOR) %>%
  filter(!YEAR %in% recent_years)

# Load tanner chela data from chela database (subsample 2)

tanner.chela <-  read.csv("./Maturity data processing/Data/snow_tanner_cheladatabase.csv") %>% #already filtered appropriately
  dplyr::select(!X) %>%
  filter(SPECIES == "TANNER", SIZE >= 55 & SIZE <= 145) %>% # filtering ! sizes without separation
  mutate(CUTOFF = BETA0 + BETA1*LN_CW, # apply cutline model
         MATURE = case_when((LN_CH > CUTOFF) ~ 1,
                            TRUE ~ 0),
         SIZE_1MM = floor(SIZE),
         BIN = cut_width(SIZE, width = 5, center = 2.5, closed = "left", dig.lab = 4),
         BIN2 = BIN) %>%
  separate(BIN2, sep = ",", into = c("LOWER", "UPPER")) %>%
  mutate(LOWER = as.numeric(sub('.', '', LOWER)),
         UPPER = as.numeric(gsub('.$', '', UPPER)),
         SIZE_BINNED = (UPPER + LOWER)/2) %>%
  st_as_sf(., coords = c("LONGITUDE", "LATITUDE"), crs = "+proj=longlat +datum=WGS84") %>%
  st_transform(., crs = "+proj=utm +zone=2") %>%
  cbind(st_coordinates(.)) %>%
  as.data.frame(.) %>%
  mutate(LATITUDE = Y/1000, # scale to km so values don't get too large
         LONGITUDE = X/1000,
         SIZE_CATEGORY = as.factor(paste0("SIZE", SIZE_BINNED)),
         YEAR_F = as.factor(YEAR),
         YEAR_SCALED = as.numeric(scale(YEAR)),
         MATURE = case_when((SIZE <=55) ~ 0,
                            (SIZE >=145) ~ 1, # or just apply differently to DISTRICTs, though that made the mature plots wonky
                            TRUE ~ MATURE)) %>%
  as.data.frame(.) %>%
  rename(SIZE_5MM = SIZE_BINNED) %>%
  dplyr::select(YEAR, DISTRICT, YEAR_F, YEAR_SCALED, STATION_ID, LATITUDE, LONGITUDE, SIZE, SIZE_5MM, SIZE_CATEGORY, MATURE) %>%
  filter(!YEAR %in% recent_years)

# Fit models ----
set.seed(1)
mod.dat <- tanner.chela
# Set params
mat.msh <- sdmTMB::make_mesh(mod.dat, c("LONGITUDE","LATITUDE"), n_knots = 200, type = "kmeans")
# Extra time
xtra.time <- c(2013, 2015, 2020) # missing years across all size bins

mod.2 <- sdmTMB(
  MATURE ~ s(SIZE, k = 10) + YEAR_SCALED, 
  spatial = "on",
  spatiotemporal = "iid",
  mesh = mat.msh,
  family = binomial(),
  time = "YEAR",
  spatial_varying = ~ 0 + SIZE,
  extra_time = xtra.time,
  anisotropy = TRUE,
  data = mod.dat
)

saveRDS(mod.2, "./Maturity data processing/Doc/Tanner models/sdmTMB_spVAR_noBIN_k200_retrospective.rda")

# Calc pop estimates ----
model <- readRDS("./Maturity data processing/Doc/Tanner models/sdmTMB_spVAR_noBIN_k200_retrospective.rda")
tanner_dat <- readRDS("./Maturity data processing/Data/tanner_survey_specimenEBS.rda")
tanner_yrs <- c(1990:2019, 2021) # omitting recent years
crab_data = tanner_dat
years = tanner_yrs
species = "TANNER"
output = NULL


t.out <- calc_maturepop_estimates(model, crab_data, years, species, output)
 
saveRDS(t.out, "./Maturity data processing/Doc/tanner_matpop_retrospective.rda")

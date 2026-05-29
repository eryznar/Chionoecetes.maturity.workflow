# PURPOSE: to fit sdmTMB models to estimate p(mat) for snow crab

# Author: Emily Ryznar

# NOTES:


# LOAD LIBS/PARAMS ---------------------------------------------------------------------------------------
source("./Scripts/Sourced scripts/load_libs_params.R")

# LOAD DATA ----------------------------------------------------------------------------------------------
# Load minima data, calculate cutline params
snow_minima <- read.csv("./Output/chela_cutline_minima.csv") %>%
  filter(SPECIES == "SNOW") %>%
  mutate(BETA0 = coef(lm(MINIMUM ~ MIDPOINT))[1],
         BETA1 = coef(lm(MINIMUM ~ MIDPOINT))[2])

BETA0 <- unique(snow_minima$BETA0)
BETA1 <- unique(snow_minima$BETA1)

# Load snow specimen data (subsample 1) and filter appropriately
snow.specimen<-  readRDS("./Data/snow_survey_specimenEBS.rda")$specimen %>%
  filter(SHELL_CONDITION == 2, SEX == 1) %>%
  st_as_sf(., coords = c("LONGITUDE", "LATITUDE"), crs = "+proj=longlat +datum=WGS84") %>%
  st_transform(., crs = "+proj=utm +zone=2") %>%
  cbind(st_coordinates(.)) %>%
  as.data.frame(.) %>%
  mutate(LATITUDE = Y/1000, # scale to km so values don't get too large
         LONGITUDE = X/1000,
         YEAR_F = as.factor(YEAR)) %>%
  dplyr::select(YEAR, YEAR_F, SIZE, STATION_ID, LATITUDE, LONGITUDE, 
                SAMPLING_FACTOR) %>%
  filter(YEAR != 2012) # very little data in this year

# Load snow chela data from chela database (subsample 2)
snow.chela <-  read.csv("./Data/snow_tanner_cheladatabase.csv") %>% #already filtered appropriately
  dplyr::select(!X) %>%
  filter(SPECIES == "SNOW", SIZE >= 35 & SIZE <= 135) %>% # filter to sizes that have variation in maturity
  mutate(CUTOFF = BETA0 + BETA1*LN_CW, # apply cutline model
         MATURE = case_when((LN_CH > CUTOFF) ~ 1,
                            TRUE ~ 0)) %>%
  st_as_sf(., coords = c("LONGITUDE", "LATITUDE"), crs = "+proj=longlat +datum=WGS84") %>% # make spatial
  st_transform(., crs = "+proj=utm +zone=2") %>% # transform to utm
  cbind(st_coordinates(.)) %>%
  as.data.frame(.) %>%
  mutate(LATITUDE = Y/1000, # scale to km so values don't get too large
         LONGITUDE = X/1000,
         YEAR_F = as.factor(YEAR),
         YEAR_SCALED = as.numeric(scale(YEAR)),
         MATURE = case_when((SIZE <=35) ~ 0,
                            (SIZE >= 135) ~ 1,
                            TRUE ~ MATURE)) %>%
  as.data.frame(.) %>%
  dplyr::select(YEAR, YEAR_F, YEAR_SCALED, STATION_ID, LATITUDE, LONGITUDE, SIZE, MATURE) %>%
  filter(YEAR != 2012) # very little data in this year


# FIT MATURITY SDMTMB ----
# Specify model data
mod.dat <- snow.chela

# Make mesh
mat.msh <- sdmTMB::make_mesh(mod.dat, c("LONGITUDE","LATITUDE"), n_knots = 300, type = "kmeans")

# Specify xtra time
xtra.time <- c(2008, 2012, 2014, 2016, 2020) # missing years across all size bins

# Fit model (model parameters have already been vetted as the "best")
snow.mod <- sdmTMB(
  MATURE ~ s(SIZE, k = 13) + YEAR_SCALED, 
  spatial = "on",
  spatiotemporal = "iid",
  mesh = mat.msh,
  family = binomial(),
  time = "YEAR",
  spatial_varying = ~ 0 + SIZE,
  extra_time = xtra.time,
  anisotropy = TRUE,
  data = mod.dat,
  silent = FALSE
)

saveRDS(snow.mod, "./Models/snow_maturity_sdmTMB.rda")


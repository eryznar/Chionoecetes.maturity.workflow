# PURPOSE: to explore fitting sdmTMB models to estimate p(mat) for Chionoecetes

# Author: Emily Ryznar

# NOTES:


# LOAD LIBS/PARAMS -----
source("./Scripts/Sourced scripts/load_libs_params.R")

# LOAD DATA -----
# Load minima data, calculate cutline params
tanner_minima <- read.csv("./Output/chela_cutline_minima.csv") %>%
  filter(SPECIES == "TANNER") %>%
  mutate(BETA0 = coef(lm(MINIMUM ~ MIDPOINT))[1],
         BETA1 = coef(lm(MINIMUM ~ MIDPOINT))[2])

BETA0 <- unique(tanner_minima$BETA0)
BETA1 <- unique(tanner_minima$BETA1)

# Load tanner specimen data (subsample 1)
tanner.specimen<-  readRDS("./Data/tanner_survey_specimenEBS.rda")$specimen %>%
  filter(SHELL_CONDITION == 2, SEX == 1) %>%
  st_as_sf(., coords = c("LONGITUDE", "LATITUDE"), crs = "+proj=longlat +datum=WGS84") %>%
  st_transform(., crs = "+proj=utm +zone=2") %>%
  cbind(st_coordinates(.)) %>%
  as.data.frame(.) %>%
  mutate(LATITUDE = Y/1000, # scale to km so values don't get too large
         LONGITUDE = X/1000,
         YEAR_F = as.factor(YEAR),
         YEAR_SCALED = scale(YEAR)) %>%
  dplyr::select(YEAR, YEAR_F, YEAR_SCALED, DISTRICT, STATION_ID, LATITUDE, LONGITUDE, 
                SAMPLING_FACTOR) 

# Load tanner chela data from chela database (subsample 2)
tanner.chela <-  read.csv("./Data/snow_tanner_cheladatabase.csv") %>% #already filtered appropriately
  dplyr::select(!X) %>%
  filter(SPECIES == "TANNER", SIZE >= 55 & SIZE <= 145) %>% # filtering ! sizes without separation
  mutate(CUTOFF = BETA0 + BETA1*LN_CW, # apply cutline model
         MATURE = case_when((LN_CH > CUTOFF) ~ 1,
                            TRUE ~ 0)) %>%
  st_as_sf(., coords = c("LONGITUDE", "LATITUDE"), crs = "+proj=longlat +datum=WGS84") %>%
  st_transform(., crs = "+proj=utm +zone=2") %>% # transform to UTM
  cbind(st_coordinates(.)) %>%
  as.data.frame(.) %>%
  mutate(LATITUDE = Y/1000, # scale to km so values don't get too large
         LONGITUDE = X/1000,
         YEAR_F = as.factor(YEAR),
         YEAR_SCALED = as.numeric(scale(YEAR)),
         MATURE = case_when((SIZE <=55) ~ 0,
                            (SIZE >=145) ~ 1, # or just apply differently to DISTRICTs, though that made the mature plots wonky
                            TRUE ~ MATURE)) %>%
  as.data.frame(.) %>%
  dplyr::select(YEAR, DISTRICT, YEAR_F, YEAR_SCALED, STATION_ID, LATITUDE, LONGITUDE, SIZE,  MATURE) 

# FIT MATURITY SDMTMB ----
# Make mesh
mat.msh <- sdmTMB::make_mesh(tanner.chela, c("LONGITUDE","LATITUDE"), n_knots = 200, type = "kmeans")

# Specify extra time
xtra.time <- c(2013, 2015, 2020) # missing years across all size bins

# Fit model (model parameters have already been vetted as the "best")
tanner.mod <- sdmTMB(
  MATURE ~ s(SIZE, k = 10) + YEAR_SCALED, 
  spatial = "on",
  spatiotemporal = "iid",
  mesh = mat.msh,
  family = binomial(),
  time = "YEAR",
  spatial_varying = ~ 0 + SIZE,
  extra_time = xtra.time,
  anisotropy = TRUE,
  data = tanner.chela,
  silent = FALSE
)

saveRDS(tanner.mod, "./Models/tanner_maturity_sdmTMB.rda")

# PURPOSE: to evaluate whether maturity exhibits spatiotemporal variation

# Author: Emily Ryznar

# NOTES:


# LOAD LIBS/PARAMS ---------------------------------------------------------------------------------------
source("./Maturity data processing/Scripts/1) load_libs_params.R")

# LOAD DATA ----------------------------------------------------------------------------------------------
# Load minima data, calculate cutline params
snow_minima <- read.csv("./Maturity data processing/Output/chela_cutline_minima.csv") %>%
  filter(SPECIES == "SNOW") %>%
  mutate(BETA0 = coef(lm(MINIMUM ~ MIDPOINT))[1],
         BETA1 = coef(lm(MINIMUM ~ MIDPOINT))[2])

BETA0 <- unique(snow_minima$BETA0)
BETA1 <- unique(snow_minima$BETA1)

# Load snow specimen data (subsample 1)
snow_specimen<-  readRDS("./Maturity data processing/Data/snow_survey_specimenEBS.rda")$specimen %>%
  filter(SHELL_CONDITION == 2, SEX == 1) %>%
  mutate(SIZE_1MM = floor(SIZE),
         BIN_5MM = cut_width(SIZE_1MM, width = 5, center = 2.5, closed = "left", dig.lab = 4),
         BIN2 = BIN_5MM) %>%
  separate(BIN2, sep = ",", into = c("LOWER", "UPPER")) %>%
  mutate(LOWER = as.numeric(sub('.', '', LOWER)),
         UPPER = as.numeric(gsub('.$', '', UPPER)),
         SIZE_5MM = (UPPER + LOWER)/2) %>%
  dplyr::select(YEAR, STATION_ID, LATITUDE, LONGITUDE, SIZE_1MM, SIZE_5MM, 
                SAMPLING_FACTOR, AREA_SWEPT, CALCULATED_WEIGHT_1MM)

# Load snow chela data from chela database (subsample 2)
snow_chela <-  read.csv("./Maturity data processing/Data/snow_tanner_cheladatabase.csv") %>% #already filtered appropriately
  filter(SPECIES == "SNOW") %>%
  mutate(CUTOFF = BETA0 + BETA1*LN_CW, # apply cutline model
         MATURE = case_when((LN_CH > CUTOFF) ~ 1,
                            TRUE ~ 0),
         SIZE_1MM = floor(SIZE),
         BIN_5MM = cut_width(SIZE_1MM, width = 5, center = 2.5, closed = "left", dig.lab = 4),
         BIN2 = BIN_5MM) %>%
  separate(BIN2, sep = ",", into = c("LOWER", "UPPER")) %>%
  mutate(LOWER = as.numeric(sub('.', '', LOWER)),
         UPPER = as.numeric(gsub('.$', '', UPPER)),
         SIZE_5MM = (UPPER + LOWER)/2,
         MATURE = case_when((SIZE <40) ~ 0,
                            (SIZE >= 135) ~ 1,
                            TRUE ~ MATURE)) %>%
  dplyr::select(YEAR, STATION_ID, LATITUDE, LONGITUDE, SIZE_1MM, SIZE_5MM, 
                LN_CH, LN_CW, MATURE, AREA_SWEPT)

# FIT GAMS -----------------------------------------------------------------------------------------------
mod.1 <- gam(MATURE ~ s(SIZE_1MM, k = 4) + s(YEAR, k = 4) + s(LATITUDE, LONGITUDE, k = 4), data = snow_chela, family = "binomial")
summary(mod.1)

mod.2 <- gam(MATURE ~ s(SIZE_1MM, k = 4) + s(YEAR, k = 4), data = snow_chela, family = "binomial")

mod.3 <- gam(MATURE ~ s(SIZE_1MM, k = 4) + + s(LATITUDE, LONGITUDE, k = 4), data = snow_chela, family = "binomial")

mod.4 <- gam(MATURE ~ s(SIZE_1MM, k = 4), data = snow_chela, family = "binomial")


AICc(mod.1, mod.2, mod.3, mod.4) %>%
  cbind(., data.frame(ST_pars = c("s(YEAR) + s(LAT*LON)", "s(YEAR)", "s(LAT*LON)", "None")))

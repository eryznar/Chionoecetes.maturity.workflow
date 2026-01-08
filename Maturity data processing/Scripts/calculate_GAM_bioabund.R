# PURPOSE: to generate "model-based" indices of mature abundance and biomass using a GAMs

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
  filter(YEAR != 2012)

# Load snow chela data (subsample2)
snow_chela <-  read.csv("./Maturity data processing/Data/snow_tanner_cheladatabase.csv") %>% #already filtered appropriately
  dplyr::select(!X) %>%
  filter(SPECIES == "SNOW", SIZE >= 35 & SIZE <= 135) %>%
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
         MATURE = case_when((SIZE <=35) ~ 0,
                            (SIZE >= 135) ~ 1,
                            TRUE ~ MATURE)) %>%
  filter(YEAR != 2012)

# FIT GAM MODELS ------------------------------------------------------------------
# Fit mods (no spatial) ----
gam.mods <- list()
yrs <- unique(snow_chela$YEAR)
aic <- data.frame()

for(ii in 1:length(yrs)){
  print(paste("Fitting ", yrs[ii]))
  
  snow_chela %>% filter(YEAR == yrs[ii]) -> mod.dat
  
  gam.mod <- gam(MATURE ~ s(SIZE_5MM, k = 4), data = mod.dat, family = "binomial")
 
  gam.mods[[as.character(yrs[ii])]] <- gam.mod

  aic <- rbind(aic, data.frame(year = yrs[ii], aic = AICc(gam.mod)))
}

saveRDS(gam.mods, "./Maturity research/Output/gam_mods.rda")

# Predict mods
gam.preds <- data.frame()

for(ii in 1:length(yrs)){
  print(paste("Predicting ", yrs[ii]))
  
  # select year model
  mod <- gam.mods[[as.character(yrs[ii])]]
  
  # filter specimen data by year
  snow_specimen %>% filter(YEAR == yrs[ii]) -> sub1
  
  # predict model
  pp <- cbind(sub1, predict(mod, sub1, , type = "response", se = TRUE))
  
  # bind predictions across year
  gam.preds <- rbind(gam.preds, pp) 
  
}

# Multiply sampling factor by prop mature in each size bin
gam.cpue<- gam.preds %>%
  rename(PROP_MATURE = fit) %>%
  mutate(SAMPLING_FACTOR = SAMPLING_FACTOR * PROP_MATURE)

# Bind prediction dataframe to crab.dat to calculate bioabund using crabpack
crab.dat <-  readRDS("./Maturity data processing/Data/snow_survey_specimenEBS.rda")

crab.dat$specimen <- gam.cpue

# Calculate mature biomass and abundance using
bioabund <- crabpack::calc_bioabund(crab_data = crab.dat, species = "SNOW", 
                                       sex = "male", shell_condition = "new hardshell") %>%
  mutate(ABUNDANCE = ABUNDANCE/1e6,
         ABUNDANCE_CI = ABUNDANCE_CI/1e6) %>%
  right_join(., data.frame(YEAR = seq(min(.$YEAR), max(.$YEAR)))) %>%
  rename(BIOMASS = BIOMASS_MT,
         BIOMASS_CI = BIOMASS_MT_CI)

# write.csv(bioabund, "./Maturity research/Output/GAM_mature_bioabund.csv")

# Fit mods (spatial) ----
gam.mods.ST <- list()
yrs <- unique(snow_chela$YEAR)
aic.ST <- data.frame()
for(ii in 1:length(yrs)){
  print(paste("Fitting ", yrs[ii]))
  
  snow_chela %>% filter(YEAR == yrs[ii]) -> mod.dat
  
  gam.mod <- gam(MATURE ~ s(SIZE_5MM, k = 4) + s(LATITUDE, LONGITUDE, k = 4), data = mod.dat, family = "binomial")
  
  gam.mods.ST[[as.character(yrs[ii])]] <- gam.mod
  
  aic.ST <- rbind(aic.ST, data.frame(year = yrs[ii], aic = AICc(gam.mod)))
  
}

saveRDS(gam.mods.ST, "./Maturity research/Output/gam_mods_spatial.rda")

# Fit GAMM (via Mullowny and Baker 2021) ----
m1 <- gam(MATURE ~ s(SIZE_5MM, k = 4) + s(YEAR, k = 4) + ti(SIZE_5MM, YEAR) + s(LATITUDE, LONGITUDE), 
                data = snow_chela, family = "binomial", method = "ML")

m2 <- gam(MATURE ~ s(SIZE_5MM, k = 4), 
          data = snow_chela, family = "binomial", method = "ML")

m3 <- gam(MATURE ~ s(SIZE_5MM, k = 4) + s(YEAR, k = 4), 
               data = snow_chela, family = "binomial", method = "ML")

m4 <- gam(MATURE ~ s(SIZE_5MM, k = 4) + s(YEAR, k = 4) + ti(SIZE_5MM, YEAR), 
               data = snow_chela, family = "binomial", method = "ML")

m5 <- gam(MATURE ~ s(SIZE_5MM, k = 4) + s(LATITUDE, LONGITUDE), 
               data = snow_chela, family = "binomial", method = "ML")

m6 <-  gam(MATURE ~ s(SIZE_5MM, k = 4) + s(YEAR, k = 4) + s(LATITUDE, LONGITUDE), 
                data = snow_chela, family = "binomial", method = "ML")

m7 <-  gamm(MATURE ~ s(SIZE_5MM, k = 4) + s(YEAR, k = 4) + s(LATITUDE, LONGITUDE) + s(YEAR, bs = "re"), 
           data = snow_chela, family = "binomial", method = "ML")


cbind(AICc(m1, m2, m3, m4, m5, m6, mod.1, mod.2), terms = c("GAM: s(SIZE_5MM) + s(YEAR) + ti(SIZE_5MM, YEAR) + s(LATITUDE, LONGITUDE)",
                                          "GAM: s(SIZE_5MM)",
                                          "GAM: s(SIZE_5MM) + s(YEAR)",
                                          "GAM: s(SIZE_5MM) + s(YEAR) + ti(SIZE_5MM, YEAR)",
                                          "GAM: s(SIZE_5MM) + s(LATITUDE, LONGITUDE)",
                                          "GAM: s(SIZE_5MM) + s(YEAR) + s(LATITUDE, LONGITUDE)",
                                          "sdmTMB: s(SIZE_5MM) + s(YEAR)",
                                          "sdmTMB: s(SIZE_5MM) + YEAR_F")) %>%
  arrange(., AICc) 
                                          
      

saveRDS(m1, "./Maturity research/Output/gam_pooled.rda")


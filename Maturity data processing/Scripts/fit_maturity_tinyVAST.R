# PURPOSE: to explore fitting tinyVAST models to estimate p(mat) by size for Chionoecetes

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
                    st_as_sf(., coords = c("LONGITUDE", "LATITUDE"), crs = "+proj=longlat +datum=WGS84") %>%
                    st_transform(., crs = "+proj=utm +zone=2") %>%
                    cbind(st_coordinates(.)) %>%
                    as.data.frame(.) %>%
                    mutate(LATITUDE = Y/1000, # scale to km so values don't get too large
                           LONGITUDE = X/1000) %>%
                    dplyr::select(YEAR, STATION_ID, LATITUDE, LONGITUDE, SIZE_1MM, SIZE_5MM, 
                                  SAMPLING_FACTOR, CALCULATED_WEIGHT_1MM)

# Load snow chela data from chela database (subsample 2)
snow_chela <-  read.csv("./Maturity data processing/Data/snow_tanner_cheladatabase.csv") %>% #already filtered appropriately
                  dplyr::select(!X) %>%
                  filter(SPECIES == "SNOW") %>%
                  mutate(CUTOFF = BETA0 + BETA1*LN_CW, # apply cutline model
                         MATURE = case_when((LN_CH > CUTOFF) ~ 1,
                                            TRUE ~ 0),
                         SIZE_1MM = floor(SIZE),
                         BIN_10MM = cut_width(SIZE_1MM, width = 10, center = 5, closed = "left", dig.lab = 4),
                         BIN2 = BIN_10MM) %>%
                  separate(BIN2, sep = ",", into = c("LOWER", "UPPER")) %>%
                  mutate(LOWER = as.numeric(sub('.', '', LOWER)),
                         UPPER = as.numeric(gsub('.$', '', UPPER)),
                         SIZE_10MM = (UPPER + LOWER)/2,
                         MATURE = case_when((SIZE <=35) ~ 0,
                                            (SIZE >= 135) ~ 1,
                                            TRUE ~ MATURE)) %>%
                  st_as_sf(., coords = c("LONGITUDE", "LATITUDE"), crs = "+proj=longlat +datum=WGS84") %>%
                  st_transform(., crs = "+proj=utm +zone=2") %>%
                  cbind(st_coordinates(.)) %>%
                  as.data.frame(.) %>%
                  filter(SIZE_1MM >= 35 & SIZE_1MM <= 135) %>%
                  #filter(N_chela > 40) %>% #necessary for consecutive size correlations to run
                  mutate(LATITUDE = Y/1000, # scale to km so values don't get too large
                         LONGITUDE = X/1000,
                         SIZE_CATEGORY = as.factor(paste0("SIZE_", SIZE_10MM))) %>%
                  as.data.frame(.) %>%
                  dplyr::select(YEAR, STATION_ID, LATITUDE, LONGITUDE, SIZE_10MM, SIZE_CATEGORY, MATURE) %>%
                  filter(YEAR != 2012)
  

#write.csv(snow_chela, "./Maturity research/Output/tinyVAST_snow_data.csv")

# FIT tinyVAST MODELS ------------------------------------------------------------------------------------
# 1) specify mesh
mat.mesh = fmesher::fm_mesh_2d(snow_chela[, c("LONGITUDE", "LATITUDE")], cutoff = 60) # higher cutoff, lower # knots

nrow(mat.mesh$loc) # N knots

# 2) specify spacetime term
sizes <- sort(unique(snow_chela$SIZE_10MM))
same_size_cor <- c() # correlations over time within age
next_size_cor <- c() # correlations among ages within a year

for(ii in seq_along(sizes)) {
  same_size_cor[ii] <- paste0("SIZE_", sizes[ii], " -> SIZE_", sizes[ii], ", 1, rho_T")
}

for(ii in 1:(length(sizes)-1)) {
  next_size_cor[ii] <- paste0("SIZE_", sizes[ii], " -> SIZE_", sizes[ii+1], ", 0, rho_S")
}

# Combine vectors and join with single newlines, then trim
spacetime_term <- paste(c(same_size_cor, next_size_cor), collapse = "\n")
spacetime_term <- trimws(spacetime_term)
cat(spacetime_term)

# Look at sample sizes
snow_chela %>%
  count(SIZE_10MM, YEAR)


# specify tinyVAST controls
control <- tinyVASTcontrol(profile = c("alpha_j","alpha2_j"),
                           trace = 1,
                           #nlminb_loops = 0,
                           #newton_loops = 0,
                           #getsd = FALSE,
                           calculate_deviance_explained = TRUE)

# Fit model
mod.1 = tinyVAST(
          data = snow_chela,
          formula = MATURE ~ s(SIZE_10MM, k = 3) + s(YEAR, k = 3) + ti(SIZE_10MM, YEAR, k = c(3,3)),
          spacetime_term = spacetime_term,
          family = binomial(),
          space_column = c("LATITUDE", "LONGITUDE"), 
          variable_column = "SIZE_CATEGORY",
          time_column = "YEAR",
          spatial_domain = mat.mesh,
          control = control) 

# PDE = 0.2316

saveRDS(mod.1, "./Maturity research/Output/tinyVAST_fullsize_correlation_no2012.rda")

mod.2 = tinyVAST(
  data = snow_chela,
  formula = MATURE ~ s(SIZE_10MM, k = 3) + s(YEAR, k = 3),
  spacetime_term = spacetime_term,
  family = binomial(),
  space_column = c("LATITUDE", "LONGITUDE"), 
  variable_column = "SIZE_CATEGORY",
  time_column = "YEAR",
  spatial_domain = mat.mesh) 

# PDE = 0.2316

saveRDS(mod.1, "./Maturity research/Output/tinyVAST_fullsize_correlation.rda")

        

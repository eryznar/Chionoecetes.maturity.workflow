# PURPOSE: to explore fitting sdmTMB models to estimate p(mat) for Chionoecetes

# Author: Emily Ryznar

# NOTES:


# LOAD LIBS/PARAMS ---------------------------------------------------------------------------------------
source("./Maturity data processing/Scripts/load_libs_params.R")

# LOAD DATA ----------------------------------------------------------------------------------------------
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
  dplyr::select(YEAR, YEAR_F, YEAR_SCALED, DISTRICT, STATION_ID, LATITUDE, LONGITUDE, SIZE_5MM, SIZE_CATEGORY, 
                SAMPLING_FACTOR) 

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
  dplyr::select(YEAR, DISTRICT, YEAR_F, YEAR_SCALED, STATION_ID, LATITUDE, LONGITUDE, SIZE_5MM, SIZE_CATEGORY, MATURE) 

tanner.chela %>% dplyr::select(YEAR, DISTRICT, LATITUDE, LONGITUDE) %>% distinct() -> rr
ggplot(rr, aes(LONGITUDE, LATITUDE, color = DISTRICT))+
  geom_point(alpha = 0.5)+
  facet_wrap(~YEAR)+
  theme_bw()

# Fit models with year fixed effect and smooth for size ----
mod.dat <- tanner.chela

# Cross validation for models with different k-smooths ----
mod.dat <- tanner.chela

# Make mesh
mat.msh <- sdmTMB::make_mesh(mod.dat, c("LONGITUDE","LATITUDE"), n_knots = 200, type = "kmeans")

# Extra time
xtra.time <- c(2013, 2015, 2020) # missing years across all size bins

kk<- 4:15

# Define folds
set.seed(1)
n_folds <- 5
mod.dat$fold <- sample.int(n_folds, nrow(mod.dat), replace = TRUE)

# Set params
mods <- list()
results <- data.frame()
for(ii in 1:length(kk)){
  
  print(paste0("Fitting k=", kk[ii]))
  mod <- sdmTMB(MATURE ~ s(SIZE_5MM, k = kk[ii]) + YEAR_F, #the 0 is there's a factor predictor for each YEAR, no intercept
                spatial = "on",
                spatiotemporal = "iid",
                mesh = mat.msh,
                family = binomial(),
                time = "YEAR",
                extra_time = xtra.time,
                anisotropy = TRUE,
                data = mod.dat)
  
  mods[[ii]] <- mod
  
  names(mods)[ii] <- paste0("k=", kk[ii])
  
  print(paste0("CVing k=", kk[ii]))
  cv <- sdmTMB_cv(
    MATURE ~ s(SIZE_5MM, k = kk[ii]) + YEAR_F,
    data = mod.dat,
    mesh = mat.msh,
    family = binomial(),
    spatial = "on",
    spatiotemporal = "iid",
    time = "YEAR",
    extra_time = xtra.time,
    anisotropy = TRUE,
    k_folds = n_folds,
    fold_ids = mod.dat$fold
  )
  
  ll <- sum(cv$fold_loglik)
  
  results <- rbind(results, data.frame(k = kk[ii], AICc(mod), logLik = ll))
  
}

results %>%
  arrange(., -logLik) -> CV.results

write.csv(CV.results, "./Maturity data processing/Output/sdmTMB_ksmooth_CV_TANNER.csv")


# Cross validation with different model parameterizations ----
xtra.time <- c(2013, 2015, 2020) # missing years across all size bins
set.seed(1)
n_folds <- 5
tanner.chela$fold <- sample.int(n_folds, nrow(tanner.chela), replace = TRUE)

# Make mesh
mat.msh <- sdmTMB::make_mesh(mod.dat, c("LONGITUDE","LATITUDE"), n_knots = 200, type = "kmeans")

mod.1 <- sdmTMB(
  MATURE ~ s(SIZE_5MM, k = 10) + YEAR_F, #the 0 is there's a factor predictor for each YEAR, no intercept
  spatial = "on",
  spatiotemporal = "iid",
  mesh = mat.msh,
  family = binomial(),
  time = "YEAR",
  extra_time = xtra.time,
  anisotropy = TRUE,
  data = tanner.chela
)
saveRDS(mod.1, "./Maturity data processing/Doc/Tanner models/sdmTMB_nospVAR_k200.rda")

cv.1 <- sdmTMB_cv(
  MATURE ~ s(SIZE_5MM, k = 10) + YEAR_F, #the 0 is there's a factor predictor for each YEAR, no intercept
  spatial = "on",
  spatiotemporal = "iid",
  mesh = mat.msh,
  family = binomial(),
  time = "YEAR",
  extra_time = xtra.time,
  anisotropy = TRUE,
  data = tanner.chela,
  k_folds = n_folds,
  fold_ids = tanner.chela$fold
)

mod.2 <- sdmTMB(
  MATURE ~ s(SIZE_5MM, k = 10) + YEAR_SCALED, #the 0 is there's a factor predictor for each YEAR, no intercept
  spatial = "on",
  spatiotemporal = "iid",
  mesh = mat.msh,
  family = binomial(),
  time = "YEAR",
  spatial_varying = ~ 0 + SIZE_5MM,
  extra_time = xtra.time,
  anisotropy = TRUE,
  data = tanner.chela
)

saveRDS(mod.2, "./Maturity data processing/Doc/Tanner models/sdmTMB_spVAR_SIZE_k200.rda")

cv.2 <- sdmTMB_cv(
  MATURE ~ s(SIZE_5MM, k = 10) + YEAR_SCALED, #the 0 is there's a factor predictor for each YEAR, no intercept
  spatial = "on",
  spatiotemporal = "iid",
  mesh = mat.msh,
  family = binomial(),
  time = "YEAR",
  extra_time = xtra.time,
  spatial_varying = ~ 0 + SIZE_5MM,
  anisotropy = TRUE,
  data = tanner.chela,
  k_folds = n_folds,
  fold_ids = tanner.chela$fold
)

mat.msh <- sdmTMB::make_mesh(mod.dat, c("LONGITUDE","LATITUDE"), n_knots = 300, type = "kmeans")
mod.3 <- sdmTMB(
  MATURE ~ s(SIZE_5MM, k = 10) + YEAR_F, #the 0 is there's a factor predictor for each YEAR, no intercept
  spatial = "on",
  spatiotemporal = "iid",
  mesh = mat.msh,
  family = binomial(),
  time = "YEAR",
  extra_time = xtra.time,
  anisotropy = TRUE,
  data = tanner.chela
)

saveRDS(mod.3, "./Maturity data processing/Doc/Tanner models/sdmTMB_nospVAR_k300.rda")

cv.3 <- sdmTMB_cv(
  MATURE ~ s(SIZE_5MM, k = 10) + YEAR_F, #the 0 is there's a factor predictor for each YEAR, no intercept
  spatial = "on",
  spatiotemporal = "iid",
  mesh = mat.msh,
  family = binomial(),
  time = "YEAR",
  extra_time = xtra.time,
  anisotropy = TRUE,
  data = tanner.chela,
  k_folds = n_folds,
  fold_ids = tanner.chela$fold
)

mod.4 <- sdmTMB(
  MATURE ~ s(SIZE_5MM, k = 10) + YEAR_SCALED, #the 0 is there's a factor predictor for each YEAR, no intercept
  spatial = "on",
  spatiotemporal = "iid",
  mesh = mat.msh,
  family = binomial(),
  time = "YEAR",
  spatial_varying = ~ 0 + SIZE_5MM,
  extra_time = xtra.time,
  anisotropy = TRUE,
  data = tanner.chela
)

saveRDS(mod.4, "./Maturity data processing/Doc/Tanner models/sdmTMB_spVAR_SIZE_k300.rda")

cv.4 <- sdmTMB_cv(
  MATURE ~ s(SIZE_5MM, k = 10) + YEAR_SCALED, #the 0 is there's a factor predictor for each YEAR, no intercept
  spatial = "on",
  spatiotemporal = "iid",
  mesh = mat.msh,
  family = binomial(),
  time = "YEAR",
  extra_time = xtra.time,
  spatial_varying = ~ 0 + SIZE_5MM,
  anisotropy = TRUE,
  data = tanner.chela,
  k_folds = n_folds,
  fold_ids = tanner.chela$fold
)




par.results <- data.frame(AICc(mod.1, mod.2, mod.3, mod.4), 
                          pass_sanity = c("Y", "Y", "Y", "Y"),
                          logLik = c(sum(cv.1$fold_loglik), sum(cv.2$fold_loglik), sum(cv.3$fold_loglik), sum(cv.4$fold_loglik)),
                          knots = c(200, 200, 300, 300),
                          terms = c("YEAR_F, no sptemp var",
                                    "YEAR_SCALED, spatial var of size",
                                    "YEAR_F, no sptemp var",
                                    "YEAR_SCALED, spatial var of size"))
write.csv(par.results, "./Maturity data processing/Doc/sdmTMB_terms_CV_TANNER.csv")


mod.1 <- readRDS("./Maturity data processing/Doc/Tanner models/sdmTMB_nospVAR_k200.rda")
mod.2 <- readRDS("./Maturity data processing/Doc/Tanner models/sdmTMB_spVAR_SIZE_k200.rda")
mod.3 <- readRDS("./Maturity data processing/Doc/Tanner models/sdmTMB_nospVAR_k300.rda")
mod.4 <- readRDS("./Maturity data processing/Doc/Tanner models/sdmTMB_spVAR_SIZE_k300.rda")



sanity(mod.1)
sanity(mod.2)
sanity(mod.3)
sanity(mod.4)



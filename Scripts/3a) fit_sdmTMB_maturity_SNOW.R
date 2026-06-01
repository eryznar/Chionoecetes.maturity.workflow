# PURPOSE: to fit sdmTMB models to estimate p(mat) for snow crab

# Author: Emily Ryznar

# NOTES:


# LOAD LIBS/PARAMS -----
source("./Scripts/Sourced scripts/load_libs_params.R")

# LOAD CHELA DATA -----
snow.chela <-  read.csv("./Data/chionoecetes_chela_withcutlines.csv") %>% #already filtered appropriately with cutline applied
                  dplyr::select(!X) %>%
                  filter(SPECIES == "SNOW") 

# FIT MATURITY SDMTMB ----
# Make mesh
mat.msh <- sdmTMB::make_mesh(snow.chela, c("LONGITUDE","LATITUDE"), n_knots = 300, type = "kmeans")

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
  data = snow.chela,
  silent = FALSE
)

saveRDS(snow.mod, "./Models/snow_maturity_sdmTMB.rda")


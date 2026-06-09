# PURPOSE: to explore fitting sdmTMB models to estimate p(mat) for Chionoecetes

# Author: Emily Ryznar

# NOTES:


# LOAD LIBS/PARAMS -----
source("./Scripts/Sourced scripts/load_libs_params.R")

# LOAD CHELA DATA -----
tanner.chela <-  read.csv("./Data/chionoecetes_chela_withcutlines.csv") %>% #already filtered appropriately with cutlines applied
                    dplyr::select(!X) %>%
                    filter(SPECIES == "TANNER") 

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

# Make sure it passes sanity check
sanity(tanner.mod)

# Save
saveRDS(tanner.mod, "./Models/tanner_maturity_sdmTMB.rda")

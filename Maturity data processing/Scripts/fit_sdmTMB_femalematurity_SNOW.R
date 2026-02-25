## ------------------------------------------------------------
# PURPOSE: to analyze maturity patterns in relation to biological and environmental drivers

# Author: Emily Ryznar

# NOTES:
# Decision points:


# LOAD LIBS/PARAMS ---------------------------------------------------------------------------------------
source("./Maturity research/Scripts/load_libs_params.R")

# LOAD DATA AND PROCESS ----------------------------------------------------------------------------------
spec.dat <- readRDS("./Maturity research/Data/snow_survey_specimenEBS.rda")


mod.dat <- spec.dat$specimen %>% 
              filter(SEX == 2) %>%
              mutate(MATURE  = case_when(CLUTCH_SIZE >0 ~ 1,
                                         TRUE ~ 0)) 

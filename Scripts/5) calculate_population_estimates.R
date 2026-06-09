# PURPOSE: to generate ogives, SAM, and mature bio/abund for sdmTMB maturity models and compare to legacy workflow

# Author: Emily Ryznar

# NOTES:


# LOAD LIBS/PARAMS -----
source("./Scripts/Sourced scripts/load_libs_params.R")

# Source function to simulate from sdmTMB model, and calculate ogives, SAM, and mature bioabund with uncertainty
source("./Scripts/Sourced scripts/calc_maturepop_estimates_function.R")

# SNOW CRAB ----
snow_mod <- readRDS("./Models/snow_maturity_sdmTMB.rda")
snow_dat <- readRDS("./Data/snow_survey_specimenEBS.rda")

snow_out <- calc_maturepop_estimates(model = snow_mod, 
                                     crab_data = snow_dat, 
                                     years = c(1989:2019, 2021:2025), 
                                     species = "SNOW", 
                                     region = "EBS",
                                     district = NULL, 
                                     size_1mm = NULL, 
                                     size_min = NULL, 
                                     size_max = NULL, 
                                     fill_missing_years = TRUE, 
                                     output = c("ogives", "SAM", "bioabund"),
                                     n_sim = 5) 
    
# Plot ogives
ggplot(snow_out$ogives, aes(SIZE_5MM, PROP_MATURE, color = DISTRICT))+
  geom_line()+
  facet_wrap(~YEAR)+
  geom_ribbon(aes(SIZE_5MM, ymax = PROP_MATURE+PROP_MATURE_SD*1.96, 
                  ymin = PROP_MATURE-PROP_MATURE_SD*1.96, fill = DISTRICT), alpha = 0.25, color = NA)+
  theme_bw()

snow_ogives <- snow_out$ogives

write.csv(snow_ogives, "./Output/snow_ogives_example.csv")

snow_ogives %>%
    filter(SIZE_5MM >= 27.5 & SIZE_5MM <= 132.5) %>%
    dplyr::select(YEAR, SIZE_5MM, PROP_MATURE) %>%
    arrange(SIZE_5MM) %>%
    tidyr::pivot_wider(names_from  = SIZE_5MM,
                       values_from = PROP_MATURE)

# Plot SAM
ggplot(snow_out$SAM, aes(YEAR, SAM))+
  geom_line()+
  facet_wrap(~DISTRICT)+
  geom_errorbar(aes(YEAR, ymax = SAM+SAM_SD*1.96, 
                    ymin = SAM-SAM_SD*1.96), alpha = 0.25)+
  theme_bw()

# Plot bioabund
ggplot(snow_out$mature_bioabund  %>% filter(CATEGORY == "Mature male", is.na(DISTRICT) == FALSE), 
       aes(YEAR, ABUNDANCE/1e6))+
  geom_line()+
  facet_wrap(~DISTRICT, nrow = 2, scales = "free")+
  theme_bw()+
  geom_ribbon(aes(YEAR, ymin = ABUNDANCE_lo/1e6, ymax = ABUNDANCE_hi/1e6), alpha = 0.25)


  
# TANNER CRAB ----
# Run function
tanner_mod <- readRDS("./Models/tanner_maturity_sdmTMB.rda")
tanner_dat <- readRDS("./Data/tanner_survey_specimenEBS.rda")

tanner_out <- calc_maturepop_estimates(model = tanner_mod, 
                                       crab_data = tanner_dat, 
                                       years = c(1990:2019, 2021:2025), 
                                       species = "TANNER", 
                                       region = "EBS",
                                       district = NULL, # using only "ALL" won't work for CPUE or bioabund until crabpack includes EBS-wide option for Tanners
                                       size_1mm = NULL, 
                                       size_min = NULL, 
                                       size_max = NULL, 
                                       fill_missing_years = TRUE, 
                                       output = c("ogives"),
                                       n_sim = 5) 
# Plot ogives
ggplot(tanner_out$ogives, aes(SIZE_5MM, PROP_MATURE, color = DISTRICT))+
  geom_line()+
  facet_wrap(~YEAR)+
  geom_ribbon(aes(SIZE_5MM, ymax = PROP_MATURE+PROP_MATURE_SD*1.96, 
                  ymin = PROP_MATURE-PROP_MATURE_SD*1.96, fill = DISTRICT), alpha = 0.25, color = NA)+
  theme_bw()

# Plot SAM
ggplot(tanner_out$SAM, aes(YEAR, SAM))+
  geom_line()+
  facet_wrap(~DISTRICT, nrow = 3)+
  geom_errorbar(aes(YEAR, ymax = SAM+SAM_SD*1.96, 
                  ymin = SAM-SAM_SD*1.96), alpha = 0.25)+
  theme_bw()

# Plot bioabund
ggplot(tanner_out$mature_bioabund  %>% filter(CATEGORY == "Mature male", is.na(DISTRICT) == FALSE), 
       aes(YEAR, ABUNDANCE/1e6))+
  geom_line()+
  facet_wrap(~DISTRICT, nrow = 2, scales = "free")+
  theme_bw()+
  geom_ribbon(aes(YEAR, ymin = ABUNDANCE_lo/1e6, ymax = ABUNDANCE_hi/1e6), alpha = 0.25)

              
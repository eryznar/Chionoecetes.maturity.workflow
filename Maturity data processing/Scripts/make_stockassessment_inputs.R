# PURPOSE: to generate maturity inputs for stock assessment models

# Author: Emily Ryznar

# NOTES:


# LOAD LIBS/PARAMS ---------------------------------------------------------------------------------------
source("./Maturity data processing/Scripts/load_libs_params.R")

# Source function to simulate from sdmTMB model, and calculate ogives, SAM, and mature bioabund with uncertainty
source("./Maturity data processing/Scripts/calc_maturepop_estimates_function.R")


# FOR BUCK ----
cp_ratio <- crabpack::get_male_maturity(species = "TANNER", region = "EBS")$male_mat_ratio


ggplot(cp_ratio, aes(SIZE_BIN, PROP_MATURE))+
  geom_line()+
  facet_wrap(~YEAR)

# Specify function parameters
tanner_mod <- readRDS("./Maturity data processing/Doc/Tanner models/sdmTMB_spVAR_noBIN_k200.rda")
tanner_dat <- readRDS("./Maturity data processing/Data/tanner_survey_specimenEBS.rda")
tanner_yrs <- c(1990:2019, 2021:2025)
species <- "TANNER"
output <- c("ogives")
region <- "EBS"
district <- "ALL"
size_1mm <- FALSE

# Run function 
calc_maturepop_estimates(tanner_mod, tanner_dat, tanner_yrs, species, region, district, size_1mm, output) -> tanner.out

saveRDS(tanner.out, "./Maturity data processing/Output/tanner.outputforBuck.rda")
t.out <- readRDS("./Maturity data processing/Output/tanner.outputforBuck.rda")

t.out$ogives %>%
  filter(YEAR >=1990) %>%
  dplyr::select(SPECIES, REGION, DISTRICT, YEAR, SIZE_5MM, NUM_IMMATURE, 
                NUM_MATURE, TOTAL_CRAB, PROP_MATURE_mean, PROP_MATURE_sd) %>%
  rename(PROP_MATURE = PROP_MATURE_mean, PROP_MATURE_SD = PROP_MATURE_sd, SIZE_BIN = SIZE_5MM) %>%
  filter(!YEAR %in% c(2008, 2012, 2014, 2016)) -> for.buck
  

write.csv(for.buck, "./Maturity data processing/Output/TANNER_male_mat_ratio.csv")

ggplot(for.buck, aes(SIZE_BIN, PROP_MATURE))+
  geom_line()+
  scale_y_continuous(breaks = seq(0, 1, by = 0.25))+
  facet_wrap(~YEAR)+
  theme_bw()


for.buck %>%
  mutate(pmat = NUM_MATURE/TOTAL_CRAB) -> pp

ggplot()+
  geom_line(pp, mapping = aes(SIZE_BIN, PROP_MATURE), color = "blue")+ # making sure PROP_MATURE is correct
  geom_line(pp, mapping = aes(SIZE_BIN, pmat), color = "green")+
  facet_wrap(~YEAR)
  
sub1 %>%
  group_by(YEAR, SIZE_5MM) %>%
  reframe(TOTAL_CRAB = sum(SAMPLING_FACTOR)) -> ll

ggplot()+
  geom_line(pp, mapping = aes(SIZE_5MM, TOTAL_CRAB), color = "blue")+ # making sure sampling factor is correct
  geom_line(ll, mapping = aes(SIZE_5MM, TOTAL_CRAB), color = "green")+
  facet_wrap(~YEAR)

# FOR CODY ----
# Specify function parameters
snow_mod <- readRDS("./Maturity data processing/Doc/Snow models/sdmTMB_spVAR_noBIN_k300.rda")
snow_dat <- readRDS("./Maturity data processing/Data/snow_survey_specimenEBS.rda")
snow_yrs <- c(1989:2019, 2021:2025)
species <- "SNOW"
output <- "ogives"
region <- "EBS"
district <- "ALL"
size_1mm <- FALSE

# Run function 
calc_maturepop_estimates(snow_mod, snow_dat, snow_yrs, species, region, district, size_1mm, output) -> snow.out

saveRDS(snow.out, "./Maturity data processing/Output/snow.outputforCody.rda")

snow.ogives <- snow.out$ogives %>%
                    filter(!YEAR %in% c(2008, 2012, 2014, 2016, 2020), # omitting no maturity years
                           SIZE_5MM < 172.5) # filtering large Tanner crab classified as snow in 2025

for.cody <- snow.ogives %>%
  filter(SIZE_5MM >= 27.5 & SIZE_5MM <= 132.5) %>%
  dplyr::select(YEAR, SIZE_5MM, PROP_MATURE_mean) %>%
  rename(PROP_MATURE = PROP_MATURE_mean) %>%
  arrange(SIZE_5MM) %>%
  tidyr::pivot_wider(
    names_from  = SIZE_5MM,
    values_from = PROP_MATURE
  )

# convert to a matrix with rownames = YEAR, colnames = sizes
ogive_mat <- as.matrix(for.cody[ , -1])            # drop YEAR column
rownames(ogive_mat) <- for.cody$YEAR              # YEARS as rownames
colnames(ogive_mat) <- as.numeric(colnames(ogive_mat))
  
write.csv(ogive_mat, "./Maturity data processing/Output/SNOW_male_pmolt_array.csv")





ggplot()+
  scale_y_continuous(breaks = seq(0, 1, by = 0.25))+
  geom_ribbon(snow.ogives, mapping = aes(SIZE_5MM, ymin = PROP_MATURE_lo, ymax = PROP_MATURE_hi), alpha = 0.5, color = NA, fill = "cadetblue")+
  geom_line(snow.ogives, mapping = aes(SIZE_5MM, PROP_MATURE_mean), color = "cadetblue", linewidth = 1)+
  geom_hline(yintercept = 0.5, linetype = "dashed")+
  facet_wrap(~YEAR)+
  theme_bw()

ggplot(snow.out$SAM, aes(YEAR, SAM_mean))+
  geom_point()+
  geom_line()+
  geom_smooth(method = "lm")

summary(lme(SAM_mean ~ YEAR, data = na.omit(snow.out$SAM), random = ~ 1 | YEAR, correlation = corAR1()))
         
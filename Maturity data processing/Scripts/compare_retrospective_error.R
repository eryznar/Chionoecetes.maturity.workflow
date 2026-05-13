# PURPOSE: to explore fitting sdmTMB models to estimate p(mat) for Chionoecetes

# Author: Emily Ryznar

# NOTES:


# LOAD LIBS/PARAMS ---------------------------------------------------------------------------------------
source("./Maturity data processing/Scripts/load_libs_params.R")

# Source function to simulate from sdmTMB model, and calculate ogives, SAM, and mature bioabund with uncertainty
source("./Maturity data processing/Scripts/calc_maturepop_estimates_function.R")


# LOAD DATA ----------------------------------------------------------------------------------------------
recent_years <- c(2022:2025)
# SNOW ----
# Load data ----
# Load minima data, calculate cutline params
snow_minima <- read.csv("./Maturity data processing/Output/chela_cutline_minima.csv") %>%
  filter(SPECIES == "SNOW") %>%
  mutate(BETA0 = coef(lm(MINIMUM ~ MIDPOINT))[1],
         BETA1 = coef(lm(MINIMUM ~ MIDPOINT))[2])

BETA0 <- unique(snow_minima$BETA0)
BETA1 <- unique(snow_minima$BETA1)

# Load snow specimen data (subsample 1)
snow.specimen<-  readRDS("./Maturity data processing/Data/snow_survey_specimenEBS.rda")$specimen %>%
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
         YEAR_F = as.factor(YEAR)) %>%
  rename(SIZE_5MM = SIZE_BINNED) %>%
  dplyr::select(YEAR, YEAR_F, STATION_ID, LATITUDE, LONGITUDE, SIZE_5MM, SIZE_CATEGORY, 
                SAMPLING_FACTOR) %>%
  filter(!YEAR %in% c(2012, recent_years))

# Load snow chela data from chela database (subsample 2)
snow.chela <-  read.csv("./Maturity data processing/Data/snow_tanner_cheladatabase.csv") %>% #already filtered appropriately
  dplyr::select(!X) %>%
  filter(SPECIES == "SNOW", SIZE >= 35 & SIZE <= 135) %>%
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
  # MATURE = case_when((SIZE <=35) ~ 0,
  #                    (SIZE >= 135) ~ 1,
  #TRUE ~ MATURE)) %>%
  st_as_sf(., coords = c("LONGITUDE", "LATITUDE"), crs = "+proj=longlat +datum=WGS84") %>%
  st_transform(., crs = "+proj=utm +zone=2") %>%
  cbind(st_coordinates(.)) %>%
  as.data.frame(.) %>%
  #filter(N_chela > 40) %>% #necessary for consecutive size correlations to run
  mutate(LATITUDE = Y/1000, # scale to km so values don't get too large
         LONGITUDE = X/1000,
         SIZE_CATEGORY = as.factor(paste0("SIZE", SIZE_BINNED)),
         YEAR_F = as.factor(YEAR),
         YEAR_SCALED = as.numeric(scale(YEAR)),
         MATURE = case_when((SIZE <=35) ~ 0,
                            (SIZE >= 135) ~ 1,
                            TRUE ~ MATURE)) %>%
  as.data.frame(.) %>%
  rename(SIZE_5MM = SIZE_BINNED) %>%
  dplyr::select(YEAR, YEAR_F, YEAR_SCALED, STATION_ID, LATITUDE, LONGITUDE, SIZE, SIZE_5MM, SIZE_CATEGORY, MATURE) %>%
  filter(!YEAR %in% c(2012, recent_years))

# Fit models ----
mod.dat <- snow.chela

set.seed(1)

# Set params
mat.msh <- sdmTMB::make_mesh(mod.dat, c("LONGITUDE","LATITUDE"), n_knots = 300, type = "kmeans")

xtra.time <- c(2008, 2012, 2014, 2016, 2020) # missing years across all size bins

mod.2 <- sdmTMB(
  MATURE ~ s(SIZE, k = 13) + YEAR_SCALED, 
  spatial = "on",
  spatiotemporal = "iid",
  mesh = mat.msh,
  family = binomial(),
  time = "YEAR",
  spatial_varying = ~ 0 + SIZE,
  extra_time = xtra.time,
  anisotropy = TRUE,
  data = mod.dat
)

saveRDS(mod.2, "./Maturity data processing/Doc/Snow models/sdmTMB_spVAR_noBIN_k300_retrospective.rda")




# Calc pop estimates ----
model <- readRDS("./Maturity data processing/Doc/Snow models/sdmTMB_spVAR_noBIN_k300_retrospective.rda")
snow_dat <- readRDS("./Maturity data processing/Data/snow_survey_specimenEBS.rda")
snow_yrs <- c(1989:2019, 2021) # omitting recent years
crab_data = snow_dat
years = snow_yrs
species = "SNOW"
output = NULL
district = "ALL"
size_1mm = NULL
region = "EBS"


s.out <- calc_maturepop_estimates(model, crab_data, years, species, region, district, size_1mm, output)
saveRDS(s.out, "./Maturity data processing/Doc/snow_matpop_retrospective.rda")


# TANNER ----
# Load data ----
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
  dplyr::select(YEAR, YEAR_F, YEAR_SCALED, DISTRICT, STATION_ID, LATITUDE, LONGITUDE, SIZE, SIZE_5MM, SIZE_CATEGORY, 
                SAMPLING_FACTOR) %>%
  filter(!YEAR %in% recent_years)

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
  dplyr::select(YEAR, DISTRICT, YEAR_F, YEAR_SCALED, STATION_ID, LATITUDE, LONGITUDE, SIZE, SIZE_5MM, SIZE_CATEGORY, MATURE) %>%
  filter(!YEAR %in% recent_years)

# Fit models ----
set.seed(1)
mod.dat <- tanner.chela
# Set params
mat.msh <- sdmTMB::make_mesh(mod.dat, c("LONGITUDE","LATITUDE"), n_knots = 200, type = "kmeans")
# Extra time
xtra.time <- c(2013, 2015, 2020) # missing years across all size bins

mod.2 <- sdmTMB(
  MATURE ~ s(SIZE, k = 10) + YEAR_SCALED, 
  spatial = "on",
  spatiotemporal = "iid",
  mesh = mat.msh,
  family = binomial(),
  time = "YEAR",
  spatial_varying = ~ 0 + SIZE,
  extra_time = xtra.time,
  anisotropy = TRUE,
  data = mod.dat
)

saveRDS(mod.2, "./Maturity data processing/Doc/Tanner models/sdmTMB_spVAR_noBIN_k200_retrospective.rda")

# Calc pop estimates ----
model <- readRDS("./Maturity data processing/Doc/Tanner models/sdmTMB_spVAR_noBIN_k200_retrospective.rda")
tanner_dat <- readRDS("./Maturity data processing/Data/tanner_survey_specimenEBS.rda")
tanner_yrs <- c(1990:2019, 2021) # omitting recent years
crab_data = tanner_dat
years = tanner_yrs
species = "TANNER"
output = c("ogives", "SAM", "bioabund")
district = "E166"
size_1mm = NULL
region = "EBS"


t.out.east <- calc_maturepop_estimates(model, crab_data, years, species, region, district, size_1mm, output)

district <- "W166"
t.out.west <- calc_maturepop_estimates(model, crab_data, years, species, region, district, size_1mm, output)

saveRDS(t.out.east, "./Maturity data processing/Doc/tannerE_matpop_retrospective.rda")
saveRDS(t.out.west, "./Maturity data processing/Doc/tannerW_matpop_retrospective.rda")

# Compare estimates----
# Snow crab
sr <- readRDS("./Maturity data processing/Doc/snow_matpop_retrospective.rda")
s <- readRDS("./Maturity data processing/Doc/snow_matpop.rda")

# Snow ogives
sr$ogives %>%
  filter(SIZE_5MM <=150) -> sor
s$ogives %>%
  filter(SIZE_5MM <=150) -> so

ggplot()+
  geom_line(sor, mapping = aes(SIZE_5MM, PROP_MATURE_mean, color = as.factor(1)), linewidth = 0.75)+
  geom_ribbon(sor, mapping = aes(SIZE_5MM, ymin = PROP_MATURE_lo, ymax = PROP_MATURE_hi, fill = as.factor(1)), color = NA, alpha = 0.35)+
  geom_line(so, mapping = aes(SIZE_5MM, PROP_MATURE_mean, color = as.factor(2)), linewidth = 0.75)+
  geom_ribbon(so , mapping = aes(SIZE_5MM, ymin = PROP_MATURE_lo, ymax = PROP_MATURE_hi, fill = as.factor(2)), color = NA, alpha = 0.35)+
  facet_wrap(~YEAR)+
  scale_fill_manual(values = c("blue", "salmon"), labels = c("≤2021", "≤2025"), name = "")+
  scale_color_manual(values = c( "blue", "salmon"), labels = c("≤2021", "≤2025"), name = "")+
  theme_bw()+
  ggtitle("Snow")+
  ylab("Proportion mature")+
  xlab("Carapace width (mm)")+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 10), 
        legend.text  = element_text(size = 12))

ggsave("./Maturity data processing/Doc/snow_retrospective_ogives.png", width = 8, height = 7)

# Snow bioabund
sr$mature_bioabund %>% filter(CATEGORY == "Mature male", YEAR >1988 & YEAR <2022) %>%
  full_join(data.frame(YEAR = 2020)) -> pp
s$mature_bioabund %>% filter(CATEGORY == "Mature male", YEAR >1988) %>%
  full_join(data.frame(YEAR = 2020)) -> tt


ggplot()+
  geom_line(pp, mapping = aes(YEAR, ABUNDANCE, color = as.factor(1)), linewidth = 0.75)+
  geom_ribbon(pp, mapping = aes(YEAR, ymin = ABUNDANCE - ABUNDANCE_CI,
                                ymax = ABUNDANCE + ABUNDANCE_CI,  fill = as.factor(1)), alpha = 0.25)+
  geom_line(tt, mapping = aes(YEAR, ABUNDANCE, color = as.factor(2)), linewidth = 0.75)+
  geom_ribbon(tt, mapping = aes(YEAR, ymin = ABUNDANCE - ABUNDANCE_CI,
                                ymax = ABUNDANCE + ABUNDANCE_CI, fill = as.factor(2)), alpha = 0.25)+
  geom_point(pp %>% filter(YEAR %in% c(2013, 2015)),
             mapping = aes(YEAR, ABUNDANCE, color = as.factor(1)))+
  geom_errorbar(pp %>% filter(YEAR %in% c(2013, 2015)),
                mapping = aes(x = YEAR, ymin = ABUNDANCE - ABUNDANCE_CI, ymax = ABUNDANCE + ABUNDANCE_CI, 
                              color= as.factor(1)), linewidth = 1)+
  geom_point(tt %>% filter(YEAR %in% c(2013, 2015)),
             mapping = aes(YEAR, ABUNDANCE, color = as.factor(2)))+
  geom_errorbar(tt %>% filter(YEAR %in% c(2013, 2015)),
                mapping = aes(x = YEAR, ymin = ABUNDANCE - ABUNDANCE_CI, ymax = ABUNDANCE + ABUNDANCE_CI, 
                              color= as.factor(2)), linewidth = 1)+
  scale_fill_manual(values = c("blue", "salmon"), labels = c("≤2021", "≤2025"), name = "")+
  scale_color_manual(values = c( "blue", "salmon"), labels = c("≤2021", "≤2025"), name = "")+
  theme_bw()+
  ggtitle("Snow")+
  ylab("Abundance (millions)")+
  xlab("Year")+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 10))

ggsave("./Maturity data processing/Doc/snow_retrospective_bioabund.png", width = 7, height = 5)


# Tanner 
twr <- readRDS("./Maturity data processing/Doc/tannerW_matpop_retrospective.rda")
tw <- readRDS("./Maturity data processing/Doc/tannerW_matpop.rda")

# Tanner ogives
twr$ogives -> or
tw$ogives -> o

ggplot()+
  geom_line(or %>% filter(DISTRICT == "W166"), mapping = aes(SIZE_5MM, PROP_MATURE_mean, color = as.factor(1)), linewidth = 0.75)+
  geom_ribbon(or %>% filter(DISTRICT == "W166"), mapping = aes(SIZE_5MM, ymin = PROP_MATURE_lo, ymax = PROP_MATURE_hi, fill = as.factor(1)), color = NA, alpha = 0.35)+
  geom_line(o %>% filter(DISTRICT == "W166"), mapping = aes(SIZE_5MM, PROP_MATURE_mean, color = as.factor(2)), linewidth = 0.75)+
  geom_ribbon(o %>% filter(DISTRICT == "W166"), mapping = aes(SIZE_5MM, ymin = PROP_MATURE_lo, ymax = PROP_MATURE_hi, fill = as.factor(2)), color = NA, alpha = 0.35)+
  facet_wrap(~YEAR)+
  scale_fill_manual(values = c("blue", "salmon"), labels = c("≤2021", "≤2025"), name = "")+
  scale_color_manual(values = c( "blue", "salmon"), labels = c("≤2021", "≤2025"), name = "")+
  theme_bw()+
  ggtitle("Tanner West")+
  ylab("Proportion mature")+
  xlab("Carapace width (mm)")+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 10), 
        legend.text  = element_text(size = 12))

ggsave("./Maturity data processing/Doc/TannerW_retrospective_ogives.png", width = 8, height = 7)

ggplot()+
  geom_line(or %>% filter(DISTRICT == "E166"), mapping = aes(SIZE_5MM, PROP_MATURE_mean, color = as.factor(1)), linewidth = 0.75)+
  geom_ribbon(or %>% filter(DISTRICT == "E166"), mapping = aes(SIZE_5MM, ymin = PROP_MATURE_lo, ymax = PROP_MATURE_hi, fill = as.factor(1)), color = NA, alpha = 0.35)+
  geom_line(o %>% filter(DISTRICT == "E166"), mapping = aes(SIZE_5MM, PROP_MATURE_mean, color = as.factor(2)), linewidth = 0.75)+
  geom_ribbon(o %>% filter(DISTRICT == "E166"), mapping = aes(SIZE_5MM, ymin = PROP_MATURE_lo, ymax = PROP_MATURE_hi, fill = as.factor(2)), color = NA, alpha = 0.35)+
  facet_wrap(~YEAR)+
  scale_fill_manual(values = c("blue", "salmon"), labels = c("≤2021", "≤2025"), name = "")+
  scale_color_manual(values = c( "blue", "salmon"), labels = c("≤2021", "≤2025"), name = "")+
  theme_bw()+
  ggtitle("Tanner East")+
  ylab("Proportion mature")+
  xlab("Carapace width (mm)")+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 10), 
        legend.text  = element_text(size = 12))

ggsave("./Maturity data processing/Doc/TannerE_retrospective_ogives.png", width = 8, height = 7)


# Tanner bioabund
twr$mature_bioabund %>% filter(CATEGORY == "Mature male", 
                               YEAR >1989 & YEAR <2022, is.na(DISTRICT)==FALSE) %>%
  full_join(data.frame(YEAR = 2020, DISTRICT = c(unique(.$DISTRICT)))) -> pp
tw$mature_bioabund %>% filter(CATEGORY == "Mature male", YEAR >1989, is.na(DISTRICT)==FALSE) %>%
  full_join(data.frame(YEAR = 2020, DISTRICT = c(unique(.$DISTRICT)))) -> tt


ggplot()+
  geom_line(pp, mapping = aes(YEAR, ABUNDANCE, color = as.factor(1)), linewidth = 0.75)+
  geom_ribbon(pp, mapping = aes(YEAR, ymin = ABUNDANCE - ABUNDANCE_CI,
                                ymax = ABUNDANCE + ABUNDANCE_CI,  fill = as.factor(1)), alpha = 0.25)+
  geom_line(tt, mapping = aes(YEAR, ABUNDANCE, color = as.factor(2)), linewidth = 0.75)+
  geom_ribbon(tt, mapping = aes(YEAR, ymin = ABUNDANCE - ABUNDANCE_CI,
                                ymax = ABUNDANCE + ABUNDANCE_CI, fill = as.factor(2)), alpha = 0.25)+
  geom_point(pp %>% filter(DISTRICT == "E166", YEAR %in% c(2012, 2014)),
             mapping = aes(YEAR, ABUNDANCE, color = as.factor(1)))+
  geom_errorbar(pp %>% filter(DISTRICT == "E166", YEAR %in% c(2012, 2014)),
                mapping = aes(x = YEAR, ymin = ABUNDANCE - ABUNDANCE_CI, ymax = ABUNDANCE + ABUNDANCE_CI, 
                              color= as.factor(1)), linewidth = 1)+
  geom_point(tt %>% filter(DISTRICT == "E166", YEAR %in% c(2012, 2014)),
             mapping = aes(YEAR, ABUNDANCE, color = as.factor(2)))+
  geom_errorbar(tt %>% filter(DISTRICT == "E166", YEAR %in% c(2012, 2014)),
                mapping = aes(x = YEAR, ymin = ABUNDANCE - ABUNDANCE_CI, ymax = ABUNDANCE + ABUNDANCE_CI, 
                              color= as.factor(2)), linewidth = 1)+
  geom_point(pp %>% filter(DISTRICT == "W166", YEAR %in% c(2014)),
             mapping = aes(YEAR, ABUNDANCE, color = as.factor(1)))+
  geom_errorbar(pp %>% filter(DISTRICT == "W166", YEAR %in% c(2014)),
                mapping = aes(x = YEAR, ymin = ABUNDANCE - ABUNDANCE_CI, ymax = ABUNDANCE + ABUNDANCE_CI, 
                              color= as.factor(1)), linewidth = 1)+
  geom_point(tt %>% filter(DISTRICT == "W166", YEAR %in% c(2014)),
             mapping = aes(YEAR, ABUNDANCE, color = as.factor(2)))+
  geom_errorbar(tt %>% filter(DISTRICT == "W166", YEAR %in% c(2014)),
                mapping = aes(x = YEAR, ymin = ABUNDANCE - ABUNDANCE_CI, ymax = ABUNDANCE + ABUNDANCE_CI, 
                              color= as.factor(2)), linewidth = 1)+
  facet_wrap(~DISTRICT, nrow = 2)+
  scale_fill_manual(values = c("blue", "salmon"), labels = c("≤2021", "≤2025"), name = "")+
  scale_color_manual(values = c( "blue", "salmon"), labels = c("≤2021", "≤2025"), name = "")+
  theme_bw()+
  ggtitle("Tanner")+
  ylab("Abundance (millions)")+
  xlab("Year")+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 10))

ggsave("./Maturity data processing/Doc/Tanner_retrospective_bioabund.png", width = 7, height = 5)

# LOAD LIBS/PARAMS -----
source("./Scripts/load_libs_params.R")

snow.spec <- readRDS("./Data/snow_survey_specimenEBS.rda")$specimen %>%
  filter(HAUL_TYPE !=17, SEX == 1, SHELL_CONDITION == 2,
         YEAR ==2019, is.na(CHELA_HEIGHT) == FALSE) %>%
  group_by(STATION_ID, SIZE_1MM) %>%
  reframe(N = n()) # make sure filter for males, sh2, not HT17


ggplot(snow.spec, aes(SIZE_1MM, N))+
  geom_bar(stat = "identity")


ggplot(snow.spec, aes(STATION_ID, N))+
  geom_bar(stat = "identity")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))

# SNOW CRAB Proportion sh 2 industry preferred abundance/shell 2 mature male abundance ----
# Calculate biomass/abundance 
ind.pref <- calc_bioabund(crab_data = readRDS("./Data/snow_survey_specimenEBS.rda"),
                          species = "SNOW",
                          years = years, # set in libs/params script
                          sex = "male",
                          shell_condition = "new_hardshell",
                          size_min = 102) %>%
  rename(indpref_abund = ABUNDANCE, indpref_bio = BIOMASS_MT)

sh2male <- calc_bioabund(crab_data = readRDS("./Data/snow_survey_specimenEBS.rda"),
                         species = "SNOW",
                         years = years, # set in libs/params script
                         sex = "male",
                         shell_condition = "new_hardshell")

# matmale_prop <- get_male_maturity(species = "SNOW",
#                                   region = "EBS",
#                                   channel = "API")$male_mat_ratio %>%
#   group_by(YEAR) %>%
#   reframe(tot_imm = sum(NUM_IMMATURE),
#           tot_mat = sum(NUM_MATURE),
#           tot_crab = sum(TOTAL_CRAB),
#           prop_mat = tot_mat/tot_crab)

tm_matabund <- data.frame(tm_matmaleabund = c(485.5, 712.3, 412.2, 366.2, 544.3, 436.2, 599.4, 627.0, 301.1, 231.0, 80.2, 92.7, 226.1, 121.9, 113.5,
                 175.3, 360.9, 202.8, 326.1, NA, 154.7, 318.6, 304.3, NA, 242.5, NA, 157.2, NA, 326.8, 1135.2, 577.1, 117.9,
                 61.6, 57.2, 243.2), YEAR = c(1989:2019, 2021:2024)) %>%
          mutate(tm_matmaleabund = tm_matmaleabund*1e6)

indpref.matmaleprop <- #sh2male %>%
  #dplyr::select(YEAR, ABUNDANCE, BIOMASS_MT) %>%
  # mutate(sh2matmale_abund= prop_mat*ABUNDANCE,
  #        sh2matmale_bio = prop_mat*BIOMASS_MT) %>%
  # dplyr::select(!c(ABUNDANCE, BIOMASS_MT)) %>%
  #right_join(., ind.pref) %>%
  ind.pref %>%
  right_join(., tm_matabund, by ="YEAR") %>%
  as.data.frame() %>%
  na.omit() %>%
  mutate(indpref_matmaleprop = indpref_abund/tm_matmaleabund)

ggplot(indpref.matmaleprop, aes(YEAR, indpref_matmaleprop))+
  geom_point()+
  geom_line()+
  theme_bw()+
  geom_smooth(method = "lm", se = FALSE)+
  ylab("prop_ind_pref")+
  ggtitle("Snow shell 2 male ≥ 102mm abundance/shell 2 mature male abundance") -> p.1


# TANNER CRAB WEST Proportion sh 2 industry preferred abundance/shell 2 mature male abundance ----
# Calculate biomass/abundance 
ind.pref <- calc_bioabund(crab_data = readRDS("./Data/tanner_survey_specimenEBS.rda"),
                          species = "TANNER",
                          years = years, # set in libs/params script
                          sex = "male",
                          district = "W166",
                          shell_condition = "new_hardshell",
                          size_min = 125) %>%
  rename(indpref_abund = ABUNDANCE, indpref_bio = BIOMASS_MT)

sh2male <- calc_bioabund(crab_data = readRDS("./Data/tanner_survey_specimenEBS.rda"),
                         species = "TANNER",
                         years = years, # set in libs/params script
                         sex = "male",
                         district = "W166",
                         shell_condition = "new_hardshell")

# matmale_prop <- get_male_maturity(species = "TANNER",
#                                   region = "EBS",
#                                   district = "W166",
#                                   channel = "API")$male_mat_ratio %>%
#   group_by(YEAR) %>%
#   reframe(tot_imm = sum(NUM_IMMATURE),
#           tot_mat = sum(NUM_MATURE),
#           tot_crab = sum(TOTAL_CRAB),
#           prop_mat = tot_mat/tot_crab)

tm_matabund <- data.frame(tm_matmaleabund = c(65.5, 37.9, 16.4, 19.8, 7.7, 4, 3.1, 1.8, 3.8, 4.6, 3.6, 6.7, 7.2, 12.5,
                                              15.5, 38.2, 35.1, 36.4, 47.3, NA, 29.9, 21.7, NA, 12.5, NA, 50.7, 25.4, 6.5,
                                              12.0, 8, 12.8, 10, 17.5, 33.0), YEAR = c(1990:2019, 2021:2024)) %>%
              mutate(tm_matmaleabund = tm_matmaleabund*1e6)

indpref.matmaleprop <- #right_join(sh2male, matmale_prop) %>%
  # dplyr::select(YEAR, ABUNDANCE, BIOMASS_MT, tot_mat, tot_crab, prop_mat) %>%
  # mutate(sh2matmale_abund= prop_mat*ABUNDANCE,
  #        sh2matmale_bio = prop_mat*BIOMASS_MT) %>%
  # dplyr::select(!c(ABUNDANCE, BIOMASS_MT)) %>%
  # right_join(., ind.pref) %>%
  ind.pref %>%
  right_join(., tm_matabund, by ="YEAR") %>%
  as.data.frame() %>%
  na.omit() %>%
  mutate(indpref_matmaleprop = indpref_abund/tm_matmaleabund)

ggplot(indpref.matmaleprop, aes(YEAR, indpref_matmaleprop))+
  geom_point()+
  geom_line()+
  theme_bw()+
  geom_smooth(method = "lm", se = FALSE)+
  ylab("prop_ind_pref")+
  ggtitle("Tanner West shell 2 male ≥ 125mm abundance/shell 2 mature male abundance") -> p.2

# TANNER CRAB EAST Proportion sh 2 industry preferred abundance/shell 2 mature male abundance ----
# Calculate biomass/abundance 
ind.pref <- calc_bioabund(crab_data = readRDS("./Data/tanner_survey_specimenEBS.rda"),
                          species = "TANNER",
                          years = years, # set in libs/params script
                          sex = "male",
                          district = "E166",
                          shell_condition = "new_hardshell",
                          size_min = 125) %>%
  rename(indpref_abund = ABUNDANCE, indpref_bio = BIOMASS_MT)

sh2male <- calc_bioabund(crab_data = readRDS("./Data/tanner_survey_specimenEBS.rda"),
                         species = "TANNER",
                         years = years, # set in libs/params script
                         sex = "male",
                         district = "E166",
                         shell_condition = "new_hardshell")

# matmale_prop <- get_male_maturity(species = "TANNER",
#                                   region = "EBS",
#                                   district = "E166",
#                                   channel = "API")$male_mat_ratio %>%
#   group_by(YEAR) %>%
#   reframe(tot_imm = sum(NUM_IMMATURE),
#           tot_mat = sum(NUM_MATURE),
#           tot_crab = sum(TOTAL_CRAB),
#           prop_mat = tot_mat/tot_crab)

tm_matabund <- data.frame(tm_matmaleabund = c(46.6, 38.3, 52.7, 24.4, 12.8, 1.1, 1.2, 1.5, 4.7, 5.9, 8, 5.9, 1.9, 4.8,
                                              6.3, 10.2, 15.4, 15.8, 31.7, NA, 5.9, NA, 21.9, NA, 39.7, NA, 
                                              8.1, 2.9, 0.7, 1.6, 6.9, 13.4, 8.0, 12.5), YEAR = c(1990:2019, 2021:2024)) %>%
  mutate(tm_matmaleabund = tm_matmaleabund*1e6)

indpref.matmaleprop <- #right_join(sh2male, matmale_prop) %>%
  # dplyr::select(YEAR, ABUNDANCE, BIOMASS_MT, tot_mat, tot_crab, prop_mat) %>%
  # mutate(sh2matmale_abund= prop_mat*ABUNDANCE,
  #        sh2matmale_bio = prop_mat*BIOMASS_MT) %>%
  # dplyr::select(!c(ABUNDANCE, BIOMASS_MT)) %>%
  # right_join(., ind.pref) %>%
  ind.pref %>%
  right_join(., tm_matabund, by ="YEAR") %>%
  as.data.frame() %>%
  na.omit() %>%
  mutate(indpref_matmaleprop = indpref_abund/tm_matmaleabund)

ggplot(indpref.matmaleprop, aes(YEAR, indpref_matmaleprop))+
  geom_point()+
  geom_line()+
  theme_bw()+
  geom_smooth(method = "lm", se = FALSE)+
  ylab("prop_ind_pref")+
  ggtitle("Tanner East shell 2 male ≥ 125mm abundance/shell 2 mature male abundance") -> p.3

p.1/p.2/p.3

ggsave("./Figures/chionoecetes_propindpref_matmale.png", width = 8, height = 11)

# SNOW ind preferred biomass vs. df catch biomass ----
# Calculate biomass/abundance 
ind.pref <- calc_bioabund(crab_data = readRDS("./Data/snow_survey_specimenEBS.rda"),
                          species = "SNOW",
                          years = years, # set in libs/params script
                          sex = "male",
                          #shell_condition = "new_hardshell",
                          size_min = 102) %>%
  rename(indpref_abund = ABUNDANCE, indpref_bio = BIOMASS_MT)

df.catch <- data.frame(YEAR = 1982:2023, catch_kt = c(11.85, 12.16, 29.94, 44.45, 46.22, 61.4, 67.79, 73.4, 149.1, 143, 104.7, 67.94,
                                                      34.13, 39.81, 54.22, 114.4, 88.09, 15.1, 11.46, 14.8, 12.84, 10.86, 11.29, 16.77, 
                                                      16.49, 28.59, 26.56, 21.78, 24.61, 40.29, 30.05, 24.49, 30.82, 18.42, 9.67, 8.6,
                                                      12.51, 15.43, 20.41, 2.48, NA, NA))

mod.dat <- right_join(ind.pref, df.catch)

mod <- lm(log(mod.dat$indpref_bio + 10) ~ log(mod.dat$catch_kt + 10), data = mod.dat)

ggplot(mod.dat, aes(log(catch_kt + 10), log(indpref_bio + 10)))+
  geom_point()+
  geom_smooth(method = "lm", se = FALSE) + 
  theme_bw()+
  ylab("log(survey biomass + 10)")+
  xlab("log(commercial catch biomass + 10)")+
  ggtitle("Snow ≥ 102mm survey biomass vs. commercial catch biomass") -> p.4

# Tanner ind preferred biomass vs. df catch biomass ----
# Calculate biomass/abundance 
ind.pref <- calc_bioabund(crab_data = readRDS("./Data/tanner_survey_specimenEBS.rda"),
                          species = "TANNER",
                          years = years, # set in libs/params script
                          sex = "male",
                          #shell_condition = "new_hardshell",
                          size_min = 125) %>%
  rename(indpref_abund = ABUNDANCE, indpref_bio = BIOMASS_MT)

df.catch <- data.frame(YEAR = c(1980:1996, 2005:2023), catch_kt = c(13426, 4990, 2390, 549, 1429, NA, NA, 998, 3180, 11113, 18189, 14424, 15921,
                                                      7666, 3538, 1919, 821, 244.5, 786.7, 861.1, 854.1,
                                                      592.4, NA, NA, NA, 1247.9, 6198, 8878, NA,
                                                      1117.6, 1103.9, NA, 655.2, 493.5, 913.3, 944.5))

mod.dat <- right_join(ind.pref, df.catch)

mod <- lm(log(mod.dat$indpref_bio + 10) ~ log(mod.dat$catch_kt + 10), data = mod.dat)

ggplot(mod.dat, aes(log(catch_kt + 10), log(indpref_bio + 10)))+
  geom_point()+
  geom_smooth(method = "lm", se = FALSE) + 
  theme_bw()+
  ylab("log(survey biomass + 10)")+
  xlab("log(commercial catch biomass + 10)")+
  ggtitle("Tanner ≥ 125mm survey biomass vs. commercial catch biomass") -> p.5

# BBRKC ind preferred biomass vs. df catch biomass ----
# Get specimen data
specimen_data <- crabpack::get_specimen_data(species = "RKC",
                                             region = "EBS",
                                             district = "BB",
                                             years = 1979:2024, # set in libs/params script
                                             channel = "API")
# Calculate biomass/abundance 
ind.pref <- calc_bioabund(crab_data = specimen_data,
                          species = "RKC",
                          years = 1979:2024, # set in libs/params script
                          sex = "male",
                          district = "BB",
                          #shell_condition = "new_hardshell",
                          size_min = 135) %>%
  rename(indpref_abund = ABUNDANCE, indpref_bio = BIOMASS_MT)

df.catch <- data.frame(YEAR = c(1979:2022), catch_kt = c(107.83, 129.95, 33.37, 2.99, NA, 4.08, 4.09, 11.31, 12.29, 7.36, 10.16, 20.44,
                                                         17.18, 8.07, 14.59, NA, NA, 8.5, 8.91, 15, 11.84, 8.24, 8.52, 9.67, 15.73, 15.45,
                                                         18.31, 15.62, 20.37, 20.33, 15.93, 14.83, 7.83, 7.85, 8.6, 9.99, 9.97, 8.47, 6.6,
                                                         4.31, 3.79, 2.65, NA, NA))

mod.dat <- right_join(ind.pref, df.catch)

mod <- lm(log(mod.dat$indpref_bio + 10) ~ log(mod.dat$catch_kt + 10), data = mod.dat)

ggplot(mod.dat, aes(log(catch_kt + 10), log(indpref_bio + 10)))+
  geom_point()+
  geom_smooth(method = "lm", se = FALSE) + 
  theme_bw()+
  ylab("log(survey biomass + 10)")+
  xlab("log(commercial catch biomass + 10)")+
  ggtitle("BBRKC ≥ 135mm survey biomass vs. commercial catch biomass") -> p.6


p.4/p.5

ggsave("./Figures/chionoecetes_indpref_exploitation.png", width = 8, height = 11)



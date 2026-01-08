# PURPOSE: to evaluate female maturity

# Author: Emily Ryznar

# NOTES:


# LOAD LIBS/PARAMS ---------------------------------------------------------------------------------------
source("./Maturity data processing/Scripts/1) load_libs_params.R")

# LOAD SPECIMEN DATA -------------------------------------------------------------------------------------
crab_dat <- readRDS("./Maturity data processing/Data/snow_survey_specimenEBS.rda")

# Add in 5mm bins
fem_spec <-  crab_dat$specimen %>%
  filter(SHELL_CONDITION == 2, SEX == 2) %>%
  mutate(BIN_5MM = cut_width(SIZE_1MM, width = 5, center = 2.5, closed = "left", dig.lab = 4),
         BIN2 = BIN_5MM) %>%
  separate(BIN2, sep = ",", into = c("LOWER", "UPPER")) %>%
  mutate(LOWER = as.numeric(sub('.', '', LOWER)),
         UPPER = as.numeric(gsub('.$', '', UPPER)),
         SIZE_5MM = (UPPER + LOWER)/2,
         MATURE = case_when((CLUTCH_SIZE >0)~ 1,
                            TRUE ~ 0)) %>%
  dplyr::select(!c(BIN_5MM, LOWER, UPPER)) %>%
  mutate(BIN_10MM = cut_width(SIZE_1MM, width = 10, center = 5, closed = "left", dig.lab = 4),
         BIN2 = BIN_10MM) %>%
  separate(BIN2, sep = ",", into = c("LOWER", "UPPER")) %>%
  mutate(LOWER = as.numeric(sub('.', '', LOWER)),
         UPPER = as.numeric(gsub('.$', '', UPPER)),
         SIZE_10MM = (UPPER + LOWER)/2) %>%
  dplyr::select(!c(BIN_10MM, LOWER, UPPER))

# Calculate mean and median size for SAM
fem_spec %>%
  filter(MATURE == 1)%>%
  group_by(YEAR) %>%
  reframe(mean_SAM = mean(SIZE),
          median_SAM = median(SIZE)) %>%
  pivot_longer(., !YEAR, names_to = "SAM_type", values_to = "SAM") %>%
  rbind(., data.frame(YEAR = c(2020, 2020), SAM_type = c("mean_SAM", "median_SAM"), SAM = c(NA, NA))) -> fem_SAM

m1 <- lme(SAM ~ YEAR, fem_SAM %>% na.omit() %>% filter(SAM_type == "mean_SAM"), random = ~ 1 | YEAR, correlation = corAR1())
m2 <- lme(SAM ~ YEAR, fem_SAM %>% na.omit() %>% filter(SAM_type == "median_SAM"), random = ~ 1 | YEAR, correlation = corAR1())

# Plot
ggplot(fem_SAM, aes(YEAR, SAM))+
  geom_line()+
  geom_point()+
  facet_wrap(~SAM_type, ncol = 1)+
  theme_bw()+
  geom_smooth(method = "lm")

ggsave("./Maturity research/Figures/female_SAM_comparison.png", width = 5, height = 6, units = "in")

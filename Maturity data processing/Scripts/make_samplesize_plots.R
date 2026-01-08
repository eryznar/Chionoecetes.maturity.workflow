# PURPOSE: to make sample size plots for Chionoecetes chela measurements

# Author: Emily Ryznar

# NOTES:

# LOAD LIBS/PARAMS ---------------------------------------------------------------------------------------
source("./Maturity data processing/Scripts/1) load_libs_params.R")

# SNOW ----
# Load minima data, calculate cutline params
snow_minima <- read.csv("./Maturity data processing/Output/chela_cutline_minima.csv") %>%
  filter(SPECIES == "SNOW") %>%
  mutate(BETA0 = coef(lm(MINIMUM ~ MIDPOINT))[1],
         BETA1 = coef(lm(MINIMUM ~ MIDPOINT))[2])

BETA0 <- unique(snow_minima$BETA0)
BETA1 <- unique(snow_minima$BETA1)


# Load snow chela data from chela database (subsample 2)
snow_chela <-  read.csv("./Maturity data processing/Data/snow_tanner_cheladatabase.csv") %>% #already filtered appropriately
  dplyr::select(!X) %>%
  filter(SPECIES == "SNOW") %>%
  mutate(CUTOFF = BETA0 + BETA1*LN_CW, # apply cutline model
         MATURE = case_when((LN_CH > CUTOFF) ~ 1,
                            TRUE ~ 0),
         SIZE_1MM = floor(SIZE),
         BIN_5MM = cut_width(SIZE_1MM, width = 5, center = 2.5, closed = "left", dig.lab = 4),
         BIN2 = BIN_5MM) %>%
  mutate(MATURE = case_when((SIZE <=35) ~ 0,
                            (SIZE >= 135) ~ 1,
                            TRUE ~ MATURE)) %>% #applying standard cutoffs
  separate(BIN2, sep = ",", into = c("LOWER", "UPPER")) %>%
  mutate(LOWER = as.numeric(sub('.', '', LOWER)),
         UPPER = as.numeric(gsub('.$', '', UPPER)),
         SIZE_5MM = (UPPER + LOWER)/2) %>%
  dplyr::select(!c(UPPER, LOWER)) %>%
  mutate(BIN_10MM = cut_width(SIZE_1MM, width = 10, center = 5, closed = "left", dig.lab = 4),
         BIN2 = BIN_10MM) %>%
  separate(BIN2, sep = ",", into = c("LOWER", "UPPER")) %>%
  mutate(LOWER = as.numeric(sub('.', '', LOWER)),
         UPPER = as.numeric(gsub('.$', '', UPPER)),
         SIZE_10MM = (UPPER + LOWER)/2) %>%
  dplyr::select(!c(UPPER, LOWER)) %>%
  mutate(BIN_15MM = cut_width(SIZE_1MM, width = 15, center = 7.5, closed = "left", dig.lab = 4),
         BIN2 = BIN_15MM) %>%
  separate(BIN2, sep = ",", into = c("LOWER", "UPPER")) %>%
  mutate(LOWER = as.numeric(sub('.', '', LOWER)),
         UPPER = as.numeric(gsub('.$', '', UPPER)),
         SIZE_15MM = (UPPER + LOWER)/2) %>%
  mutate(BIN_20MM = cut_width(SIZE_1MM, width = 20, center = 10, closed = "left", dig.lab = 4),
         BIN2 = BIN_20MM) %>%
  separate(BIN2, sep = ",", into = c("LOWER", "UPPER")) %>%
  mutate(LOWER = as.numeric(sub('.', '', LOWER)),
         UPPER = as.numeric(gsub('.$', '', UPPER)),
         SIZE_20MM = (UPPER + LOWER)/2) %>%
  dplyr::select(!c(UPPER, LOWER)) %>%
  group_by(YEAR, SIZE_5MM) %>%
  mutate(N_chela_5mm = n()) %>%
  ungroup() %>%
  group_by(YEAR, SIZE_10MM) %>%
  mutate(N_chela_10mm = n()) %>%
  ungroup() %>%
  group_by(YEAR, SIZE_15MM) %>%
  mutate(N_chela_15mm = n()) %>%
  ungroup() %>%
  group_by(YEAR, SIZE_20MM) %>%
  mutate(N_chela_20mm = n()) %>%
  ungroup()


length(unique(snow_chela$SIZE_CATEGORY))

ggplot(snow_chela %>% dplyr::filter(YEAR != 2012), aes(YEAR, SIZE_5MM, fill = N_chela_5mm))+
  geom_tile()+
  theme_bw()+
  ggtitle("Snow")+
  scale_fill_viridis_c(option = "inferno") -> s.heatmap


ggplot(snow_chela, aes(SIZE_5MM, fill = as.factor(MATURE))) +
  geom_histogram(position = "identity", binwidth = 5, alpha = 0.5, color = "darkgrey") +
  #ggtitle("5MM bins")+
  scale_fill_manual(values = c("darkturquoise", "salmon"), labels = c("Immature", "Mature"), name = "")+
  theme_bw()+
 # facet_wrap(~YEAR, nrow = 2)+
  xlab("Carapace width (mm)")+
  ylab("Chela measurements")+
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 12)) -> fig.2


spat.dat <- snow_chela %>%
              group_by(YEAR, LATITUDE, LONGITUDE) %>%
              reframe(N = n(),
                      MATURE = sum(MATURE),
                      PROP_MATURE = MATURE/N) 
              
ggplot()+
  geom_point(spat.dat %>% dplyr::filter(YEAR != 2012), mapping = aes(LONGITUDE, LATITUDE, fill = PROP_MATURE), stroke = NA, pch=21, color = "black", size = 1.5)+
  facet_wrap(~YEAR)+
  theme_bw()+
  #ggtitle("Snow")+
  scale_fill_viridis_c(option = "mako", direction = -1) +
  theme(legend.position = "bottom",
        legend.direction = "horizontal") -> s.spmat

ggsave(plot = s.spmat, "./Maturity data processing/Doc/snow_spatial_pmat.png", width = 8.5, height = 9)


# TANNER ----
# Load minima data, calculate cutline params
tanner_minima <- read.csv("./Maturity data processing/Output/chela_cutline_minima.csv") %>%
  filter(SPECIES == "TANNER") %>%
  mutate(BETA0 = coef(lm(MINIMUM ~ MIDPOINT))[1],
         BETA1 = coef(lm(MINIMUM ~ MIDPOINT))[2])

BETA0 <- unique(tanner_minima$BETA0)
BETA1 <- unique(tanner_minima$BETA1)

# Load tanner specimen data 
tanner_specimen<-  readRDS("./Maturity data processing/Data/tanner_survey_specimenEBS.rda")$specimen %>%
  filter(SHELL_CONDITION == 2, SEX == 1) %>%
  mutate(SIZE_1MM = floor(SIZE),
         BIN = cut_width(SIZE, width = 5, center = 2.5, closed = "left", dig.lab = 4),
         BIN2 = BIN) %>%
  separate(BIN2, sep = ",", into = c("LOWER", "UPPER")) %>%
  mutate(LOWER = as.numeric(sub('.', '', LOWER)),
         UPPER = as.numeric(gsub('.$', '', UPPER)),
         SIZE_BINNED = (UPPER + LOWER)/2) %>%
  rename(SIZE_5MM = SIZE_BINNED) %>%
  dplyr::select(YEAR, STATION_ID, CHELA_HEIGHT,LATITUDE, LONGITUDE, SIZE_5MM,
                SAMPLING_FACTOR) 

# Load snow chela data from chela database (subsample 2)
tanner_chela <-  read.csv("./Maturity data processing/Data/snow_tanner_cheladatabase.csv") %>% #already filtered appropriately
  dplyr::select(!X) %>%
  filter(SPECIES == "TANNER") %>%
  mutate(CUTOFF = BETA0 + BETA1*LN_CW, # apply cutline model
         MATURE = case_when((LN_CH > CUTOFF) ~ 1,
                            TRUE ~ 0),
         SIZE_1MM = floor(SIZE),
         BIN_5MM = cut_width(SIZE_1MM, width = 5, center = 2.5, closed = "left", dig.lab = 4),
         BIN2 = BIN_5MM, 
        MATURE = case_when((SIZE <=55) ~ 0,
                           (SIZE >=145) ~ 1,
                     TRUE ~ MATURE)) %>%
  separate(BIN2, sep = ",", into = c("LOWER", "UPPER")) %>%
  mutate(LOWER = as.numeric(sub('.', '', LOWER)),
         UPPER = as.numeric(gsub('.$', '', UPPER)),
         SIZE_5MM = (UPPER + LOWER)/2) %>%
  dplyr::select(!c(UPPER, LOWER)) %>%
  mutate(BIN_10MM = cut_width(SIZE_1MM, width = 10, center = 5, closed = "left", dig.lab = 4),
         BIN2 = BIN_10MM) %>%
  separate(BIN2, sep = ",", into = c("LOWER", "UPPER")) %>%
  mutate(LOWER = as.numeric(sub('.', '', LOWER)),
         UPPER = as.numeric(gsub('.$', '', UPPER)),
         SIZE_10MM = (UPPER + LOWER)/2) %>%
  dplyr::select(!c(UPPER, LOWER)) %>%
  mutate(BIN_15MM = cut_width(SIZE_1MM, width = 15, center = 7.5, closed = "left", dig.lab = 4),
         BIN2 = BIN_15MM) %>%
  separate(BIN2, sep = ",", into = c("LOWER", "UPPER")) %>%
  mutate(LOWER = as.numeric(sub('.', '', LOWER)),
         UPPER = as.numeric(gsub('.$', '', UPPER)),
         SIZE_15MM = (UPPER + LOWER)/2) %>%
  mutate(BIN_20MM = cut_width(SIZE_1MM, width = 20, center = 10, closed = "left", dig.lab = 4),
         BIN2 = BIN_20MM) %>%
  separate(BIN2, sep = ",", into = c("LOWER", "UPPER")) %>%
  mutate(LOWER = as.numeric(sub('.', '', LOWER)),
         UPPER = as.numeric(gsub('.$', '', UPPER)),
         SIZE_20MM = (UPPER + LOWER)/2) %>%
  dplyr::select(!c(UPPER, LOWER)) %>%
  group_by(YEAR, SIZE_5MM) %>%
  mutate(N_chela_5mm = n()) %>%
  ungroup() %>%
  group_by(YEAR, SIZE_10MM) %>%
  mutate(N_chela_10mm = n()) %>%
  ungroup() %>%
  group_by(YEAR, SIZE_15MM) %>%
  mutate(N_chela_15mm = n()) %>%
  ungroup() %>%
  group_by(YEAR, SIZE_20MM) %>%
  mutate(N_chela_20mm = n()) %>%
  ungroup()


ggplot(tanner_chela, aes(YEAR, SIZE_5MM, fill = N_chela_5mm))+
  geom_tile()+
  theme_bw()+
  ggtitle("Tanner")+
  scale_fill_viridis_c(option = "inferno") -> t.heatmap

s.heatmap/t.heatmap  + plot_layout(axes = "collect")

ggsave("./Maturity data processing/Figures/nchela_heatmap.png", height = 8, width =7)


ggplot(tanner_chela, aes(SIZE_5MM, fill = as.factor(MATURE))) +
  geom_histogram(position = "identity", binwidth = 5, alpha = 0.5) +
  ggtitle("5MM bins")+
  scale_fill_manual(values = c("darkturquoise", "salmon"), labels = c("Immature", "Mature"), name = "")+
  theme_bw()+
  facet_wrap(~YEAR, scales = "free_y")+
  ggtitle("Tanner")+
  xlab("Carapace width (mm)")+
  ylab("Chela measurements")+
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 12))


ggplot(tanner_chela, aes(SIZE_5MM, fill = as.factor(MATURE))) +
  geom_histogram(position = "identity", binwidth = 5, alpha = 0.5, color = "darkgrey") +
  ggtitle("5MM bins")+
  scale_fill_manual(values = c("darkturquoise", "salmon"), labels = c("Immature", "Mature"), name = "")+
  theme_bw()+
  ggtitle("Tanner")+
  xlab("Carapace width (mm)")+
  ylab("Chela measurements")+
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 12))





spat.dat <- tanner_chela %>%
  group_by(YEAR, LATITUDE, LONGITUDE) %>%
  reframe(N = n(),
          MATURE = sum(MATURE),
          PROP_MATURE = MATURE/N)

ggplot()+
  geom_point(spat.dat, mapping = aes(LONGITUDE, LATITUDE, fill = PROP_MATURE), stroke = NA, pch=21, color = "black", size = 1.5)+
  facet_wrap(~YEAR)+
  #ggtitle("Tanner")+
  theme_bw()+
  scale_fill_viridis_c(option = "mako", direction = -1) +
  theme(legend.position = "bottom",
        legend.direction = "horizontal") -> t.spmat

ggsave(plot = t.spmat, "./Maturity data processing/Doc/tanner_spatial_pmat.png", width = 8.5, height = 9)


s.spmat/t.spmat + plot_layout(axes = "collect", guides= "collect") & theme(legend.position = "bottom",
                                                                           legend.direction = "horizontal")

ggsave("./Maturity data processing/Doc/spatial_pmat.png", height = 11, width =8.5)



spat.spec <- tanner_specimen %>%
  #filter(is.na(CHELA_HEIGHT) == "FALSE") %>%
  group_by(YEAR, LATITUDE, LONGITUDE) %>%
  reframe(N = n()) %>%
  filter(YEAR %in% unique(tanner_chela$YEAR))

ggplot()+
  geom_point(spat.spec, mapping = aes(LONGITUDE, LATITUDE), size = 1)+
  facet_wrap(~YEAR)+
  theme_bw()




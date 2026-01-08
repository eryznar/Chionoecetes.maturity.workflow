# PURPOSE: to make a MWK ("maturity-width-key"; similar to age-length-key for fish) for Chionoecetes and apply
# that to CPUE - these are "designed" based indices

# Author: Emily Ryznar

# NOTES:


# LOAD LIBS/PARAMS ---------------------------------------------------------------------------------------
source("./Maturity data processing/Scripts/1) load_libs_params.R")

# LOAD DATA ----------------------------------------------------------------------------------------------
# Load minima data, calculate cutline params
snow_minima <- read.csv("./Maturity data processing/Output/chela_cutline_minima.csv") %>%
  filter(SPECIES == "SNOW") %>%
  mutate(BETA0 = coef(lm(MINIMUM ~ MIDPOINT))[1],
         BETA1 = coef(lm(MINIMUM ~ MIDPOINT))[2])

BETA0 <- unique(snow_minima$BETA0)
BETA1 <- unique(snow_minima$BETA1)

# Load snow specimen data (subsample 1)
snow_specimen<-  readRDS("./Maturity data processing/Data/snow_survey_specimenEBS.rda")$specimen %>%
                    filter(SHELL_CONDITION == 2, SEX == 1) %>%
                    mutate(SIZE_1MM = floor(SIZE),
                           BIN_5MM = cut_width(SIZE_1MM, width = 5, center = 2.5, closed = "left", dig.lab = 4),
                           BIN2 = BIN_5MM) %>%
                    separate(BIN2, sep = ",", into = c("LOWER", "UPPER")) %>%
                    mutate(LOWER = as.numeric(sub('.', '', LOWER)),
                           UPPER = as.numeric(gsub('.$', '', UPPER)),
                           SIZE_5MM = (UPPER + LOWER)/2) %>%
                    dplyr::select(YEAR, STATION_ID, LATITUDE, LONGITUDE, SIZE_1MM, SIZE_5MM, 
                                  SAMPLING_FACTOR, AREA_SWEPT, CALCULATED_WEIGHT_1MM)

# Load snow chela data from chela database (subsample 2)
snow_chela <-  read.csv("./Maturity data processing/Data/snow_tanner_cheladatabase.csv") %>% #already filtered appropriately
                    filter(SPECIES == "SNOW") %>%
                    mutate(CUTOFF = BETA0 + BETA1*LN_CW, # apply cutline model
                           MATURE = case_when((LN_CH > CUTOFF) ~ 1,
                                              TRUE ~ 0),
                           SIZE_1MM = floor(SIZE),
                           BIN_5MM = cut_width(SIZE_1MM, width = 5, center = 2.5, closed = "left", dig.lab = 4),
                           BIN2 = BIN_5MM) %>%
                    separate(BIN2, sep = ",", into = c("LOWER", "UPPER")) %>%
                    mutate(LOWER = as.numeric(sub('.', '', LOWER)),
                           UPPER = as.numeric(gsub('.$', '', UPPER)),
                           SIZE_5MM = (UPPER + LOWER)/2,
                           MATURE = case_when((SIZE <40) ~ 0,
                                              (SIZE >= 135) ~ 1,
                                              TRUE ~ MATURE)) %>%
                    dplyr::select(YEAR, STATION_ID, LATITUDE, LONGITUDE, SIZE_1MM, SIZE_5MM, 
                                  LN_CH, LN_CW, MATURE, AREA_SWEPT)
 

# MAKE MWK ----------------------------------------------------------------------------------------
##   Aggregate sub 2 (chela data) information: the `p_ysw` df will contain a field
##   msrd_crab that is the total number of chela crab collected in year-y, at station-s, and at
##   width-w (SIZE).
##   Calculate distribution of maturity proportions (prop_mat) for a given width and station, `p_ysw`.
##   This is the non-global maturity-width key.

p_ysw <- snow_chela %>%
          group_by(YEAR, SIZE_5MM) %>%
          reframe(MSRD_CRAB = n(),
                  TOTAL_MATURE = sum(MATURE),
                  PROP_MATURE = TOTAL_MATURE/MSRD_CRAB) %>%
          mutate(GLOBAL_FILLED = "N") # these have not been filled by a global mwk
          
## Calculate every combination of width (from sub 1) for global MWK (## i.e., an ALK that is aggregated over years
## so we can fill in years where a width-bin is missing)
p_ysw <- right_join(p_ysw, snow_specimen)

missing <- p_ysw %>% filter(is.na(TOTAL_MATURE) == TRUE)

p_sw <- snow_chela %>%
          group_by(SIZE_5MM) %>%
          reframe(MSRD_CRAB = n(),
          TOTAL_MATURE = sum(MATURE),
          PROP_MATURE = TOTAL_MATURE/MSRD_CRAB)

## For each record in `missing`, merge all the mature probabilities
## from `p_sw` for that width bin.
missing2 <- right_join(p_sw, missing %>% dplyr::select(!c(MSRD_CRAB, TOTAL_MATURE, PROP_MATURE))) %>%
              mutate(PROP_MATURE = case_when(((is.na(PROP_MATURE) == TRUE) & SIZE_1MM < 25) ~ 0, # filling in the rest of the blanks
                                             ((is.na(PROP_MATURE) == TRUE) & SIZE_1MM > 130) ~ 1,
                                             TRUE ~ PROP_MATURE),
                     GLOBAL_FILLED = "Y") # these have been filled by a global mwk

missing2 %>% filter(is.na(PROP_MATURE) == TRUE) -> ll

## Append the globally-filled widths with the the non-global `p_ysw` mwk
## to get a now global mwk. 

p_ysw2 <- rbind(p_ysw %>% filter(is.na(PROP_MATURE) == FALSE), missing2) 


## Clean and save
p_ysw2 %>% dplyr::select(YEAR, STATION_ID, LATITUDE, LONGITUDE, SIZE_1MM, SIZE_5MM, PROP_MATURE, GLOBAL_FILLED) %>%
  distinct() %>%
  write.csv(., "./Maturity research/Output/MWK.csv")

## Plot
plot.dat <- p_ysw2 %>% group_by(YEAR, SIZE_5MM) %>% reframe(PROP_MATURE = mean(PROP_MATURE))

ggplot(plot.dat, aes(SIZE_5MM, PROP_MATURE))+
  geom_line(linewidth = 1)+
  facet_wrap(~YEAR)+
  theme_bw()

# Plot chela sample size by maturity and bin
snow_chela %>%
  rename("1MM" = "SIZE_1MM", "5MM" = "SIZE_5MM") %>%
  dplyr::select(!c(STATION_ID, LATITUDE, LONGITUDE)) %>%
  pivot_longer(., c(2:3), values_to = "size", names_to = "bin") %>%
  mutate(MATURE = case_when((MATURE == 0) ~ "Immature",
                            TRUE ~ "Mature")) -> plot.dat2

ggplot(plot.dat2 %>% filter(bin == "1MM"), aes(size, fill = MATURE)) +
  geom_histogram(position = "identity", binwidth = 1, alpha = 0.5, color = "darkgrey") +
  ggtitle("1MM bins")+
  scale_fill_manual(values = c("darkturquoise", "salmon"), labels = c("Immature", "Mature"), name = "")+
  theme_bw()+
  xlab("Carapace width (mm)")+
  ylab("Chela measurements")+
  theme(legend.position = "bottom",
        legend.direction = "horizontal")

ggplot(plot.dat2 %>% filter(bin == "5MM"), aes(size, fill = MATURE)) +
  geom_histogram(position = "identity", binwidth = 5, alpha = 0.5, color = "darkgrey") +
  ggtitle("5MM bins")+
  scale_fill_manual(values = c("darkturquoise", "salmon"), labels = c("Immature", "Mature"), name = "")+
  theme_bw()+
  xlab("Carapace width (mm)")+
  ylab("Chela measurements")+
  theme(legend.position = "bottom",
        legend.direction = "horizontal")


ggplot(plot.dat2 %>% filter(bin == "5MM"), aes(size, fill = MATURE)) +
  geom_histogram(position = "identity", binwidth = 5, alpha = 0.5) +
  facet_wrap(~YEAR, scales = "free_y") +
  ggtitle("5MM bins")+
  scale_fill_manual(values = c("darkturquoise", "salmon"), labels = c("Immature", "Mature"), name = "")+
  theme_bw()+
  xlab("Carapace width (mm)")+
  ylab("Chela measurements")+
  theme(legend.position = "bottom",
        legend.direction = "horizontal")


ggplot(plot.dat2 %>% filter(bin == "1MM"), aes(size, fill = MATURE)) +
  geom_histogram(position = "identity", binwidth = 1, alpha = 0.5) +
  facet_wrap(~YEAR, scales = "free_y") +
  ggtitle("1MM bins")+
  scale_fill_manual(values = c("darkturquoise", "salmon"), labels = c("Immature", "Mature"), name = "")+
  theme_bw()+
  xlab("Carapace width (mm)")+
  ylab("Chela measurements")+
  theme(legend.position = "bottom",
        legend.direction = "horizontal")

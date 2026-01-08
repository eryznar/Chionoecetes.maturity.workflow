# PURPOSE: to fit binomial GAM to mature (0,1) weighted by CPU, log CPUE, or unweighted

# Author: Emily Ryznar

# NOTES:
# Make sure you're using updated chela database compiled by shannon, though can't use it for everything bc
# it doesn't have sampling factor...

# LOAD LIBS/PARAMS ---------------------------------------------------------------------------------------
source("./Scripts/load_libs_params.R")

# SNOW CRAB ----
# Get survey cpue by 1mm bin for ALL shell 2 males, not just chela msrd
snow_cpue <- calc_cpue(crab_data = readRDS("./Data/snow_survey_specimenEBS.rda"),
                       species = "SNOW",
                       years = years,
                       sex = "male",
                       shell_condition = "new_hardshell",
                       bin_1mm = TRUE) %>%
  dplyr::select(YEAR, STATION_ID, LATITUDE, LONGITUDE, SIZE_1MM, CPUE)


# Join 1mm cpue to chela crab
mod.dat <- snow_chela %>%
  mutate(SIZE_1MM = floor(SIZE),
         MAT_TEXT = case_when((MATURE == 1) ~ "Mature",
                              TRUE ~ "Immature")) %>% # could change this to 10mm or 5mm bins instead of 1mm
  right_join(., snow_cpue) %>% # adding in haul-level CPUE for corresponding 1mm size bin
  mutate(YEAR = as.factor(YEAR)) %>%
  mutate(MATURE = case_when((SIZE <40) ~ 0,
                            (SIZE >= 115) ~ 1,
                            TRUE ~ MATURE),
         log.CPUE = as.integer(round(log(CPUE+10))),
         CPUE = as.integer(round(CPUE))) %>%
  na.omit() 

# Binomial GAM for loop to fit and predict by year
yrs <- unique(mod.dat$YEAR)
plot.dat <- data.frame()

for(ii in 1:length(yrs)){
  mod.dat %>%
    filter(YEAR %in% yrs[ii]) -> mod.dat2
  
  print(paste0("Fitting ", yrs[ii]))
  
  new.dat <- data.frame(YEAR = unique(mod.dat2$YEAR), SIZE = mod.dat2$SIZE, CPUE = mod.dat2$CPUE)
  
  mod.1 <- bam(MATURE ~ s(SIZE), family = binomial(link = "logit"), data = mod.dat2)
  prd.1 <- predict(mod.1, new.dat, type = "response", se.fit = TRUE)
  
  mod.2 <- bam(MATURE ~ s(SIZE), family = binomial(link = "logit"), data = mod.dat2, weights = CPUE)
  prd.2 <- predict(mod.2, new.dat, type = "response", se.fit = TRUE)
  
  
  mod.3 <- bam(MATURE ~ s(SIZE), family = binomial(link = "logit"), data = mod.dat2, weights = log.CPUE)
  prd.3 <- predict(mod.3, new.dat, type = "response", se.fit = TRUE)
  
  
  out <- rbind(data.frame(new.dat, prd = prd.1, type = "Unweighted"),
               data.frame(new.dat, prd = prd.2, type = "CPUE-weighted"),
               data.frame(new.dat, prd = prd.3, type = "log(CPUE)-weighted")) %>%
    mutate(prd.fit = case_when((SIZE <=44) ~ 0,
                               (SIZE >= 115) ~ 1,
                               TRUE ~ prd.fit),
           N = nrow(mod.dat2))
  
  plot.dat <- rbind(plot.dat, out)
}

labelz <- paste0(unique(plot.dat$YEAR), " \n(N=", unique(plot.dat$N), ")")
names(labelz) <- unique(plot.dat$YEAR)

ggplot()+
  geom_line(plot.dat %>% filter(type == "Unweighted"), mapping = aes(SIZE, prd.fit, color = as.factor(1), group = YEAR), linewidth = 1)+
  geom_line(plot.dat %>% filter(type == "CPUE-weighted"), mapping = aes(SIZE, prd.fit, color = as.factor(2), group = YEAR), linewidth = 1)+
  geom_line(plot.dat %>% filter(type == "log(CPUE)-weighted"), mapping = aes(SIZE, prd.fit, color = as.factor(3), group = YEAR), linewidth = 1)+
  scale_color_manual(values = c("darkgoldenrod", "cadetblue", "darksalmon"), labels = c("Unweighted", "CPUE-weighted", "log(CPUE)-weighted"), name = "")+
  scale_fill_manual(values = c("darkgoldenrod", "cadetblue", "darksalmon"), labels = c("Unweighted", "CPUE-weighted", "log(CPUE)-weighted"), name = "")+
  theme_bw()+
  facet_wrap(~YEAR, labeller = labeller(YEAR = labelz))+
  # geom_ribbon(plot.dat %>% filter(type == "Unweighted"),
  #           mapping = aes(SIZE, ymin = prd - prd.se, ymax = prd + prd.se, fill = as.factor(1)), alpha = 0.25)+
  # geom_ribbon(plot.dat %>% filter(type == "CPUE-weighted"),
  #             mapping = aes(SIZE, ymin = prd - prd.se, ymax = prd + prd.se, fill = as.factor(2)), alpha = 0.25)+
  ylab("P(molt)")+
  xlab("Carapace width (mm)")+
  ggtitle("Snow crab")+
  scale_x_continuous(limits = c(25, 138), breaks = seq(25, 140, by = 25))+
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 13),
        legend.text = element_text(size = 13),
        strip.text = element_text(size = 11), 
        legend.position = "bottom",
        legend.direction = "horizontal")

ggsave("./Figures/snow_compareweights_binomialGAM.png", width = 10, height = 9, units = "in")


ggplot(plot.dat %>% filter(type == "Unweighted", SIZE >40 & SIZE < 95), aes(CPUE, prd.fit))+
  geom_point(size = 1)+
  ylab("p(molt)")+
  ggtitle("Size 40-95mm")+
  facet_wrap(~YEAR)+
  ylim(0, 1)+
  #geom_smooth(se = FALSE)+
  theme_bw()

# TANNER CRAB ----
# Get survey cpue by 1mm bin for ALL shell 2 males, not just chela msrd
tanner_cpue <- calc_cpue(crab_data = readRDS("./Data/tanner_survey_specimenEBS.rda"),
                       species = "TANNER",
                       years = years,
                       sex = "male",
                       shell_condition = "new_hardshell",
                       bin_1mm = TRUE) %>%
  dplyr::select(YEAR, STATION_ID, LATITUDE, LONGITUDE, SIZE_1MM, CPUE)


# Join 1mm cpue to chela crab
mod.dat <- tanner_chela %>%
  filter(SPECIES == "TANNER") %>%
  mutate(SIZE_1MM = floor(SIZE),
         MAT_TEXT = case_when((MATURE == 1) ~ "Mature",
                              TRUE ~ "Immature")) %>% # could change this to 10mm or 5mm bins instead of 1mm
  right_join(., tanner_cpue) %>% # adding in haul-level CPUE for corresponding 1mm size bin
  mutate(YEAR = as.factor(YEAR)) %>%
  mutate(MATURE = case_when((SIZE <=55) ~ 0,
                            (SIZE >= 145) ~ 1,
                            TRUE ~ MATURE),
         log.CPUE = as.integer(round(log(CPUE+10))),
         CPUE = as.integer(round(CPUE))) %>%
  na.omit() 

# Binomial GAM for loop to fit and predict by year
yrs <- unique(mod.dat$YEAR)
plot.dat <- data.frame()

for(ii in 1:length(yrs)){
  mod.dat %>%
    filter(YEAR %in% yrs[ii]) -> mod.dat2
  
  print(paste0("Fitting ", yrs[ii]))
  
  new.dat <- data.frame(YEAR = unique(mod.dat2$YEAR), SIZE = seq(min(mod.dat2$SIZE), max(mod.dat2$SIZE), by = 1))
  
  mod.1 <- bam(MATURE ~ s(SIZE), family = binomial(link = "logit"), data = mod.dat2)
  prd.1 <- predict(mod.1, new.dat, type = "response", se.fit = TRUE)
  
  mod.2 <- bam(MATURE ~ s(SIZE), family = binomial(link = "logit"), data = mod.dat2, weights = CPUE)
  prd.2 <- predict(mod.2, new.dat, type = "response", se.fit = TRUE)
  
  
  mod.3 <- bam(MATURE ~ s(SIZE), family = binomial(link = "logit"), data = mod.dat2, weights = log.CPUE)
  prd.3 <- predict(mod.3, new.dat, type = "response", se.fit = TRUE)
  
  
  out <- rbind(data.frame(new.dat, prd = prd.1, type = "Unweighted"),
               data.frame(new.dat, prd = prd.2, type = "CPUE-weighted"),
               data.frame(new.dat, prd = prd.3, type = "log(CPUE)-weighted")) %>%
    mutate(prd.fit = case_when((SIZE <=55) ~ 0,
                               (SIZE >= 145) ~ 1,
                               TRUE ~ prd.fit),
           N = nrow(mod.dat2))
  
  plot.dat <- rbind(plot.dat, out)
}

labelz <- paste0(unique(plot.dat$YEAR), " \n(N=", unique(plot.dat$N), ")")
names(labelz) <- unique(plot.dat$YEAR)

ggplot()+
  geom_line(plot.dat %>% filter(type == "Unweighted"), mapping = aes(SIZE, prd.fit, color = as.factor(1), group = YEAR), linewidth = 1)+
  geom_line(plot.dat %>% filter(type == "CPUE-weighted"), mapping = aes(SIZE, prd.fit, color = as.factor(2), group = YEAR), linewidth = 1)+
  geom_line(plot.dat %>% filter(type == "log(CPUE)-weighted"), mapping = aes(SIZE, prd.fit, color = as.factor(3), group = YEAR), linewidth = 1)+
  scale_color_manual(values = c("darkgoldenrod", "cadetblue", "darksalmon"), labels = c("Unweighted", "CPUE-weighted", "log(CPUE)-weighted"), name = "")+
  scale_fill_manual(values = c("darkgoldenrod", "cadetblue", "darksalmon"), labels = c("Unweighted", "CPUE-weighted", "log(CPUE)-weighted"), name = "")+
  theme_bw()+
  facet_wrap(~YEAR, labeller = labeller(YEAR = labelz))+
  # geom_ribbon(plot.dat %>% filter(type == "Unweighted"),
  #           mapping = aes(SIZE, ymin = prd - prd.se, ymax = prd + prd.se, fill = as.factor(1)), alpha = 0.25)+
  # geom_ribbon(plot.dat %>% filter(type == "CPUE-weighted"),
  #             mapping = aes(SIZE, ymin = prd - prd.se, ymax = prd + prd.se, fill = as.factor(2)), alpha = 0.25)+
  ylab("P(molt)")+
  xlab("Carapace width (mm)")+
  ggtitle("Tanner crab")+
  #scale_x_continuous(limits = c(25, 138), breaks = seq(25, 140, by = 25))+
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 13),
        legend.text = element_text(size = 13),
        strip.text = element_text(size = 11), 
        legend.position = "bottom",
        legend.direction = "horizontal")

ggsave("./Figures/tanner_compareweights_binomialGAM.png", width = 10, height = 9, units = "in")

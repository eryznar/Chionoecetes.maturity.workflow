# SNOW CRAB
snow_ogives <- read.csv("./Maturity data processing/Data/snow_legacy_matmale_ratio.csv") 

snow.ogive.spec <-  readRDS("./Maturity data processing/Data/snow_survey_specimenEBS.rda")$specimen %>%
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
  filter(YEAR %in% unique(snow_mod$data$YEAR))

snow_ogives
yy <- unique(snow_ogives$YEAR)
preds <- data.frame()
for(ii in 1:length(yy)){
  m.dat <- snow_ogives %>% filter(YEAR == yy[ii]) %>% mutate(NUM_MATURE = round(NUM_MATURE), NUM_IMMATURE = round(NUM_IMMATURE))
  n.dat <- snow.ogive.spec %>% filter(YEAR == yy[ii]) %>% rename(SIZE_BIN = SIZE_5MM) # using same 5mm size bins in spec dat
  
  mod <- gam(
    cbind(NUM_MATURE, NUM_IMMATURE) ~ s(SIZE_BIN),
    family = binomial(link = "logit"),
    data   = m.dat
  )
  
  p.dat <- data.frame(n.dat, PROP_MATURE = predict(mod, newdata= n.dat, type = "response"))
  
  preds <- rbind(preds, p.dat)
}

preds2 <- preds %>%
  group_by(YEAR, SIZE_BIN) %>%
  reframe(PROP_MATURE = mean(PROP_MATURE)) %>%
            mutate(PROP_MATURE = case_when(SIZE_BIN <=40 ~ 0,
                                           SIZE_BIN >=115 & PROP_MATURE <0.95 ~ 1, # enforcing good behavior
                                           TRUE ~ PROP_MATURE))

write.csv(preds2, "./Maturity data processing/Doc/snow_pmat_5mminterp.csv")


ggplot(preds2, aes(SIZE_BIN, PROP_MATURE))+
  geom_line(linewidth = 1)+
  theme_bw()+
    facet_wrap(~YEAR)

ggplot(preds2, aes(SIZE_BIN, PROP_MATURE, color = YEAR, group = YEAR))+
  geom_line(linewidth = 1)+
  theme_bw()


# TANNER E ----
tannerE_ogives <- read.csv("./Maturity data processing/Data/tanner_legacy_matmale_ratio.csv") %>%
  filter(DISTRICT == "E166")

tannerE.ogive.spec <-  readRDS("./Maturity data processing/Data/tanner_survey_specimenEBS.rda")$specimen %>%
  filter(SHELL_CONDITION == 2, SEX == 1, DISTRICT == "E166") %>%
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
  filter(YEAR %in% unique(tanner_mod$data$YEAR))


yy <- unique(tannerE_ogives$YEAR)
preds <- data.frame()
for(ii in 1:length(yy)){
  m.dat <- tannerE_ogives %>% filter(YEAR == yy[ii]) %>% mutate(NUM_MATURE = round(NUM_MATURE), NUM_IMMATURE = round(NUM_IMMATURE))
  n.dat <- tannerE.ogive.spec %>% filter(YEAR == yy[ii]) %>% rename(SIZE_BIN = SIZE_5MM) # using same 5mm size bins in spec dat
  
  mod <- gam(
    cbind(NUM_MATURE, NUM_IMMATURE) ~ s(SIZE_BIN),
    family = binomial(link = "logit"),
    data   = m.dat
  )
  
  #n.dat = data.frame(YEAR = yy[ii], SIZE_BIN = seq(0, 200, by = 5))
  p.dat <- data.frame(n.dat, PROP_MATURE = predict(mod, newdata= n.dat, type = "response"))
  
  preds <- rbind(preds, p.dat)
}

tanE_preds <- preds %>%
  group_by(YEAR, SIZE_BIN) %>%
  reframe(PROP_MATURE = mean(PROP_MATURE)) %>%
  #filter(YEAR != 1994) %>%
  mutate(PROP_MATURE = case_when(SIZE_BIN <=65 & PROP_MATURE > 0.05 ~ 0, # enforcing good behavior
                                 TRUE ~ PROP_MATURE))

write.csv(tanE_preds, "./Maturity data processing/Doc/tanE_pmat_5mminterp.csv")


ggplot(tanE_preds, aes(SIZE_BIN, PROP_MATURE))+
  geom_line(linewidth = 1)+
  theme_bw()+
  facet_wrap(~YEAR)

ggplot(tanE_preds, aes(SIZE_BIN, PROP_MATURE, color = YEAR, group = YEAR))+
  geom_line(linewidth = 1)+
  theme_bw()

# TANNER W ----
tannerW_ogives <- read.csv("./Maturity data processing/Data/tanner_legacy_matmale_ratio.csv") %>%
  filter(DISTRICT == "W166")


tanner_chela <-  read.csv("./Maturity data processing/Data/snow_tanner_cheladatabase.csv") %>% #already filtered appropriately
  dplyr::select(!X) %>%
  filter(SPECIES == "TANNER") %>%
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
  #filter(N_chela > 40) %>% #necessary for consecutive size correlations to run
  mutate(LATITUDE = Y/1000, # scale to km so values don't get too large
         LONGITUDE = X/1000,
         SIZE_CATEGORY = as.factor(paste0("SIZE", SIZE_BINNED)),
         YEAR_F = as.factor(YEAR),
         YEAR_SCALED = as.numeric(scale(YEAR)),
         MATURE = case_when((SIZE <=55) ~ 0,
                            (SIZE >= 145) ~ 1,
                            TRUE ~ MATURE)) %>%
  as.data.frame(.) %>%
  rename(SIZE_5MM = SIZE_BINNED) %>%
  dplyr::select(YEAR, YEAR_F, YEAR_SCALED, STATION_ID, LATITUDE, LONGITUDE, SIZE_5MM, SIZE_CATEGORY, MATURE)



tannerW.ogive.spec <-  readRDS("./Maturity data processing/Data/tanner_survey_specimenEBS.rda")$specimen %>%
  filter(SHELL_CONDITION == 2, SEX == 1, DISTRICT == "W166") %>%
  mutate(SIZE_1MM = floor(SIZE),
         BIN = cut_width(SIZE, width = 5, center = 2.5, closed = "left", dig.lab = 4),
         BIN2 = BIN, 
         LN_CW = log(SIZE),
         LN_CH = log(CHELA_HEIGHT)) %>%
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
  filter(YEAR %in% unique(tanner_mod$data$YEAR))


yy <- unique(tannerW_ogives$YEAR)
preds <- data.frame()
for(ii in 1:length(yy)){
  m.dat <- tannerW_ogives %>% filter(YEAR == yy[ii]) %>% mutate(NUM_MATURE = round(NUM_MATURE), NUM_IMMATURE = round(NUM_IMMATURE))
  n.dat <- tannerW.ogive.spec %>% filter(YEAR == yy[ii]) %>% rename(SIZE_BIN = SIZE_5MM) # using same 5mm size bins in spec dat
  
  mod <- gam(
    cbind(NUM_MATURE, NUM_IMMATURE) ~ s(SIZE_BIN),
    family = binomial(link = "logit"),
    data   = m.dat
  )
  
  #n.dat = data.frame(YEAR = yy[ii], SIZE_BIN = seq(0, 200, by = 5))
  p.dat <- data.frame(n.dat, PROP_MATURE = predict(mod, newdata= n.dat, type = "response"))
  
  preds <- rbind(preds, p.dat)
}

tanW_preds <- preds %>%
  group_by(YEAR, SIZE_BIN) %>%
  reframe(PROP_MATURE = mean(PROP_MATURE)) %>%
  #filter(YEAR != 1994) %>%
  mutate(PROP_MATURE = case_when(SIZE_BIN <=65 & PROP_MATURE > 0.05 ~ 0, # enforcing good behavior
                                 SIZE_BIN >=140 & PROP_MATURE <0.75 ~ 1,
                                 TRUE ~ PROP_MATURE))

write.csv(tanW_preds, "./Maturity data processing/Doc/tanW_pmat_5mminterp.csv")


# ggplot(tanW_preds, aes(SIZE_BIN, PROP_MATURE))+
#   geom_line(linewidth = 1)+
#   theme_bw()+
#   facet_wrap(~YEAR)

ggplot(tanW_preds, aes(SIZE_BIN, PROP_MATURE, color = YEAR, group = YEAR))+
  geom_line(linewidth = 1)+
  theme_bw()



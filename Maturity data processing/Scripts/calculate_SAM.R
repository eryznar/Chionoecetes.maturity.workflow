crab_dat <- readRDS("./Maturity data processing/Data/snow_survey_specimenEBS.rda")

crab_dat$specimen <-  crab_dat$specimen %>%
  filter(YEAR %in% years, SHELL_CONDITION == 2, SEX == 1) %>%
  mutate(BIN_5MM = cut_width(SIZE_1MM, width = 5, center = 2.5, closed = "left", dig.lab = 4),
         BIN2 = BIN_5MM) %>%
  separate(BIN2, sep = ",", into = c("LOWER", "UPPER")) %>%
  mutate(LOWER = as.numeric(sub('.', '', LOWER)),
         UPPER = as.numeric(gsub('.$', '', UPPER)),
         SIZE_5MM = (UPPER + LOWER)/2) %>%
  dplyr::select(!c(BIN_5MM, LOWER, UPPER)) %>%
  mutate(BIN_10MM = cut_width(SIZE_1MM, width = 10, center = 5, closed = "left", dig.lab = 4),
         BIN2 = BIN_10MM) %>%
  separate(BIN2, sep = ",", into = c("LOWER", "UPPER")) %>%
  mutate(LOWER = as.numeric(sub('.', '', LOWER)),
         UPPER = as.numeric(gsub('.$', '', UPPER)),
         SIZE_10MM = (UPPER + LOWER)/2, 
         YEAR_F = as.factor(YEAR)) %>%
  dplyr::select(!c(BIN_10MM, LOWER, UPPER))

sdmTMB.spec <- crab_dat
mod <- readRDS(paste0(remote_dir, "/SNOW/sdmTMB/s(SIZE, k=13)_iid_200_sdmTMB.rda"))

ss <- readRDS("./Maturity research/Output/sdmTMB_specdat.csv") %>%
  mutate(SEL = predict(s.gam, newdata= ., type = "response"), # predict size-specific selectivity
         SAMPLING_FACTOR_SEL = SAMPLING_FACTOR/SEL, # account for size specific selectivity in abundance
         SAMPLING_FACTOR_MATURE = SAMPLING_FACTOR_SEL * PROP_MATURE,  
         SAMPLING_FACTOR_IMMATURE = SAMPLING_FACTOR_SEL - SAMPLING_FACTOR_MATURE)

# sdmTMB.spec$specimen %>%
#   filter(YEAR %in% years) %>%
#   st_as_sf(., coords = c("LONGITUDE", "LATITUDE"), crs = "+proj=longlat +datum=WGS84") %>%
#   st_transform(., crs = "+proj=utm +zone=2") %>%
#   cbind(st_coordinates(.)) %>%
#   as.data.frame(.) %>%
#   mutate(LATITUDE = Y/1000, # scale to km so values don't get too large
#          LONGITUDE = X/1000,
#          YEAR_F = as.factor(YEAR)) %>%
#   mutate(PROP_MATURE = predict(mod, newdata = ., type = "response")$est) -> ss

odat <- ss %>%
  mutate(CPUE = SAMPLING_FACTOR/AREA_SWEPT,
         CPUE_SEL = SAMPLING_FACTOR_SEL/AREA_SWEPT) %>%
  group_by(YEAR, SIZE_5MM) %>%
  reframe(
    PROP_MATURE_WTD   = sum(PROP_MATURE * SAMPLING_FACTOR) / sum(SAMPLING_FACTOR),
    PROP_MATURE_WTDSEL = sum(PROP_MATURE * SAMPLING_FACTOR_SEL) / sum(SAMPLING_FACTOR_SEL),
    PROP_MATURE_UNWTD = mean(PROP_MATURE)
  )

ggplot()+
  geom_line(odat, mapping = aes(SIZE_5MM, PROP_MATURE_WTD, color = as.factor(1)), linewidth = 1)+
  geom_line(odat, mapping = aes(SIZE_5MM, PROP_MATURE_UNWTD, color = as.factor(2)), linewidth = 1)+
  facet_wrap(~YEAR)+
  ylab("Probability of having undergone terminal molt")+
  theme_bw()+
  scale_color_manual(values = c("violet", "cadetblue"), labels = c("Weighted", "Unweighted"), name = "")

odatsel <- odat %>% dplyr::select(!PROP_MATURE_WTD) %>% rename(PROP_MATURE_WTD = PROP_MATURE_WTDSEL)
wtdSAM_by_year <- odatsel %>%
  arrange(YEAR, SIZE_5MM) %>%
  group_by(YEAR) %>%
  group_modify(~{
    df <- .x
    
    # first bin where weighted proportion ≥ 0.5
    i_upper <- which(df$PROP_MATURE_WTD >= 0.5)[1]
    i_lower <- i_upper - 1
    
    if (!is.na(i_upper) && i_upper > 1) {
      size_lower <- df$SIZE_5MM[i_lower]
      p_lower    <- df$PROP_MATURE_WTD[i_lower]
      size_upper <- df$SIZE_5MM[i_upper]
      p_upper    <- df$PROP_MATURE_WTD[i_upper]
      
      SAM <- size_lower + ((0.5 - p_lower) / (p_upper - p_lower)) *
        (size_upper - size_lower)
      
      tibble(SAM = SAM)
    } else {
      tibble(SAM = NA_real_)   # never reaches 0.5 or 1st bin ≥ 0.5
    }
  }) %>%
  ungroup()

unwtdSAM_by_year <- odatsel %>%
  arrange(YEAR, SIZE_5MM) %>%
  group_by(YEAR) %>%
  group_modify(~{
    df <- .x
    
    # first bin where weighted proportion ≥ 0.5
    i_upper <- which(df$PROP_MATURE_UNWTD >= 0.5)[1]
    i_lower <- i_upper - 1
    
    if (!is.na(i_upper) && i_upper > 1) {
      size_lower <- df$SIZE_5MM[i_lower]
      p_lower    <- df$PROP_MATURE_UNWTD[i_lower]
      size_upper <- df$SIZE_5MM[i_upper]
      p_upper    <- df$PROP_MATURE_UNWTD[i_upper]
      
      SAM <- size_lower + ((0.5 - p_lower) / (p_upper - p_lower)) *
        (size_upper - size_lower)
      
      tibble(SAM = SAM)
    } else {
      tibble(SAM = NA_real_)   # never reaches 0.5 or 1st bin ≥ 0.5
    }
  }) %>%
  ungroup()

SAM.dat <- rbind(wtdSAM_by_year %>% mutate(WEIGHT = "Y"),
                 unwtdSAM_by_year %>% mutate(WEIGHT = "N")) %>%
  right_join(., data.frame(YEAR = rep(seq(min(.$YEAR), max(.$YEAR), by = 1),2), WEIGHT = c("Y", "N")))

ggplot(SAM.dat, aes(YEAR, SAM, color = WEIGHT)) +
  geom_point(size = 2) +
  geom_line(linewidth = 1) +
  scale_color_manual(
    values = c("Y" = "violet", "N" = "cadetblue"),
    labels = c("Y" = "Weighted", "N" = "Unweighted"),
    name = ""
  ) +
  theme_bw()#+
  #geom_smooth(method = "lm")

ll <- lme(SAM ~ YEAR, data = SAM.dat %>% filter(WEIGHT == "Y") %>% na.omit(),  random = ~ 1 | YEAR, correlation = corAR1())

write.csv(SAM.dat, "./Maturity data processing/SAM_comparison.csv")

# Spatial SAM
odat_sp <- ss %>%
  group_by(YEAR, STATION_ID, LATITUDE, LONGITUDE, SIZE_5MM) %>%
  reframe(
    PROP_MATURE_WTD = sum(PROP_MATURE * SAMPLING_FACTOR) /
      sum(SAMPLING_FACTOR),
    PROP_MATURE_UNWTD = mean(PROP_MATURE),  # simple average
    N_CRAB = sum(SAMPLING_FACTOR),
    .groups = "drop"
  )

SAM_spatial <- odat_sp %>%
  arrange(YEAR, STATION_ID, LATITUDE, LONGITUDE, SIZE_5MM) %>%
  group_by(YEAR, STATION_ID, LATITUDE, LONGITUDE) %>%
  group_modify(~{
    df <- .x
    
    # weighted SAM
    i_upper_w <- which(df$PROP_MATURE_WTD >= 0.5)[1]
    i_lower_w <- i_upper_w - 1
    
    SAM_WTD <- NA_real_
    if (!is.na(i_upper_w) && i_upper_w > 1) {
      size_lower <- df$SIZE_5MM[i_lower_w]
      p_lower    <- df$PROP_MATURE_WTD[i_lower_w]
      size_upper <- df$SIZE_5MM[i_upper_w]
      p_upper    <- df$PROP_MATURE_WTD[i_upper_w]
      
      SAM_WTD <- size_lower + ((0.5 - p_lower) / (p_upper - p_lower)) *
        (size_upper - size_lower)
    }
    
    # unweighted SAM
    i_upper_u <- which(df$PROP_MATURE_UNWTD >= 0.5)[1]
    i_lower_u <- i_upper_u - 1
    
    SAM_UNWTD <- NA_real_
    if (!is.na(i_upper_u) && i_upper_u > 1) {
      size_lower <- df$SIZE_5MM[i_lower_u]
      p_lower    <- df$PROP_MATURE_UNWTD[i_lower_u]
      size_upper <- df$SIZE_5MM[i_upper_u]
      p_upper    <- df$PROP_MATURE_UNWTD[i_upper_u]
      
      SAM_UNWTD <- size_lower + ((0.5 - p_lower) / (p_upper - p_lower)) *
        (size_upper - size_lower)
    }
    
    tibble(
      SAM_WTD    = SAM_WTD,
      SAM_UNWTD  = SAM_UNWTD
    )
  }) %>%
  ungroup() %>%
  distinct()

write.csv(SAM_spatial, "./Maturity research/Output/SAM_spatial.csv")


ggplot(SAM_spatial2, aes(LATITUDE, LONGITUDE, fill = SAM))+
  geom_point(size =1.5, shape = 21, stroke = NA)+
  facet_wrap(~YEAR)+
  theme_bw()+
  scale_fill_gradient2(midpoint = mean(SAM_spatial2$SAM))



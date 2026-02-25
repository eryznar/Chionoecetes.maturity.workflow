model <- readRDS("./Maturity research/Models/snowfemale_sdmTMB_spVAR_k300.rda")

# filter specimen data by year and transform to sdmTMB coordinates
readRDS("./Maturity research/Data/snow_survey_specimenEBS.rda")$specimen %>%
  filter( SEX == 2, SHELL_CONDITION == 2) %>%
  mutate(BIN_5MM = cut_width(SIZE_1MM, width = 5, center = 2.5, closed = "left", dig.lab = 4),
         BIN2 = BIN_5MM) %>%
  separate(BIN2, sep = ",", into = c("LOWER", "UPPER")) %>%
  mutate(LOWER = as.numeric(sub('.', '', LOWER)),
         UPPER = as.numeric(gsub('.$', '', UPPER)),
         SIZE_5MM = (UPPER + LOWER)/2,
         YEAR_F = as.factor(YEAR),
         YEAR_SCALED = scale(YEAR)) %>%
  dplyr::select(!c(BIN_5MM, LOWER, UPPER)) %>%
  st_as_sf(., coords = c("LONGITUDE", "LATITUDE"), crs = "+proj=longlat +datum=WGS84") %>%
  st_transform(., crs = "+proj=utm +zone=2") %>%
  cbind(st_coordinates(.)) %>%
  as.data.frame(.) %>%
  mutate(LATITUDE = Y/1000, # scale to km so values don't get too large
         LONGITUDE = X/1000)  -> sub1


pmat <- sub1 %>%
          mutate(MATURE = case_when(CLUTCH_SIZE >0 ~ 1,
                                    TRUE ~ 0)) %>%
        group_by(YEAR, SIZE_5MM) %>%
        reframe(TOT_CRAB = sum(SAMPLING_FACTOR),
               TOT_MAT = sum(SAMPLING_FACTOR[MATURE == 1]),
               TOT_IMM = sum(SAMPLING_FACTOR[MATURE == 0]),
               PMAT = TOT_MAT/TOT_CRAB) %>%
  mutate(PMAT = case_when(SIZE_5MM >=72 ~ 1,
                                 SIZE_5MM <=33 ~ 0,
                                 TRUE ~ PMAT))

ggplot(pmat, aes(SIZE_5MM, PMAT))+
  geom_line(color = "blue")+
  facet_wrap(~YEAR)+
  theme_bw()

pmat.sim <- predict(model, sub1, type = "response")


# Ogives ----
ogives <- pmat.sim %>%
  group_by(YEAR, SPECIES, DISTRICT, SIZE_5MM) %>%
  reframe(
    denom = sum(SAMPLING_FACTOR, na.rm = TRUE),
    num   = sum(est * SAMPLING_FACTOR, na.rm = TRUE),
    PROP_MATURE = ifelse(denom > 0, num / denom, 0)) %>%
   mutate(PROP_MATURE = case_when(SIZE_5MM >=72 ~ 1,
                       SIZE_5MM <=33 ~ 0,
                       TRUE ~ PROP_MATURE))


# Plot
ggplot(ogives, aes(SIZE_5MM, PROP_MATURE))+
  geom_line()+
  facet_wrap(~YEAR)+
  theme_bw()+
  geom_rug()

#write.csv(ogives, "./Maturity research/SNOW_femaleSAM.csv")


# SAM ----
# function
get_sam <- function(size, p) {
  o <- order(size)
  size <- size[o]; p <- p[o]
  
  # if curve never crosses 0.5, return NA
  if (all(p < 0.5) || all(p > 0.5)) return(NA_real_)
  
  # first index where p >= 0.5
  i_upper <- which(p >= 0.5)[1]
  i_lower <- i_upper - 1
  
  # guard for edge cases
  if (is.na(i_upper) || i_upper <= 1) return(NA_real_)
  
  size_lower <- size[i_lower]; size_upper <- size[i_upper]
  prop_lower <- p[i_lower];    prop_upper <- p[i_upper]
  
  # linear interpolation between the two size bins bracketing p = 0.5
  size_lower + ((0.5 - prop_lower) / (prop_upper - prop_lower)) *
    (size_upper - size_lower)
}


# run function
SAM <- ogives %>%
  group_by(YEAR, SPECIES, DISTRICT) %>%
  summarise(
    SAM = get_sam(SIZE_5MM, PROP_MATURE),
    .groups = "drop"
  )

# plot
ggplot(SAM, aes(YEAR, SAM))+
  geom_line()+
  geom_point()

write.csv(SAM, "./Maturity research/SNOW_femaleSAM.csv")

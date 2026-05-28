# LOAD LIBS/PARAMS ---------------------------------------------------------------------------------------
source("./Maturity data processing/Scripts/load_libs_params.R")

# SNOW DATA ----
# Load minima data, calculate cutline params
snow_minima <- read.csv("./Maturity data processing/Output/chela_cutline_minima.csv") %>%
  filter(SPECIES == "SNOW") %>%
  mutate(BETA0 = coef(lm(MINIMUM ~ MIDPOINT))[1],
         BETA1 = coef(lm(MINIMUM ~ MIDPOINT))[2])

BETA0 <- unique(snow_minima$BETA0)
BETA1 <- unique(snow_minima$BETA1)

# Specimen data
crab_dat <- readRDS("./Maturity data processing/Data/snow_survey_specimenEBS.rda")

# Add in 5mm bins to crab_data
spec <-  crab_dat$specimen %>%
  filter(YEAR != 2012, SHELL_CONDITION == 2, SEX == 1) %>%
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
         YEAR_SCALED = as.numeric(scale(YEAR))) %>%
  dplyr::select(!c(BIN_10MM, LOWER, UPPER))


snow_chela <-  read.csv("./Maturity data processing/Data/snow_tanner_cheladatabase.csv") %>% #already filtered appropriately
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
  st_as_sf(., coords = c("LONGITUDE", "LATITUDE"), crs = "+proj=longlat +datum=WGS84") %>%
  st_transform(., crs = "+proj=utm +zone=2") %>%
  cbind(st_coordinates(.)) %>%
  as.data.frame(.) %>%
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
  dplyr::select(YEAR, YEAR_F, YEAR_SCALED, STATION_ID, LATITUDE, LONGITUDE, SIZE_5MM, SIZE_CATEGORY, MATURE) %>%
  filter(YEAR != 2012)

# TANNER DATA ----
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
         YEAR_F = as.factor(YEAR)) %>%
  rename(SIZE_5MM = SIZE_BINNED) %>%
  dplyr::select(YEAR, YEAR_F, STATION_ID, LATITUDE, LONGITUDE, SIZE_5MM, SIZE_CATEGORY, 
                SAMPLING_FACTOR) 

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
  dplyr::select(YEAR, YEAR_F, YEAR_SCALED, STATION_ID, LATITUDE, LONGITUDE, SIZE_5MM, SIZE_CATEGORY, MATURE) 


# DIAGOSTIC FUNCTION ----
plot.resids <- function(model, species, model_name){
  resids <- simulate(model, nsim = 300, type= "mle-mvn")|>
    dharma_residuals(model, return_DHARMa = TRUE)
  
  dat <- cbind(model$data, DHARMa_resid = resids$scaledResiduals)
  
  if(species == "TANNER"){
    dat <- dat %>% filter(YEAR != 2011)
  } else{
    dat <- dat
  }
  
  rr_yr  <- dat %>%
    group_by(YEAR) %>%
    arrange(DHARMa_resid, .by_group = TRUE) %>%
    mutate(
      n = n(),
      expected = ppoints(n),         # uniform quantiles
      observed = sort(DHARMa_resid)  # sort residuals for QQ
    ) %>%
    ungroup() %>%
    mutate(model = model_name)
  
  #  QQ plot with ggplot2
  ggplot()+
    theme_bw()+
    geom_point(rr_yr, mapping = aes(expected, observed), size = 1, fill = "black")+ #theoretical uniform quantiles vs. empirical residual quantiles
    geom_abline(slope = 1, intercept = 0, color = "red", linewidth = 1)+
    ylab("observed")+
    xlab("expected")+
    facet_wrap(~YEAR)+
    scale_x_continuous(breaks = c(0, 0.5, 1))+
    scale_y_continuous(breaks = c(0, 0.5, 1))+
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          strip.text = element_text(size = 12)) +
    ggtitle(paste0(species, " ", model_name)) -> by_yr
 
  rr_size <- dat %>%
    group_by(SIZE_5MM) %>%
    arrange(DHARMa_resid, .by_group = TRUE) %>%
    mutate(
      n = n(),
      expected = ppoints(n),         # uniform quantiles
      observed = sort(DHARMa_resid)  # sort residuals for QQ
    ) %>%
    ungroup() %>%
    mutate(model = model_name)
  
  #  QQ plot with ggplot2
  ggplot()+
    theme_bw()+
    geom_point(rr_size, mapping = aes(expected, observed), size = 1, fill = "black")+ #theoretical uniform quantiles vs. empirical residual quantiles
    geom_abline(slope = 1, intercept = 0, color = "red", linewidth = 1)+
    ylab("observed")+
    xlab("expected")+
    facet_wrap(~SIZE_5MM)+
    scale_x_continuous(breaks = c(0, 0.5, 1))+
    scale_y_continuous(breaks = c(0, 0.5, 1))+
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          strip.text = element_text(size = 12)) +
    ggtitle(ggtitle(paste0(species, " ", model_name))) -> by_size
  
  dat2 <- dat %>%
     group_by(STATION_ID) %>%
     mutate(LONGITUDE = mean(LONGITUDE), LATITUDE = mean(LATITUDE)) %>%
     ungroup()
  
  ggplot(dat2, aes(LONGITUDE, LATITUDE, fill = DHARMa_resid))+
    geom_point(shape = 21, size = 1.75, stroke = NA)+
    facet_wrap(~YEAR)+
    scale_fill_gradient2(midpoint = 0.5)+
    theme_bw() +
    theme(legend.position = "bottom",
          legend.direction = "horizontal",
          strip.text = element_text(size = 10)) +
    ggtitle(ggtitle(paste0(species, " ", model_name))) -> by_yr_sp
  
  ggsave(paste0("./Maturity data processing/Figures/", species, "_", model_name, "spatialDHARMa_byYEAR.png"), width = 10, height = 9)
  
  
  ggplot(dat2, aes(LONGITUDE, LATITUDE, fill = DHARMa_resid))+
    geom_point(shape = 21, size = 1.75, stroke = NA)+
    facet_wrap(~SIZE_5MM)+
    scale_fill_gradient2(midpoint = 0.5)+
    theme_bw() +
    theme(legend.position = "bottom",
      legend.direction = "horizontal",
      strip.text = element_text(size = 10)) +
    ggtitle(ggtitle(paste0(species, " ", model_name)))-> by_size_sp
  
  ggsave(paste0("./Maturity data processing/Figures/", species, "_", model_name, "spatialDHARMa_bySIZE.png"), width = 10, height = 9)
  
  return(list(by_yr = by_yr, by_size = by_size, by_yr_sp = by_yr_sp, by_size_sp = by_size_sp,
         rr_yr = rr_yr, rr_size = rr_size))
}

# SNOW RESIDUALS ----
snow_spvar <- readRDS("./Maturity data processing/Doc/Snow models/sdmTMB_spVAR_SIZE_k300.rda")

plot.resids(snow_spvar, "Snow", "sdmTMB_spatialvar_SIZE_k300") -> out.spvar.s

ggsave(plot = out.spvar.s$by_yr, "./Maturity data processing/Doc/snow_QQ_year.png", width = 10, height =10)
ggsave(plot = out.spvar.s$by_size, "./Maturity data processing/Doc/snow_QQ_size.png", width = 10, height =10)
ggsave(plot = out.spvar.s$by_yr_sp, "./Maturity data processing/Doc/snow_spatialDHARMa_year.png", width = 10, height =10)
ggsave(plot = out.spvar.s$by_size_sp, "./Maturity data processing/Doc/snow_spatialDHARMa_size.png", width = 10, height =10)




snow_spvar300 <- readRDS(paste0(remote_dir, "SNOW/sdmTMB/sdmTMB_spVAR_SIZE_k300.rda"))

plot.resids(snow_spvar300, "Snow", "sdmTMB_spatialvar_SIZE_k300") -> out.spvar300.s

snow_novar <- readRDS("./Maturity data processing/Doc/Snow models/s(SIZE, k=13)_iid_200_sdmTMB.rda")

plot.resids(snow_novar, "Snow", "sdmTMB_novar_k200") -> out.novar.s

rS_yr <- rbind(out.spvar.s$rr_yr %>% dplyr::select(!fold), out.spvar300.s$rr_yr %>% dplyr::select(!fold), out.novar.s$rr_yr)
rS_size <- rbind(out.spvar.s$rr_size %>% dplyr::select(!fold), out.spvar300.s$rr_size %>% dplyr::select(!fold), out.novar.s$rr_size)

#  QQ plot with ggplot2
ggplot()+
  theme_bw()+
  geom_line(rS_yr, mapping = aes(expected, observed,color = model), linewidth = 1, alpha = 0.75)+ #theoretical uniform quantiles vs. empirical residual quantiles
  geom_abline(slope = 1, intercept = 0, color = "black", linewidth = 1, linetype = "dashed")+
  ylab("observed")+
  xlab("expected")+
  facet_wrap(~YEAR)+
  scale_x_continuous(breaks = c(0, 0.5, 1))+
  scale_y_continuous(breaks = c(0, 0.5, 1))+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12))+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.direction = "horizontal")+
  ggtitle("Snow DHARMa by year")


ggsave("./Maturity data processing/Figures/snow_DHARMAbyYEAR.png", width = 10, height = 10)

ggplot()+
  theme_bw()+
  geom_line(rS_size, mapping = aes(expected, observed,color = model), linewidth = 1, alpha = 0.75)+ #theoretical uniform quantiles vs. empirical residual quantiles
  geom_abline(slope = 1, intercept = 0, color = "black", linewidth = 1, linetype = "dashed")+
  ylab("observed")+
  xlab("expected")+
  facet_wrap(~SIZE_5MM)+
  scale_x_continuous(breaks = c(0, 0.5, 1))+
  scale_y_continuous(breaks = c(0, 0.5, 1))+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.direction = "horizontal")+
  ggtitle("Snow DHARMa by size")

ggsave("./Maturity data processing/Figures/snow_DHARMAbySIZE.png", width = 10, height = 10)



# TANNER RESIDUALS ----
tanner_spvar <- readRDS("./Maturity data processing/Doc/Tanner models/sdmTMB_spVAR_SIZE_k200.rda")

plot.resids(tanner_spvar, "Tanner", "sdmTMB_spatialvar_SIZE_k200") -> out.spvar.t


ggsave(plot = out.spvar.t$by_yr, "./Maturity data processing/Doc/tanner_QQ_year.png", width = 10, height =10)
ggsave(plot = out.spvar.t$by_size, "./Maturity data processing/Doc/tanner_QQ_size.png", width = 10, height =10)
ggsave(plot = out.spvar.t$by_yr_sp, "./Maturity data processing/Doc/tanner_spatialDHARMa_year.png", width = 10, height =10)
ggsave(plot = out.spvar.t$by_size_sp, "./Maturity data processing/Doc/tanner_spatialDHARMa_size.png", width = 10, height =10)






tanner_spvar300 <- readRDS(paste0(remote_dir, "TANNER/sdmTMB/sdmTMB_spVAR_SIZE_k300.rda"))

plot.resids(tanner_spvar300, "Tanner", "sdmTMB_spatialvar_SIZE_k300") -> out.spvar300.t

tanner_novar <- readRDS(paste0(remote_dir, "TANNER/sdmTMB/sdmTMB_nospVAR.rda"))

plot.resids(tanner_novar, "Tanner", "sdmTMB_novar_k200") -> out.novar.t



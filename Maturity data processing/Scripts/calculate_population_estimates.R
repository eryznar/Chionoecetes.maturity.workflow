# PURPOSE: to compare biomass and abundance estimates across multiple estimation methods (GAMs, sdmTMB, tinyVAST)

# Author: Emily Ryznar

# NOTES:


# LOAD LIBS/PARAMS ---------------------------------------------------------------------------------------
source("./Maturity data processing/Scripts/load_libs_params.R")

# TANNER ---------------------------------------------------------------------------------------------------
 # Load sdmTMB models ----
 tanner_mod <- readRDS("./Maturity data processing/Doc/Tanner models/sdmTMB_spVAR_SIZE_k200.rda")
   
 # Load specimen data ----
 tanner_dat <- readRDS("./Maturity data processing/Data/tanner_survey_specimenEBS.rda")
 
 # Add in 5mm bins
 tanner_dat$specimen <-  tanner_dat$specimen %>%
   filter(SHELL_CONDITION == 2, SEX == 1) %>%
   mutate(BIN_5MM = cut_width(SIZE_1MM, width = 5, center = 2.5, closed = "left", dig.lab = 4),
          BIN2 = BIN_5MM) %>%
   separate(BIN2, sep = ",", into = c("LOWER", "UPPER")) %>%
   mutate(LOWER = as.numeric(sub('.', '', LOWER)),
          UPPER = as.numeric(gsub('.$', '', UPPER)),
          SIZE_5MM = (UPPER + LOWER)/2) %>%
   dplyr::select(!c(BIN_5MM, LOWER, UPPER)) %>%
   mutate(YEAR_F = as.factor(YEAR),
          YEAR_SCALED = scale(YEAR)) %>%
   filter(YEAR >1989 & !YEAR %in% c(2013, 2015))
 
 # Load minima data, calculate cutline params
 tanner_minima <- read.csv("./Maturity data processing/Output/chela_cutline_minima.csv") %>%
   filter(SPECIES == "TANNER") %>%
   mutate(BETA0 = coef(lm(MINIMUM ~ MIDPOINT))[1],
          BETA1 = coef(lm(MINIMUM ~ MIDPOINT))[2])
 
 BETA0 <- unique(tanner_minima$BETA0)
 BETA1 <- unique(tanner_minima$BETA1)
 
 # Load tanner chela data from chela database (subsample 2)
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
 
 # Calculate biomass and abundance ----
 # Legacy
 tanner_legacy_spec <- tanner_dat

 tanner_pars <- read.csv("./Maturity data processing/Data/tanner_legacy_modelpars.csv") 
 
 tanner_legacy_spec$specimen <-  tanner_legacy_spec$specimen %>% 
   right_join(., tanner_pars) %>%
   mutate(PROP_MATURE = (1/(1 + exp(-A_EST * (SIZE_5MM - B_EST))))) %>%
   mutate(SAMPLING_FACTOR = SAMPLING_FACTOR * PROP_MATURE) %>%
   dplyr::select(!X)
 
 # Use crabpack to calculate bioabund
 tanner_legacy_bioabund <-  crabpack::calc_bioabund(crab_data = tanner_legacy_spec, species = "TANNER", 
                                                  size_min = NULL, size_max = NULL,  sex = "male", 
                                                  shell_condition = "new_hardshell") %>%
   right_join(., data.frame(YEAR = seq(min(.$YEAR), max(.$YEAR)))) %>%
   mutate(Estimator = "Legacy",
          ABUNDANCE = ABUNDANCE/1e6,
          ABUNDANCE_CI = ABUNDANCE_CI/1e6)
 
 
 # sdmTMB 
 tanner_sdmTMB_spec <- tanner_dat
 
 # filter specimen data by year and transform to sdmTMB coordinates
 tanner_sdmTMB_spec$specimen %>%
   st_as_sf(., coords = c("LONGITUDE", "LATITUDE"), crs = "+proj=longlat +datum=WGS84") %>%
   st_transform(., crs = "+proj=utm +zone=2") %>%
   cbind(st_coordinates(.)) %>%
   as.data.frame(.) %>%
   mutate(LATITUDE = Y/1000, # scale to km so values don't get too large
          LONGITUDE = X/1000,
          YEAR_F = as.factor(YEAR)) %>%
   filter(YEAR %in% unique(tanner_mod$data$YEAR)) -> sub1
 

 # Simulate model predictive uncertainty 
 tanner.pmat.sim <- predict(tanner_mod, sub1, type = "response", nsim = 500)
 nsim <- ncol(tanner.pmat.sim)
 bioabund.df <- data.frame()
 
 for(ii in 1:nsim){
   print(paste0("Calculating bioabund sim ", ii))
   
   fit.sim <- tanner.pmat.sim[,ii]
   
   # replace PROP_MATURE with each model simulation draw
   tanner_sdmTMB_spec$specimen <- cbind(sub1, fit.sim) %>%
     rename(PROP_MATURE = fit.sim) %>%
     filter(SPECIES == "TANNER") %>%
     mutate(SAMPLING_FACTOR = SAMPLING_FACTOR * PROP_MATURE)
   
   # calculate bioabund for each simulation
   tanner_sdmTMB_bioabund_sim <-  crabpack::calc_bioabund(crab_data = tanner_sdmTMB_spec, species = "TANNER", 
                                                        size_min = NULL, size_max = NULL,  sex = "male", 
                                                        shell_condition = "new hardshell")
   
   # Bind
   bioabund.df <- rbind(bioabund.df,  tanner_sdmTMB_bioabund_sim %>% mutate(sim = ii))
 }
 
# Now propagate both model and survey uncertainty into biomass/abundance via MCMC
 n_mc <- 1000 # number of resamples per simulation
 set.seed(1)
 
tanner_results_df <-
   bioabund.df %>%
   group_by(YEAR, DISTRICT) %>%
   group_modify(~{
     df <- .x
     n_sim <- nrow(df)
     
     abund_samples    <- numeric(n_sim * n_mc)
     biom_mt_samples  <- numeric(n_sim * n_mc)
     biom_lbs_samples <- numeric(n_sim * n_mc)
     
     for (ii in seq_len(n_sim)) {
       
       abundance_mean <- df$ABUNDANCE[ii]
       abundance_sd   <- df$ABUNDANCE_CI[ii] / 1.96
       
       biomass_mt_mean <- df$BIOMASS_MT[ii]
       biomass_mt_sd   <- df$BIOMASS_MT_CI[ii] / 1.96
       
       biomass_lbs_mean <- df$BIOMASS_LBS[ii]
       biomass_lbs_sd   <- df$BIOMASS_LBS_CI[ii] / 1.96
       
       idx <- ((ii - 1) * n_mc + 1):(ii * n_mc)
       
       abund_samples[idx]    <- rnorm(n_mc, abundance_mean,    abundance_sd)
       biom_mt_samples[idx]  <- rnorm(n_mc, biomass_mt_mean,  biomass_mt_sd)
       biom_lbs_samples[idx] <- rnorm(n_mc, biomass_lbs_mean, biomass_lbs_sd)
     }
     
     tibble(
       ABUNDANCE_MEAN = mean(abund_samples),
       ABUNDANCE_SD   = sd(abund_samples),
       ABUNDANCE_CI   = sd(abund_samples) * 1.96,
       BIOMASS_MT_MEAN = mean(biom_mt_samples),
       BIOMASS_MT_SD   = sd(biom_mt_samples),
       BIOMASS_MT_CI   = sd(biom_mt_samples) * 1.96,
       BIOMASS_LBS_MEAN = mean(biom_lbs_samples),
       BIOMASS_LBS_SD   = sd(biom_lbs_samples),
       BIOMASS_LBS_CI   = sd(biom_lbs_samples) * 1.96
     )
   }) %>%
 ungroup() # this marginally increases the CI after incorporating model uncertainty

 tanner_sdmTMB_bioabund <- tanner_results_df %>%
   dplyr::select(!c(ABUNDANCE_SD, BIOMASS_MT_SD, BIOMASS_LBS_SD)) %>%
   rename(ABUNDANCE = ABUNDANCE_MEAN,
          BIOMASS_MT = BIOMASS_MT_MEAN,
          BIOMASS_LBS = BIOMASS_LBS_MEAN) %>%
   right_join(., data.frame(YEAR = seq(min(.$YEAR), max(.$YEAR)))) %>%
   mutate(Estimator = "sdmTMB",
          ABUNDANCE = ABUNDANCE/1e6,
          ABUNDANCE_CI = ABUNDANCE_CI/1e6) 
 
 
 # Bind and plot 
 tanner_bioabund_dat <- rbind(tanner_legacy_bioabund%>%
                              dplyr::select(names(.)[names(.) %in% colnames(tanner_sdmTMB_bioabund)]), tanner_sdmTMB_bioabund) %>%
   filter(YEAR >= 1990) %>%
   mutate(
     ABUNDANCE = case_when(YEAR %in% c(2011, 2013, 2015, 2020) & DISTRICT == "E166" ~ NA,
                           YEAR %in% c(2013, 2015, 2020) & DISTRICT == "W166" ~ NA,
                           TRUE ~ ABUNDANCE),
     ABUNDANCE_CI = case_when(YEAR %in% c(2011, 2013, 2015, 2020) & DISTRICT == "E166" ~ NA,
                              YEAR %in% c(2013, 2015, 2020) & DISTRICT == "W166" ~ NA,
                              TRUE ~ ABUNDANCE_CI),
     BIOMASS_MT = case_when(YEAR %in% c(2011, 2013, 2015, 2020) & DISTRICT == "E166" ~ NA,
                            YEAR %in% c(2013, 2015, 2020) & DISTRICT == "W166" ~ NA,
                            TRUE ~ BIOMASS_MT),
     BIOMASS_MT_CI = case_when(YEAR %in% c(2011, 2013, 2015, 2020) & DISTRICT == "E166" ~ NA,
                               YEAR %in% c(2013, 2015, 2020) & DISTRICT == "W166" ~ NA,
                               TRUE ~ BIOMASS_MT_CI),
     BIOMASS_LBS = case_when(YEAR %in% c(2011, 2013, 2015, 2020) & DISTRICT == "E166" ~ NA,
                             YEAR %in% c(2013, 2015, 2020) & DISTRICT == "W166" ~ NA,
                             TRUE ~ BIOMASS_LBS),
     BIOMASS_LBS_CI = case_when(YEAR %in% c(2011, 2013, 2015, 2020) & DISTRICT == "E166" ~ NA,
                                YEAR %in% c(2013, 2015, 2020) & DISTRICT == "W166" ~ NA,
                                TRUE ~ BIOMASS_LBS_CI))
 
 # Plot
 e1 <- ggplot(tanner_bioabund_dat %>% filter(DISTRICT == "E166"), aes(YEAR, ABUNDANCE))+
   geom_point(tanner_bioabund_dat %>% filter(DISTRICT == "E166", YEAR %in% c(2012, 2014)),
              mapping = aes(color = Estimator))+
   geom_errorbar(tanner_bioabund_dat %>% filter(DISTRICT == "E166", YEAR %in% c(2012, 2014)),
                 mapping = aes(x = YEAR, ymin = ABUNDANCE - ABUNDANCE_CI, ymax = ABUNDANCE + ABUNDANCE_CI, 
                               color= Estimator), linewidth = 1)+
   geom_ribbon(mapping = aes(x = YEAR, ymin = ABUNDANCE - ABUNDANCE_CI, ymax = ABUNDANCE + ABUNDANCE_CI, 
                             fill= Estimator), alpha = 0.15)+
   geom_line(mapping = aes(color = Estimator), linewidth = 0.75)+
   scale_color_manual(values = c("darkgoldenrod", "cadetblue"), name = "")+
   scale_fill_manual(values = c("darkgoldenrod", "cadetblue"), name = "")+
   theme_bw()+
   ylab("Abundance (millions)")+
   ggtitle("Tanner East morphometric mature males (newshell)")+
   xlab("Year")+
   theme(axis.text.x = element_blank(),
         axis.title.x = element_blank(),
         axis.text.y = element_text(size = 12), 
         axis.title.y = element_text(size = 12),
         legend.position = "bottom",
         legend.direction = "horizontal")
 
 e2 <- ggplot(tanner_bioabund_dat %>% filter(DISTRICT == "E166"), aes(YEAR, BIOMASS_MT))+
   geom_point(tanner_bioabund_dat %>% filter(DISTRICT == "E166", YEAR %in% c(2012, 2014)),
              mapping = aes(color = Estimator))+
   geom_errorbar(tanner_bioabund_dat %>% filter(DISTRICT == "E166", YEAR %in% c(2012, 2014)),
                 mapping = aes(x = YEAR, ymin = BIOMASS_MT - BIOMASS_MT_CI, ymax = BIOMASS_MT + BIOMASS_MT_CI, 
                               color= Estimator), linewidth = 1)+
   geom_ribbon(mapping = aes(x = YEAR, ymin = BIOMASS_MT - BIOMASS_MT_CI, ymax = BIOMASS_MT + BIOMASS_MT_CI, 
                             fill= Estimator), alpha = 0.15)+
   geom_line(mapping = aes(color = Estimator), linewidth = 0.75)+
   scale_color_manual(values = c("darkgoldenrod", "cadetblue"), name = "")+
   scale_fill_manual(values = c("darkgoldenrod", "cadetblue"), name = "")+
   theme_bw()+
   ylab("Biomass (tons)")+
   #ggtitle("Morphometric mature male biomass (SH2)")+
   xlab("Year")+
   theme(axis.text = element_text(size = 12),
         legend.position = "bottom",
         legend.direction = "horizontal",
         axis.title = element_text(size = 12))
 
 e1/e2 + plot_layout(guides = "collect") & 
   theme(legend.position = "bottom", legend.text = element_text(size = 12))
 
 ggsave("./Maturity data processing/Doc/tannerE_bioabund_comparison_legacy.sdmTMB.png", width = 7, height = 6, units = "in")
 
 w1 <- ggplot(tanner_bioabund_dat %>% filter(DISTRICT == "W166"), aes(YEAR, ABUNDANCE))+
   geom_point(tanner_bioabund_dat %>% filter(DISTRICT == "W166", YEAR %in% c(2014)),
              mapping = aes(color = Estimator))+
   geom_errorbar(tanner_bioabund_dat %>% filter(DISTRICT == "W166", YEAR %in% c(2014)),
                 mapping = aes(x = YEAR, ymin = ABUNDANCE - ABUNDANCE_CI, ymax = ABUNDANCE + ABUNDANCE_CI, 
                             color= Estimator), linewidth = 1)+
   geom_ribbon(mapping = aes(x = YEAR, ymin = ABUNDANCE - ABUNDANCE_CI, ymax = ABUNDANCE + ABUNDANCE_CI, 
                             fill= Estimator), alpha = 0.15)+
   geom_line(mapping = aes(color = Estimator), linewidth = 0.75)+
   scale_color_manual(values = c("darkgoldenrod", "cadetblue"), name = "")+
   scale_fill_manual(values = c("darkgoldenrod", "cadetblue"), name = "")+
   theme_bw()+
   ylab("Abundance (millions)")+
   ggtitle("Tanner West morphometric mature males (newshell)")+
   xlab("Year")+
   theme(axis.text.x = element_blank(),
         axis.title.x = element_blank(),
         axis.text.y = element_text(size = 12), 
         axis.title.y = element_text(size = 12),
         legend.position = "bottom",
         legend.direction = "horizontal")
 
 w2 <- ggplot(tanner_bioabund_dat %>% filter(DISTRICT == "W166"), aes(YEAR, BIOMASS_MT))+
   geom_point(tanner_bioabund_dat %>% filter(DISTRICT == "W166", YEAR %in% c(2014)),
              mapping = aes(color = Estimator))+
   geom_errorbar(tanner_bioabund_dat %>% filter(DISTRICT == "W166", YEAR %in% c(2014)),
                 mapping = aes(x = YEAR, ymin = BIOMASS_MT - BIOMASS_MT_CI, ymax = BIOMASS_MT + BIOMASS_MT_CI, 
                               color= Estimator), linewidth = 1)+
   geom_ribbon(mapping = aes(x = YEAR, ymin = BIOMASS_MT - BIOMASS_MT_CI, ymax = BIOMASS_MT + BIOMASS_MT_CI, 
                             fill= Estimator), alpha = 0.15)+
   geom_line(mapping = aes(color = Estimator), linewidth = 0.75)+
   scale_color_manual(values = c("darkgoldenrod", "cadetblue"), name = "")+
   scale_fill_manual(values = c("darkgoldenrod", "cadetblue"), name = "")+
   theme_bw()+
   ylab("Biomass (tons)")+
   #ggtitle("Morphometric mature male biomass (SH2)")+
   xlab("Year")+
   theme(axis.text = element_text(size = 12),
         legend.position = "bottom",
         legend.direction = "horizontal",
         axis.title = element_text(size = 12))
 
 w1/w2 + plot_layout(guides = "collect") & 
   theme(legend.position = "bottom", legend.text = element_text(size = 12))
 
 ggsave("./Maturity data processing/Doc/tannerW_bioabund_comparison_legacy.sdmTMB.png", width = 7, height = 6, units = "in")
 
 # Calculate ogives ----
  # legacy using gam prop_mature interpolation of 10mm down to 5mm bins
 tanner_legacy_ogives <- rbind((read.csv("./Maturity data processing/Doc/tanE_pmat_5mminterp.csv") %>% 
                                  rename(SIZE = SIZE_BIN) %>% mutate(DISTRICT = "E166")),
                               (read.csv("./Maturity data processing/Doc/tanW_pmat_5mminterp.csv") %>%
                                 rename(SIZE = SIZE_BIN) %>% mutate(DISTRICT = "W166")))%>%
   dplyr::select(!X)  %>%
   mutate(Estimator = "Legacy")
 
 tanner.ogive.spec <-  readRDS("./Maturity data processing/Data/tanner_survey_specimenEBS.rda")$specimen %>%
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
   filter(YEAR >1989 & !YEAR %in% c(2013, 2015))
 
 
 tanner_sdmTMB_ogives <- cbind(tanner.ogive.spec, tanner.pmat.sim) %>% # calculate mean prop_mature across simulations
   group_by(YEAR, SIZE_5MM, DISTRICT,) %>%
   summarise(
     across(
       matches("^[0-9]+$"),   # selects "1","2",...,"500"
       ~ sum(.x * SAMPLING_FACTOR) / sum(SAMPLING_FACTOR),
       .names = "wmean_{.col}"
     ),
     .groups = "drop"
   ) %>%
   pivot_longer(
     starts_with("wmean_"),
     names_to = "sim",
     values_to = "wmean"
   ) %>%
   group_by(YEAR, SIZE_5MM, DISTRICT) %>%
   summarise(
     PROP_MATURE = mean(wmean),
     PROP_MATURE_SE   = sd(wmean) / sqrt(n()),
     .groups = "drop"
   ) %>%
   rename(SIZE = SIZE_5MM) %>%
   mutate(Estimator = "sdmTMB")
 
 # Bind and plot
 ogive.dat <- rbind(tanner_legacy_ogives, tanner_sdmTMB_ogives %>% dplyr::select(!PROP_MATURE_SE))
 
 
 
 ggplot(ogive.dat %>% filter(DISTRICT == "W166" & !YEAR %in% c(2013, 2015)),  aes(SIZE, PROP_MATURE, color = Estimator))+
   geom_line(linewidth = 1)+
   scale_color_manual(values = c(
     "Legacy"       = "darkgoldenrod",
     "sdmTMB"       = "cadetblue",name = ""))+
   scale_fill_manual(values = c(
     "Legacy"       = "darkgoldenrod",
     "sdmTMB"       = "cadetblue"), name = "")+
   facet_wrap(~YEAR)+
   geom_rug(aes(x = SIZE), sides = "b", alpha = 0.4, inherit.aes = FALSE) +
   geom_hline(yintercept = 0.5, linetype = "dashed")+
   theme_bw()+
   ggtitle("Tanner West")+
   #scale_x_continuous(breaks = seq(0, 175, by = 10))+
   ylab("Proportion mature")+
   xlab("Carapace width (mm)")+
   theme(legend.position = "bottom", legend.direction = "horizontal",
         legend.text = element_text(size = 12),
         legend.title = element_blank(),
         axis.title = element_text(size = 12),
         strip.text = element_text(size = 10),
         axis.text = element_text(size = 10))
 
 ggsave("./Maturity data processing/Doc/tannerW166_ogive_comparison_legacy.sdmTMB.png", width  = 10, height = 11, units = "in")
 

 ggplot(ogive.dat %>% filter(DISTRICT == "E166" & !YEAR %in% c(2011, 2013, 2015)),  aes(SIZE, PROP_MATURE, color = Estimator))+
   geom_line(linewidth = 1)+
   scale_color_manual(values = c(
     "Legacy"       = "darkgoldenrod",
     "sdmTMB"       = "cadetblue",name = ""))+
   scale_fill_manual(values = c(
     "Legacy"       = "darkgoldenrod",
     "sdmTMB"       = "cadetblue"), name = "")+
   facet_wrap(~YEAR)+
   geom_hline(yintercept = 0.5, linetype = "dashed")+
   geom_rug(aes(x = SIZE), sides = "b", alpha = 0.4, inherit.aes = FALSE) +
   theme_bw()+
   ggtitle("Tanner East")+
   #scale_x_continuous(breaks = seq(0, 175, by = 10))+
   ylab("Proportion mature")+
   xlab("Carapace width (mm)")+
   theme(legend.position = "bottom", legend.direction = "horizontal",
         legend.text = element_text(size = 12),
         legend.title = element_blank(),
         axis.title = element_text(size = 12),
         strip.text = element_text(size = 10),
         axis.text = element_text(size = 10))
 
 ggsave("./Maturity data processing/Doc/tannerE166_ogive_comparison_legacy.sdmTMB.png", width  = 10, height = 11, units = "in")
 
 

 # Calculate/plot ogives ----
 ## Get specimen data
 tanner.ogive.spec <-  readRDS("./Maturity data processing/Data/tanner_survey_specimenEBS.rda")$specimen %>%
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
   filter(YEAR >1989 & !YEAR %in% c(2013, 2015))
 
 ## Calculate Weighted mean proportion mature per sim, YEAR, DISTRICT, SIZE_5MM
 tanner_long <- cbind(tanner.ogive.spec, tanner.pmat.sim) %>%
   group_by(YEAR, DISTRICT, SIZE_5MM) %>%
   summarise(
     across(
       matches("^[0-9]+$"),   # "1","2",...,"500"
       ~ sum(.x * SAMPLING_FACTOR) / sum(SAMPLING_FACTOR),
       .names = "wmean_{.col}"
     ),
     .groups = "drop"
   ) %>%
   pivot_longer(
     starts_with("wmean_"),
     names_to  = "sim",
     values_to = "wmean"
   ) %>%
   mutate(
     sim = as.integer(gsub("wmean_", "", sim))
   )
 
 ## 2. Same SAM interpolation function
 
 get_sam <- function(size, p) {
   o <- order(size)
   size <- size[o]; p <- p[o]
   if (all(p < 0.5) || all(p > 0.5)) return(NA_real_)
   i_upper <- which(p >= 0.5)[1]
   i_lower <- i_upper - 1
   if (is.na(i_upper) || i_upper <= 1) return(NA_real_)
   size_lower <- size[i_lower]; size_upper <- size[i_upper]
   prop_lower <- p[i_lower];    prop_upper <- p[i_upper]
   size_lower + ((0.5 - prop_lower) / (prop_upper - prop_lower)) *
     (size_upper - size_lower)
 }
 
 ## 3. SAM for each simulation, YEAR, DISTRICT
 
 tanner_SAM_sims <- tanner_long %>%
   group_by(YEAR, DISTRICT, sim) %>%
   summarise(
     SAM = get_sam(SIZE_5MM, wmean),
     .groups = "drop"
   )
 

 
 
 
 ## 4. Summaries by YEAR, DISTRICT (mean, SD, SE, CI)
 
 tanner_SAM_summary <- tanner_SAM_sims %>%
   group_by(YEAR, DISTRICT) %>%
   summarise(
     SAM_mean = mean(SAM, na.rm = TRUE),
     SAM_sd   = sd(SAM,   na.rm = TRUE),
     SAM_se   = sd(SAM,   na.rm = TRUE) / sqrt(n()),
     SAM_lo   = quantile(SAM, 0.025, na.rm = TRUE),
     SAM_hi   = quantile(SAM, 0.975, na.rm = TRUE),
     .groups  = "drop"
   ) %>%
   filter(!(YEAR == 2011 & DISTRICT == "E166")) %>%
   right_join(., expand.grid(YEAR = seq(min(.$YEAR), max(.$YEAR)), DISTRICT = c("E166", "W166"),
                             Estimator = "sdmTMB")) 
 
 
 
 
 
 
 
 
 
 
 legacy_SAM <- read.csv("./Maturity data processing/Data/tanner_legacy_ogives.csv") %>% 
   filter(DISTRICT != "ALL") %>%
   dplyr::select(YEAR, DISTRICT, B_EST, B_SE) %>% 
   distinct() %>%
   rename(SAM = B_EST) %>%
   right_join(., expand.grid(YEAR = seq(min(.$YEAR), max(.$YEAR)), DISTRICT = c("E166", "W166"),
                             Estimator = "Legacy")) %>%
   distinct() 
 
 
 
 ggplot()+
   #facet_wrap(~SAM_type, ncol = 1)+
   geom_line(legacy_SAM, mapping = aes(YEAR, SAM, color = as.factor(1)), linewidth = 1)+
   geom_point(legacy_SAM, mapping = aes(YEAR, SAM, color = as.factor(1)), size = 2)+
   geom_errorbar(legacy_SAM, mapping = aes(YEAR, ymin = SAM-B_SE*1.96, ymax = SAM + B_SE*1.96, color = as.factor(1)), linewidth = 1)+
   # geom_point(tanner_sdmTMB_SAM, mapping = aes(YEAR, SAM, color = as.factor(2)), size = 2)+
   # geom_line(tanner_sdmTMB_SAM, mapping = aes(YEAR, SAM, color = as.factor(2)), linewidth = 1)+
   geom_line(tanner_SAM_summary, mapping = aes(YEAR, SAM_mean, color = as.factor(2)), linewidth = 1)+
   geom_point(tanner_SAM_summary, mapping = aes(YEAR, SAM_mean, color = as.factor(2)), size = 2)+
   geom_errorbar(tanner_SAM_summary, mapping = aes(YEAR, ymin = SAM_lo, ymax = SAM_hi, color = as.factor(2)), linewidth = 1)+
   scale_color_manual(values = c("darkgoldenrod", "cadetblue"), labels = c("Legacy", "sdmTMB"), name = "")+
   scale_fill_manual(values = c("darkgoldenrod", "cadetblue"), labels = c("Legacy", "sdmTMB"), name = "")+
   facet_wrap(~DISTRICT, nrow = 2, scales = "free_y")+
   theme_bw()+
   xlab("Year")+
   ylab("Size at 50% maturity (mm)")+
   theme(legend.position = "bottom", legend.direction = "horizontal",
         axis.text = element_text(size = 12),
         axis.title = element_text(size = 12),
         legend.text = element_text(size = 12),
         strip.text = element_text(size = 10))
 
 ggsave("./Maturity data processing/Doc/tanner_SAM_comparison.png", width = 8, height = 8, units = "in")
 
 
# SNOW ---------------------------------------------------------------------------------------------------
 # Load sdmTMB models ----
 snow_mod <- readRDS("./Maturity data processing/Doc/Snow models/sdmTMB_spVAR_SIZE_k300.rda")
 
 # Load specimen data ----
 snow_dat <- readRDS("./Maturity data processing/Data/snow_survey_specimenEBS.rda")
 
 # Add in 5mm bins
 snow_dat$specimen <-  snow_dat$specimen %>%
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
          YEAR_F = as.factor(YEAR),
          YEAR_SCALED = scale(YEAR)) %>%
   dplyr::select(!c(BIN_10MM, LOWER, UPPER))
 
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
 
 
 # Calculate biomass and abundance ----
 # Legacy
 snow_legacy_spec <- snow_dat
 snow_pars <- read.csv("./Maturity data processing/Data/snow_legacy_modelpars.csv") 
 
 snow_legacy_spec$specimen <-  snow_legacy_spec$specimen %>% 
   right_join(., snow_pars) %>%
   mutate(PROP_MATURE = (1/(1 + exp(-A_EST * (SIZE_5MM - B_EST))))) %>%
   mutate(SAMPLING_FACTOR = SAMPLING_FACTOR * PROP_MATURE) %>%
   dplyr::select(!X)
 
 # Use crabpack to calculate bioabund
 snow_legacy_bioabund <-  crabpack::calc_bioabund(crab_data = snow_legacy_spec, species = "SNOW", 
                                                    size_min = NULL, size_max = NULL,  sex = "male", 
                                                    shell_condition = "new_hardshell") %>%
   right_join(., data.frame(YEAR = seq(min(.$YEAR), max(.$YEAR)))) %>%
   mutate(Estimator = "Legacy",
          ABUNDANCE = ABUNDANCE/1e6,
          ABUNDANCE_CI = ABUNDANCE_CI/1e6)
 
 
 # sdmTMB 
 snow_sdmTMB_spec <-snow_dat
 
 # filter specimen data by year and transform to sdmTMB coordinates
 snow_sdmTMB_spec$specimen %>%
   st_as_sf(., coords = c("LONGITUDE", "LATITUDE"), crs = "+proj=longlat +datum=WGS84") %>%
   st_transform(., crs = "+proj=utm +zone=2") %>%
   cbind(st_coordinates(.)) %>%
   as.data.frame(.) %>%
   mutate(LATITUDE = Y/1000, # scale to km so values don't get too large
          LONGITUDE = X/1000) %>%
   filter(YEAR %in% unique(snow_mod$data$YEAR)) -> sub1
 
 
 # Simulate model predictive uncertainty 
 snow.pmat.sim <- predict(snow_mod, sub1, type = "response", nsim = 500)
 nsim <- ncol(snow.pmat.sim)
 bioabund.df <- data.frame()
 
 for(ii in 1:nsim){
   print(paste0("Calculating bioabund sim ", ii))
   
   fit.sim <- snow.pmat.sim[,ii]
   
   # replace PROP_MATURE with each model simulation draw
   snow_sdmTMB_spec$specimen <- cbind(sub1, fit.sim) %>%
     rename(PROP_MATURE = fit.sim) %>%
     filter(SPECIES == "SNOW") %>%
     mutate(SAMPLING_FACTOR = SAMPLING_FACTOR * PROP_MATURE)
   
   # calculate bioabund for each simulation
   snow_sdmTMB_bioabund_sim <-  crabpack::calc_bioabund(crab_data = snow_sdmTMB_spec, species = "SNOW", 
                                                          size_min = NULL, size_max = NULL,  sex = "male", 
                                                          shell_condition = "new hardshell")
   
   # Bind
   bioabund.df <- rbind(bioabund.df,  snow_sdmTMB_bioabund_sim %>% mutate(sim = ii))
 }
 
 # Now propagate both model and survey uncertainty into biomass/abundance via MCMC
 n_mc <- 1000 # number of resamples per simulation
 set.seed(1)
 
 snow_results_df <-
   bioabund.df %>%
   group_by(YEAR, DISTRICT) %>%
   group_modify(~{
     df <- .x
     n_sim <- nrow(df)
     
     abund_samples    <- numeric(n_sim * n_mc)
     biom_mt_samples  <- numeric(n_sim * n_mc)
     biom_lbs_samples <- numeric(n_sim * n_mc)
     
     for (ii in seq_len(n_sim)) {
       
       abundance_mean <- df$ABUNDANCE[ii]
       abundance_sd   <- df$ABUNDANCE_CI[ii] / 1.96
       
       biomass_mt_mean <- df$BIOMASS_MT[ii]
       biomass_mt_sd   <- df$BIOMASS_MT_CI[ii] / 1.96
       
       biomass_lbs_mean <- df$BIOMASS_LBS[ii]
       biomass_lbs_sd   <- df$BIOMASS_LBS_CI[ii] / 1.96
       
       idx <- ((ii - 1) * n_mc + 1):(ii * n_mc)
       
       abund_samples[idx]    <- rnorm(n_mc, abundance_mean,    abundance_sd)
       biom_mt_samples[idx]  <- rnorm(n_mc, biomass_mt_mean,  biomass_mt_sd)
       biom_lbs_samples[idx] <- rnorm(n_mc, biomass_lbs_mean, biomass_lbs_sd)
     }
     
     tibble(
       ABUNDANCE_MEAN = mean(abund_samples),
       ABUNDANCE_SD   = sd(abund_samples),
       ABUNDANCE_CI   = sd(abund_samples) * 1.96,
       BIOMASS_MT_MEAN = mean(biom_mt_samples),
       BIOMASS_MT_SD   = sd(biom_mt_samples),
       BIOMASS_MT_CI   = sd(biom_mt_samples) * 1.96,
       BIOMASS_LBS_MEAN = mean(biom_lbs_samples),
       BIOMASS_LBS_SD   = sd(biom_lbs_samples),
       BIOMASS_LBS_CI   = sd(biom_lbs_samples) * 1.96
     )
   }) %>%
   ungroup() # this marginally increases the CI after incorporating model uncertainty
 
 snow_sdmTMB_bioabund <- snow_results_df %>%
   dplyr::select(!c(ABUNDANCE_SD, BIOMASS_MT_SD, BIOMASS_LBS_SD)) %>%
   rename(ABUNDANCE = ABUNDANCE_MEAN,
          BIOMASS_MT = BIOMASS_MT_MEAN,
          BIOMASS_LBS = BIOMASS_LBS_MEAN) %>%
   right_join(., data.frame(YEAR = seq(min(.$YEAR), max(.$YEAR)))) %>%
   mutate(Estimator = "sdmTMB",
          ABUNDANCE = ABUNDANCE/1e6,
          ABUNDANCE_CI = ABUNDANCE_CI/1e6) 
 
 
 # Bind and plot 
 snow_bioabund_dat <- rbind(snow_legacy_bioabund%>%
                                dplyr::select(names(.)[names(.) %in% colnames(snow_sdmTMB_bioabund)]), snow_sdmTMB_bioabund) %>%
   filter(YEAR >= 1989) %>%
   mutate(
     ABUNDANCE = case_when(YEAR %in% c(2008, 2012, 2014, 2016, 2020) ~ NA,
       TRUE ~ ABUNDANCE),
     ABUNDANCE_CI = case_when(YEAR %in% c(2008, 2012, 2014, 2016, 2020) ~ NA,
       TRUE ~ ABUNDANCE_CI),
     BIOMASS_MT = case_when(YEAR %in% c(2008, 2012, 2014, 2016, 2020) ~ NA,
       TRUE ~ BIOMASS_MT),
     BIOMASS_MT_CI = case_when(YEAR %in% c(2008, 2012, 2014, 2016, 2020) ~ NA,
       TRUE ~ BIOMASS_MT_CI),
     BIOMASS_LBS = case_when(YEAR %in% c(2008, 2012, 2014, 2016, 2020) ~ NA,
       TRUE ~ BIOMASS_LBS),
     BIOMASS_LBS_CI = case_when(YEAR %in% c(2008, 2012, 2014, 2016, 2020) ~ NA,
       TRUE ~ BIOMASS_LBS_CI)) %>%
   arrange(Estimator, YEAR)
 
 snow_bioabund_dat %>% filter(YEAR %in% c(2013, 2015))
   
 # Plot
 p1 <- ggplot()+
   geom_point(snow_bioabund_dat %>% filter(YEAR %in% c(2013, 2015)), mapping = aes(YEAR, ABUNDANCE, color = Estimator))+
   geom_errorbar(snow_bioabund_dat %>% filter(YEAR %in% c(2013, 2015)), mapping = aes(x = YEAR, ymin = ABUNDANCE - ABUNDANCE_CI, ymax = ABUNDANCE + ABUNDANCE_CI, 
                                                                                      color= Estimator), linewidth = 1)+
   geom_ribbon(snow_bioabund_dat, mapping = aes(x = YEAR, ymin = ABUNDANCE - ABUNDANCE_CI, ymax = ABUNDANCE + ABUNDANCE_CI, 
                                                fill= Estimator), alpha = 0.15)+
   geom_line(snow_bioabund_dat, mapping = aes(YEAR, ABUNDANCE, color = Estimator), linewidth = 0.75)+
   scale_color_manual(values = c("darkgoldenrod", "cadetblue"), name = "")+
   scale_fill_manual(values = c("darkgoldenrod", "cadetblue"), name = "")+
   theme_bw()+
   ylab("Abundance (millions)")+
   ggtitle("Snow morphometric mature males (newshell)")+
   xlab("Year")+
   theme(axis.text.x = element_blank(),
         axis.title.x = element_blank(),
         axis.text.y = element_text(size = 12), 
         axis.title.y = element_text(size = 12),
         legend.position = "bottom",
         legend.direction = "horizontal")
 
 p2 <- ggplot()+
   geom_point(snow_bioabund_dat %>% filter(YEAR %in% c(2013, 2015)), mapping = aes(YEAR, BIOMASS_MT, color = Estimator))+
   geom_errorbar(snow_bioabund_dat %>% filter(YEAR %in% c(2013, 2015)), mapping = aes(x = YEAR, ymin = BIOMASS_MT - BIOMASS_MT_CI, ymax = BIOMASS_MT + BIOMASS_MT_CI, 
                                                color= Estimator) , linewidth = 1)+
   geom_ribbon(snow_bioabund_dat, mapping = aes(x = YEAR, ymin = BIOMASS_MT - BIOMASS_MT_CI, ymax = BIOMASS_MT + BIOMASS_MT_CI, 
                             fill= Estimator), alpha = 0.15)+
   geom_line(snow_bioabund_dat, mapping = aes(YEAR, BIOMASS_MT, color = Estimator), linewidth = 0.75)+
   scale_color_manual(values = c("darkgoldenrod", "cadetblue"), name = "")+
   scale_fill_manual(values = c("darkgoldenrod", "cadetblue"), name = "")+
   theme_bw()+
   ylab("Biomass (tons)")+
   #ggtitle("Morphometric mature male biomass (SH2)")+
   xlab("Year")+
   theme(axis.text = element_text(size = 12),
         legend.position = "bottom",
         legend.direction = "horizontal",
         axis.title = element_text(size = 12))
 
 p1/p2 + plot_layout(guides = "collect") & 
   theme(legend.position = "bottom", legend.text = element_text(size = 12))
 
 ggsave("./Maturity data processing/Doc/snow_bioabund_comparison_legacy.sdmTMB.png", width = 7, height = 6, units = "in")
 
 
 # Calculate ogives ----
 # legacy using logistic curves
 snow_legacy_ogives <- snow_legacy_spec$specimen %>%
   group_by(YEAR,SIZE_5MM) %>%
   summarise(PROP_MATURE_SE = sd(PROP_MATURE) * 1.96,
             PROP_MATURE = mean(PROP_MATURE), .groups = "drop") %>%
   mutate(Estimator = "Legacy") %>%
   rename(SIZE = SIZE_5MM) %>%
   na.omit() 
 
 # legacy using gam prop_mature interpolation of 10mm down to 5mm bins
 snow_legacy_ogives <- read.csv("./Maturity data processing/Doc/snow_pmat_5mminterp.csv") %>%
   dplyr::select(!X) %>%
   rename(SIZE = SIZE_BIN) %>%
   mutate(Estimator = "Legacy")
 
 # sdmTMB
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
 
   
 snow_sdmTMB_ogives <- cbind(snow.ogive.spec, snow.pmat.sim) %>% # calculate mean prop_mature across simulations
   group_by(YEAR, SIZE_5MM) %>%
   summarise(
     across(
       matches("^[0-9]+$"),   # selects "1","2",...,"500"
       ~ sum(.x * SAMPLING_FACTOR) / sum(SAMPLING_FACTOR),
       .names = "wmean_{.col}"
     ),
     .groups = "drop"
   ) %>%
   pivot_longer(
     starts_with("wmean_"),
     names_to = "sim",
     values_to = "wmean"
   ) %>%
   group_by(YEAR, SIZE_5MM) %>%
   summarise(
     PROP_MATURE = mean(wmean),
     PROP_MATURE_SE   = sd(wmean) / sqrt(n()),
     .groups = "drop"
   ) %>%
   rename(SIZE = SIZE_5MM) %>%
   mutate(Estimator = "sdmTMB")
 

 # Bind and plot
 snow.ogive.dat <- rbind(snow_legacy_ogives, snow_sdmTMB_ogives %>% dplyr::select(!PROP_MATURE_SE))
 
 
 
 ggplot(snow.ogive.dat,  aes(SIZE, PROP_MATURE, color = Estimator))+
   geom_line(linewidth = 1)+
   scale_color_manual(values = c(
     "Legacy"       = "darkgoldenrod",
     "sdmTMB"       = "cadetblue",name = ""))+
   scale_fill_manual(values = c(
     "Legacy"       = "darkgoldenrod",
     "sdmTMB"       = "cadetblue"), name = "")+
   facet_wrap(~YEAR)+
   geom_rug(snow.ogive.dat %>% filter(Estimator == "sdmTMB"), mapping =aes(x = SIZE), sides = "b", inherit.aes = TRUE) +
   geom_hline(yintercept = 0.5, linetype = "dashed")+
   theme_bw()+
   xlim(0, 150)+
   #scale_x_continuous(breaks = seq(0, 175, by = 10))+
   ylab("Proportion mature")+
   xlab("Carapace width (mm)")+
   theme(legend.position = "bottom", legend.direction = "horizontal",
         legend.text = element_text(size = 12),
         legend.title = element_blank(),
         axis.title = element_text(size = 12),
         strip.text = element_text(size = 12),
         axis.text = element_text(size = 10))
 
 ggsave("./Maturity data processing/Doc/snow_ogive_comparison_legacy.sdmTMB.png", width  = 10, height = 11, units = "in")
 

 # SAM ----
 snow_long <- cbind(snow.ogive.spec, snow.pmat.sim) %>%
   group_by(YEAR, SIZE_5MM) %>%
   summarise(
     across(
       matches("^[0-9]+$"),
       ~ sum(.x * SAMPLING_FACTOR) / sum(SAMPLING_FACTOR),
       .names = "wmean_{.col}"
     ),
     .groups = "drop"
   ) %>%
   pivot_longer(
     starts_with("wmean_"),
     names_to  = "sim",
     values_to = "wmean"
   ) %>%
   mutate(
     sim = as.integer(gsub("wmean_", "", sim))
   )
 
 snow_mean <- snow_long %>%
   group_by(YEAR, SIZE_5MM) %>%
   summarise(wmean_mean = mean(wmean), .groups = "drop")
 
 ggplot() +
   # all simulation ogives in faded grey
   geom_line(
     data = snow_long,
     aes(x = SIZE_5MM, y = wmean, group = sim),
     color = "darkgrey",
     alpha = 0.2
   ) +
   # mean ogive in solid black
   geom_line(
     data = snow_mean,
     aes(x = SIZE_5MM, y = wmean_mean),
     color = "cadetblue",
     linewidth = 1
   ) +
   facet_wrap(~ YEAR) +
   labs(x = "Size (mm)", y = "Proportion mature") +
   geom_rug(snow.ogive.dat %>% filter(Estimator == "sdmTMB"), mapping =aes(x = SIZE), sides = "b", inherit.aes = TRUE) +
   theme_bw()+
   xlim(0, 150)+
   #scale_x_continuous(breaks = seq(0, 175, by = 10))+
   ylab("Proportion mature")+
   xlab("Carapace width (mm)")+
   theme(legend.position = "bottom", legend.direction = "horizontal",
         legend.text = element_text(size = 12),
         legend.title = element_blank(),
         axis.title = element_text(size = 12),
         strip.text = element_text(size = 12),
         axis.text = element_text(size = 10))
 
 snow_ogive_sdmTMB <- snow_long %>%
   group_by(YEAR, SIZE_5MM) %>%
   summarise(
     PROP_MATURE_mean = mean(wmean),
     PROP_MATURE_lo   = quantile(wmean, 0.025),
     PROP_MATURE_hi   = quantile(wmean, 0.975),
     .groups = "drop"
   ) %>%
   mutate(Estimator = "sdmTMB")
 
 ogive_plot_dat <- rbind(snow_ogive_sdmTMB %>% rename(PROP_MATURE = PROP_MATURE_mean, SIZE = SIZE_5MM), snow_legacy_ogives %>% mutate(PROP_MATURE_lo = NA, PROP_MATURE_hi = NA))
 
 ggplot(ogive_plot_dat,  aes(SIZE, PROP_MATURE, color = Estimator))+
   geom_ribbon(aes(SIZE, ymin = PROP_MATURE_lo, ymax = PROP_MATURE_hi, fill = Estimator), color = NA, alpha = 0.35)+
   geom_line(linewidth = 1)+
   scale_color_manual(values = c(
     "Legacy"       = "darkgoldenrod",
     "sdmTMB"       = "cadetblue",name = ""))+
   scale_fill_manual(values = c(
     "Legacy"       = "darkgoldenrod",
     "sdmTMB"       = "cadetblue"), name = "", guide= "none")+
   facet_wrap(~YEAR)+
   geom_rug(snow.ogive.dat %>% filter(Estimator == "sdmTMB"), mapping =aes(x = SIZE), sides = "b", inherit.aes = TRUE) +
   geom_hline(yintercept = 0.5, linetype = "dashed")+
   theme_bw()+
   xlim(0, 150)+
   #scale_x_continuous(breaks = seq(0, 175, by = 10))+
   ylab("Proportion mature")+
   xlab("Carapace width (mm)")+
   theme(legend.position = "bottom", legend.direction = "horizontal",
         legend.text = element_text(size = 12),
         legend.title = element_blank(),
         axis.title = element_text(size = 12),
         strip.text = element_text(size = 12),
         axis.text = element_text(size = 10))
 
 
 ggsave("./Maturity data processing/Doc/snow_ogive_comparison_legacy.sdmTMBwitherror.png", width  = 10, height = 11, units = "in")
 
 
 
 
 get_sam <- function(size, p) {
   o <- order(size)
   size <- size[o]; p <- p[o]
   if (all(p < 0.5) || all(p > 0.5)) return(NA_real_)
   i_upper <- which(p >= 0.5)[1]
   i_lower <- i_upper - 1
   if (is.na(i_upper) || i_upper <= 1) return(NA_real_)
   size_lower <- size[i_lower]; size_upper <- size[i_upper]
   prop_lower <- p[i_lower];    prop_upper <- p[i_upper]
   size_lower + ((0.5 - prop_lower) / (prop_upper - prop_lower)) *
     (size_upper - size_lower)
 }
 
 snow_SAM_sims <- snow_long %>%
   group_by(YEAR, sim) %>%
   summarise(
     SAM = get_sam(SIZE_5MM, wmean),
     .groups = "drop"
   )
 
 
 snow_SAM_summary <- snow_SAM_sims %>%
   group_by(YEAR) %>%
   summarise(
     SAM_mean = mean(SAM, na.rm = TRUE),
     SAM_sd   = sd(SAM,   na.rm = TRUE),
     SAM_se  = sd(SAM,   na.rm = TRUE) / sqrt(n()),
     SAM_lo   = quantile(SAM, 0.025, na.rm = TRUE),
     SAM_hi   = quantile(SAM, 0.975, na.rm = TRUE),
     .groups = "drop"
   ) %>%
   right_join(., data.frame(YEAR = seq(min(.$YEAR), max(.$YEAR))))
 
 ggplot(snow_SAM_summary, aes(YEAR, SAM_mean))+
   geom_errorbar(aes(YEAR, ymin = SAM_lo, ymax = SAM_hi))+
   geom_point()+
   geom_line()
 
 
 pp <- lme(SAM ~ YEAR, data = snow_sdmTMB_SAM %>% na.omit(), random = ~ 1 | YEAR, correlation = corAR1())
 
 write.csv(snow_sdmTMB_SAM, "./Maturity data processing/Output/snow_SAM.csv")
 
 legacy_SAM <- read.csv("./Maturity data processing/Data/snow_legacy_ogives.csv") %>% 
   dplyr::select(YEAR, B_EST, B_SE) %>% 
   distinct() %>%
   rename(SAM = B_EST) %>%
   right_join(., data.frame(YEAR = seq(min(.$YEAR), max(.$YEAR)))) 
 
 
 
 ggplot()+
   #facet_wrap(~SAM_type, ncol = 1)+
   geom_line(legacy_SAM, mapping = aes(YEAR, SAM, color = as.factor(1)), linewidth = 1)+
   geom_point(legacy_SAM, mapping = aes(YEAR, SAM, color = as.factor(1)), size = 2)+
   geom_errorbar(legacy_SAM, mapping = aes(YEAR, ymin = SAM-B_SE*1.96, ymax = SAM + B_SE*1.96, color = as.factor(1)), linewidth = 1)+
   # geom_point(snow_sdmTMB_SAM, mapping = aes(YEAR, SAM, color = as.factor(2)), size = 2)+
   # geom_line(snow_sdmTMB_SAM, mapping = aes(YEAR, SAM, color = as.factor(2)), linewidth = 1)+
   geom_line(snow_SAM_summary, mapping = aes(YEAR, SAM_mean, color = as.factor(2)), linewidth = 1)+
   geom_point(snow_SAM_summary, mapping = aes(YEAR, SAM_mean, color = as.factor(2)), size = 2)+
   geom_errorbar(snow_SAM_summary, mapping = aes(YEAR, ymin = SAM_lo, ymax = SAM_hi, color = as.factor(2)), linewidth = 1)+
   scale_color_manual(values = c("darkgoldenrod", "cadetblue"), labels = c("Legacy", "sdmTMB"), name = "")+
   scale_fill_manual(values = c("darkgoldenrod", "cadetblue"), labels = c("Legacy", "sdmTMB"), name = "")+
   #facet_wrap(~DISTRICT)+
   theme_bw()+
   xlab("Year")+
   ylab("Size at 50% maturity (mm)")+
   theme(legend.position = "bottom", legend.direction = "horizontal",
         axis.text = element_text(size = 12),
         axis.title = element_text(size = 12),
         legend.text = element_text(size = 12))
 
 ggsave("./Maturity data processing/Doc/snow_SAM_comparison.png", width = 8, height = 6, units = "in")
 

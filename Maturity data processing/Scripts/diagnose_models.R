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


snow.chela <-  read.csv("./Maturity data processing/Data/snow_tanner_cheladatabase.csv") %>% #already filtered appropriately
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
  # MATURE = case_when((SIZE <=35) ~ 0,
  #                    (SIZE >= 135) ~ 1,
  #TRUE ~ MATURE)) %>%
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
         MATURE = case_when((SIZE <=35) ~ 0,
                            (SIZE >= 135) ~ 1,
                            TRUE ~ MATURE)) %>%
  as.data.frame(.) %>%
  rename(SIZE_5MM = SIZE_BINNED) %>%
  dplyr::select(YEAR, YEAR_F, YEAR_SCALED, STATION_ID, LATITUDE, LONGITUDE, SIZE_5MM, SIZE_CATEGORY, MATURE) %>%
  filter(YEAR != 2012)



# Models
legacy.mods <- list.files("./Maturity research/Output/Legacy models/")
legacy.dat <- read.csv("./Maturity data processing/Data/legacy_matmale_ratio.csv") 
gam.mod <- readRDS("./Maturity research/Output/GAM/gam_pooled.rda")
sdmTMB.mod.90 <- readRDS("./Maturity research/Output/sdmTMB/s(SIZE)_iid_90_sdmTMB.rda")
sdmTMB.mod.90.spvar <- readRDS("./Maturity research/Output/sdmTMB/s(SIZE)_spvar_iid_90_sdmTMB.rda")
sdmTMB.mod.iid.200 <- readRDS(paste0(remote_dir, "/SNOW/sdmTMB/s(SIZE, k=13)_iid_200_sdmTMB.rda"))
sdmTMB.mod.ar1.200 <- readRDS(paste0(remote_dir, "/SNOW/sdmTMB/s(SIZE)_ar1_200_sdmTMB.rda"))
tinyVAST.mod <- readRDS("./Maturity research/Output/tinyVAST/tinyVAST_fullsize_correlation.rda")
sdmTMB.mod.200.sizevar <- readRDS(paste0(remote_dir, "/SNOW/sdmTMB/s(SIZE)_sizevar_iid_200_sdmTMB.rda"))
# Residuals ----
# Legacy ----
legacy.resids <- data.frame()
for(ii in 1:length(legacy.mods)){
  
  mod <- readRDS(paste0("./Maturity research/Output/Legacy models/", legacy.mods[ii]))
  rr <- nlstools::nlsResiduals(mod)$resi2 # standardized residuals
  ff <- fitted(mod)
  
  ll <- logLik(mod)
  
  obs.dat <- legacy.dat %>% 
                  filter(YEAR == years[ii]) %>% 
                  cbind(., rr) %>%
                  cbind(., ll) %>%
                rename(observed = PROP_MATURE,
                       fitted = `Fitted values`,
                       resids = `Standardized residuals`)
    
  
  legacy.resids <- rbind(legacy.resids, obs.dat)
  
}

legacy.resids %>%
  group_by(YEAR) %>%
  arrange(resids, .by_group = TRUE) %>%
  mutate(n = n(),
         theoretical = qnorm(ppoints(n)),
         ordered_resid = resids) %>%
  ungroup() -> quantile.resid

ggplot()+
  theme_bw()+
  geom_point(quantile.resid, mapping = aes(theoretical, ordered_resid), size = 1, fill = "black")+
  geom_abline(slope = 1, intercept = 0, color = "red", linewidth = 1)+
  ylab("observed")+
  xlab("expected")+
  facet_wrap(~YEAR)+
  # scale_x_continuous(breaks = c(0, 0.5, 1))+
  # scale_y_continuous(breaks = c(0, 0.5, 1))+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12)) +
  ggtitle("Legacy (NLS)")

ggsave("./Maturity research/Figures/legacy_QQ.png", width = 9, height = 7, units = "in")


# GAM ----
dat <- snow_chela
resids <- simulateResiduals(gam.mod, n = 300)

dat$resids <- resids$scaledResiduals

# Extract residuals for custom plotting
gam_resids <- dat %>%
  group_by(YEAR) %>%
  arrange(resids, .by_group = TRUE) %>%
  mutate(n = n(),
         expected = ppoints(n), ,
         observed = resids) %>%
  ungroup()

#  QQ plot with ggplot2
ggplot()+
  theme_bw()+
  geom_point(gam_resids, mapping = aes(expected, observed), size = 1, fill = "black")+ #theoretical uniform quantiles vs. empirical residual quantiles
  geom_abline(slope = 1, intercept = 0, color = "red", linewidth = 1)+
  ylab("observed")+
  xlab("expected")+
  facet_wrap(~YEAR)+
  scale_x_continuous(breaks = c(0, 0.5, 1))+
  scale_y_continuous(breaks = c(0, 0.5, 1))+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12)) +
  ggtitle("GAM (pooled)")

ggsave("./Maturity research/Figures/pooledGAM_QQ.png", width = 9, height = 7, units = "in")



# sdmTMB iid 90 ----
dat <- snow_chela
resids <- simulate(sdmTMB.mod.90, nsim = 300, type= "mle-mvn")|>
            dharma_residuals(sdmTMB.mod.90, return_DHARMa = TRUE, plot = FALSE)

dat$DHARMa_resid <- resids$scaledResiduals

sdmTMB_resids_90 <- dat %>%
  group_by(YEAR) %>%
  arrange(DHARMa_resid, .by_group = TRUE) %>%
  mutate(
    n = n(),
    expected = ppoints(n),         # uniform quantiles
    observed = sort(DHARMa_resid)  # sort residuals for QQ
  ) %>%
  ungroup()

#  QQ plot with ggplot2
ggplot()+
  theme_bw()+
  geom_point(sdmTMB_resids_90, mapping = aes(expected, observed), size = 1, fill = "black")+ #theoretical uniform quantiles vs. empirical residual quantiles
  geom_abline(slope = 1, intercept = 0, color = "red", linewidth = 1)+
  ylab("observed")+
  xlab("expected")+
  facet_wrap(~YEAR)+
  scale_x_continuous(breaks = c(0, 0.5, 1))+
  scale_y_continuous(breaks = c(0, 0.5, 1))+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12)) +
  ggtitle("sdmTMB (knots =90; spatial_varying = off)")

ggsave("./Maturity research/Figures/sdmTMB_90_QQ.png", width = 9, height = 7, units = "in")

# sdmTMB iid 200 ----
dat <- snow_chela
resids <- simulate(sdmTMB.mod.iid.200, nsim = 300, type= "mle-mvn")|>
  dharma_residuals(sdmTMB.mod.iid.200, return_DHARMa = TRUE)

dat <- cbind(dat, DHARMa_resid = resids$scaledResiduals)


sdmTMB_resids_200 <- dat %>%
  group_by(YEAR) %>%
  arrange(DHARMa_resid, .by_group = TRUE) %>%
  mutate(
    n = n(),
    expected = ppoints(n),         # uniform quantiles
    observed = sort(DHARMa_resid)  # sort residuals for QQ
  ) %>%
  ungroup()

#  QQ plot with ggplot2
ggplot()+
  theme_bw()+
  geom_point(sdmTMB_resids_200, mapping = aes(expected, observed), size = 1, fill = "black")+ #theoretical uniform quantiles vs. empirical residual quantiles
  geom_abline(slope = 1, intercept = 0, color = "red", linewidth = 1)+
  ylab("observed")+
  xlab("expected")+
  facet_wrap(~YEAR)+
  scale_x_continuous(breaks = c(0, 0.5, 1))+
  scale_y_continuous(breaks = c(0, 0.5, 1))+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12)) +
  ggtitle("sdmTMB (knots = 200; spatial_varying = off)")

ggsave("./Maturity research/Figures/sdmTMB_200_byYEARQQ.png", width = 9, height = 7, units = "in")


sdmTMB_resids_200 <- dat %>%
  group_by(SIZE_5MM) %>%
  arrange(DHARMa_resid, .by_group = TRUE) %>%
  mutate(
    n = n(),
    expected = ppoints(n),         # uniform quantiles
    observed = sort(DHARMa_resid)  # sort residuals for QQ
  ) %>%
  ungroup()

#  QQ plot with ggplot2
ggplot()+
  theme_bw()+
  geom_point(sdmTMB_resids_200, mapping = aes(expected, observed), size = 1, fill = "black")+ #theoretical uniform quantiles vs. empirical residual quantiles
  geom_abline(slope = 1, intercept = 0, color = "red", linewidth = 1)+
  ylab("observed")+
  xlab("expected")+
  facet_wrap(~SIZE_5MM)+
  scale_x_continuous(breaks = c(0, 0.5, 1))+
  scale_y_continuous(breaks = c(0, 0.5, 1))+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12)) +
  ggtitle("sdmTMB (knots = 200; spatial_varying = off)")

ggsave("./Maturity research/Figures/sdmTMB_200_bySIZEQQ.png", width = 9, height = 7, units = "in")


ggplot(dat, aes(LONGITUDE, LATITUDE, fill = DHARMa_resid))+
  geom_point(shape = 21, size = 1.75, stroke = NA)+
  facet_wrap(~YEAR)+
  scale_fill_gradient2(midpoint = 0.5)+
  theme_bw()

ggsave("./Maturity research/Figures/sdmTMB_200_iid_spatialDHARMabyYEAR.png", width = 9, height = 7, units = "in")

ggplot(dat, aes(LONGITUDE, LATITUDE, fill = DHARMa_resid))+
  geom_point(shape = 21, size = 1.75, stroke = NA)+
  facet_wrap(~SIZE_5MM)+
  scale_fill_gradient2(midpoint = 0.5)+
  theme_bw()

ggsave("./Maturity research/Figures/sdmTMB_200_iid_spatialDHARMabySIZE.png", width = 10, height = 7, units = "in")

# sdmTMB iid spatially varying effect of size ----
dat <- snow_chela
mod <- sdmTMB.mod.200.sizevar
resids <- simulate(mod, nsim = 300, type= "mle-mvn")|>
  dharma_residuals(mod, return_DHARMa = TRUE, plot = FALSE)

dat$DHARMa_resid <- resids$scaledResiduals

ggplot(dat, aes(LONGITUDE, LATITUDE, fill = DHARMa_resid))+
  geom_point(shape = 21, size = 2)+
  facet_wrap(~YEAR)+
  scale_fill_gradient2(midpoint = 0.5)+
  theme_bw()


resids <- dat %>%
  group_by(YEAR) %>%
  arrange(DHARMa_resid, .by_group = TRUE) %>%
  mutate(
    n = n(),
    expected = ppoints(n),         # uniform quantiles
    observed = sort(DHARMa_resid)  # sort residuals for QQ
  ) %>%
  ungroup()

#  QQ plot with ggplot2
ggplot()+
  theme_bw()+
  geom_point(resids, mapping = aes(expected, observed), size = 1, fill = "black")+ #theoretical uniform quantiles vs. empirical residual quantiles
  geom_abline(slope = 1, intercept = 0, color = "red", linewidth = 1)+
  ylab("observed")+
  xlab("expected")+
  facet_wrap(~YEAR)+
  scale_x_continuous(breaks = c(0, 0.5, 1))+
  scale_y_continuous(breaks = c(0, 0.5, 1))+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12)) +
  ggtitle("sdmTMB (200; iid; spatial_varying = SIZE)")

ggsave("./Maturity research/Figures/sdmTMB_200iidsizevar_QQ.png", width = 9, height = 7, units = "in")

# sdmTMB ar1 200 ----
dat <- snow_chela
resids <- simulate(sdmTMB.mod.ar1.200, nsim = 300, type= "mle-mvn")|>
  dharma_residuals(sdmTMB.mod.ar1.200, return_DHARMa = TRUE, plot = FALSE)

dat$DHARMa_resid <- resids$scaledResiduals

sdmTMB_resids_200ar1 <- dat %>%
  group_by(YEAR) %>%
  arrange(DHARMa_resid, .by_group = TRUE) %>%
  mutate(
    n = n(),
    expected = ppoints(n),         # uniform quantiles
    observed = sort(DHARMa_resid)  # sort residuals for QQ
  ) %>%
  ungroup()

#  QQ plot with ggplot2
ggplot()+
  theme_bw()+
  geom_point(sdmTMB_resids_200ar1, mapping = aes(expected, observed), size = 1, fill = "black")+ #theoretical uniform quantiles vs. empirical residual quantiles
  geom_abline(slope = 1, intercept = 0, color = "red", linewidth = 1)+
  ylab("observed")+
  xlab("expected")+
  facet_wrap(~YEAR)+
  scale_x_continuous(breaks = c(0, 0.5, 1))+
  scale_y_continuous(breaks = c(0, 0.5, 1))+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12)) +
  ggtitle("sdmTMB (knots = 200; spatial_varying = off)")

ggsave("./Maturity research/Figures/sdmTMB_ar1200_QQ.png", width = 9, height = 7, units = "in")


ggplot(dat, aes(LONGITUDE, LATITUDE, fill = DHARMa_resid))+
  geom_point(shape = 21, size = 2)+
  facet_wrap(~YEAR)+
  scale_fill_gradient2(midpoint = 0.5)+
  theme_bw()

# tinyVAST ----
y_ir <- replicate(n = 300, expr = tinyVAST.mod$obj$simulate()$y_i)
dat <- tinyVAST.mod$data
dh_obj <- DHARMa::createDHARMa(
  simulatedResponse = y_ir,
  observedResponse = dat$MATURE,  # use your actual column for the observed response
  fittedPredictedResponse = fitted(tinyVAST.mod),
  integerResponse = TRUE
)

dat$DHARMa_resid <- dh_obj$scaledResiduals

tinyVAST_resids <- dat %>%
  filter(YEAR != 2012) %>%
  group_by(YEAR) %>%
  arrange(DHARMa_resid, .by_group = TRUE) %>%
  mutate(
    n = n(),
    expected = ppoints(n),
    observed = sort(DHARMa_resid)
  ) %>%
  ungroup()

ggplot() +
  geom_point(data = tinyVAST_resids, mapping = aes(expected, observed), size = 1, color = "black") +
  geom_abline(slope = 1, intercept = 0, color = "red", linewidth = 1) +
  facet_wrap(~YEAR) +
  scale_x_continuous(breaks = c(0, 0.5, 1)) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  theme_bw() +
  ylab("observed") +
  xlab("expected") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12)) +
  ggtitle("tinyVAST")

ggsave("./Maturity research/Figures/tinyVAST_QQ.png", width = 9, height = 7, units = "in")



# Log-likelihood (all must be fit with "ML") ----
gam.LL <- logLik(gam.mod)
sdmTMB.90.LL <- logLik(sdmTMB.mod.90)
sdmTMB.200.LL <- logLik(sdmTMB.mod.200)
sdmTMB.90spvar.LL <- logLik(sdmTMB.mod.90.spvar)
tinyVAST.LL <- logLik(tinyVAST.mod)

# Join diagnostics with AIC
diagnostics <- data.frame(MODEL = c("GAM", "sdmTMB", "sdmTMB", "sdmTMB", "tinyVAST"),
           KNOTS = c("NA", "90", "200", "90", "90"),
           SIZE_BIN = c("5MM", "5MM", "5MM", "5MM", "10MM"),
           TERMS = c("s(SIZE_5MM) + s(YEAR) + ti(SIZE_5MM, YEAR) + s(LATITUDE, LONGITUDE)",
                     "s(SIZE_5MM, k = 4) + YEAR_F); spatiotemporal RF",
                     "s(SIZE_5MM, k = 4) + YEAR_F); spatiotemporal RF",
                     "s(SIZE_5MM, k = 4) + YEAR_F); spatiotemporal RF, spatially varying size effect",
                     "s(SIZE_10MM, k = 4) + s(YEAR, k = 4); intersize correlations"),
           AIC(gam.mod, sdmTMB.mod.90, sdmTMB.mod.200, sdmTMB.mod.90.spvar, tinyVAST.mod),
           LOGLIK = c(gam.LL, sdmTMB.90.LL, sdmTMB.200.LL, sdmTMB.90spvar.LL, tinyVAST.LL)) %>%
  arrange(AIC)

write.csv(diagnostics, "./Maturity research/Output/candidate_model_diagnostics.csv")



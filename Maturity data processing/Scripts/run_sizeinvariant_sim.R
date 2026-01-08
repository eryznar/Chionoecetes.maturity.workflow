# PURPOSE: to run a simulation to compare maturity (0,1) weighted by CPUE, log CPUE, or unweighted

# Author: Emily Ryznar

# NOTES:

# LOAD LIBS/PARAMS ---------------------------------------------------------------------------------------
source("./Scripts/load_libs_params.R")
set.seed(123)

yr <- 2019

# SNOW CRAB ----
# Get survey cpue by 1mm bin for ALL shell 2 males, not just chela msrd ----
snow_cpue <- calc_cpue(crab_data = readRDS("./Data/snow_survey_specimenEBS.rda"),
                       species = "SNOW",
                       years = years,
                       sex = "male",
                       shell_condition = "new_hardshell",
                       bin_1mm = TRUE) %>%
  dplyr::select(YEAR, STATION_ID, LATITUDE, LONGITUDE, SIZE_1MM, CPUE)

# Set up modeling data using chela database ----
chela.data <- snow_chela %>%
  mutate(SIZE_1MM = floor(SIZE),
         MAT_TEXT = case_when((MATURE == 1) ~ "Mature",
                              TRUE ~ "Immature")) %>% # could change this to 10mm or 5mm bins instead of 1mm
  right_join(snow_cpue, .) %>% # adding in haul-level CPUE for corresponding 1mm size bin
  mutate(log.CPUE = as.integer(round(log(CPUE+10))),
         CPUE = as.integer(round(CPUE))) %>%
  na.omit() %>%
  st_as_sf(., coords = c("LONGITUDE", "LATITUDE"), crs = in.crs) %>%
  st_transform(., map.crs) %>%
  dplyr::select(YEAR, STATION_ID, SIZE_1MM, CPUE, log.CPUE, MATURE, geometry) %>%
  dplyr::rename(year = YEAR, station = STATION_ID, size = SIZE_1MM, cpue = CPUE, log.cpue = log.CPUE, mature = MATURE) %>%
  #filter(year == yr) %>%
  # mutate(mature = case_when((size <40) ~ 0,
  #                           (size >= 115) ~ 1,
  #                           TRUE ~ mature)) %>%
  group_by(station) %>%
  mutate(N = n()) %>%
  ungroup() 

# Plot distribution of CPUE by size for 60-90mm across stations
chela.data %>%
  # group_by(size) %>%
  # reframe(mean.cpue = mean(cpue)) %>%
  filter(size >= 60, size <=90) -> cpue.dat

# cpue.dat %>%
#   filter(cpue > 40000) -> pp

ggplot()+
  geom_histogram(cpue.dat, mapping = aes(cpue), binwidth = 30)+
  facet_wrap(~size, scales = "free")+
  theme_bw()


# Size invariant simulation ----
sim.fun <- function(yr, chela.dat, n.sim, n.crab, biased.sampling, iter){
  
  chela.dat %>%
    filter(year == yr)  -> mod.dat
  
  # 1) fit binomial model of mature versus cpue
  true.mod <- gam(mature ~ cpue, mod.dat, family= binomial)
  
  true.beta0 <- coef(true.mod)[1]
  true.beta1 <- coef(true.mod)[2]
  
  # 2) simulate population by sampling CPUE values from survey data and then calculating pmolt and maturity status
  cpue.sim <- sample(mod.dat$cpue, n.sim, replace = TRUE) # simulate CPUE by sampling from the data
  
  pmolt.sim <- predict(true.mod, data.frame(cpue = cpue.sim), type = "response") # calculate pmolt
  
  mat.sim <- rbinom(n.sim, 1, pmolt.sim)
  
  sim.dat <- data.frame(pmolt = pmolt.sim, mature = mat.sim, cpue = cpue.sim)
  
  pars <- data.frame()
  preds <- data.frame()
  for(ii in 1:length(iter)){
    # 4) simulate survey sampling 
    if(biased.sampling == TRUE){
      samples <- sample(nrow(sim.dat), size = n.crab, replace = TRUE, prob = sim.dat$cpue) # which rows are sampled from sim.dat?
    } else{
      samples <- sample(nrow(sim.dat), size = n.crab, replace = TRUE) # which rows are sampled from sim.dat?
    }
    sample.data <- sim.dat[samples,] 
    sample.cpue <- sample.data$cpue
    
    # 5) fit and predict models
    wtd.mod <- gam(mature ~ cpue, sample.data, family = binomial, weights = cpue)
    wtd.beta0 <- coef(wtd.mod)[1]
    wtd.beta1 <- coef(wtd.mod)[2]
    
    unwtd.mod <- gam(mature ~ cpue, sample.data, family = binomial)
    unwtd.beta0 <- coef(unwtd.mod)[1]
    unwtd.beta1 <- coef(unwtd.mod)[2]
    
    # param df
    params <- suppressWarnings(data.frame(model = rep(c("weighted", "unweighted"), each = 2), par = rep(c("beta0", "beta1"), 2), 
                                          estimate = c(wtd.beta0, wtd.beta1, unwtd.beta0, unwtd.beta1),
                                          true.beta1 = true.beta1, true.beta0 = true.beta0))
    
    params$iteration <- iter[ii]
    
    pars <- rbind(pars, params)
    
    # pred df
    pred.dat <- data.frame(cpue = 0:20000)
    pred.dat$wtd.pmolt <- predict(wtd.mod, pred.dat, type = "response")
    pred.dat$unwtd.pmolt <- predict(unwtd.mod, pred.dat, type = "response")
    pred.dat$true.pmolt <- predict(true.mod, pred.dat, type = "response")
    pred.dat$iteration <- iter[ii]
    
    preds <- rbind(preds, pred.dat)
  }
 
  return(list(pars = pars, preds = preds))
}

# Run sim
sim.fun(2019, chela.data, 1000, 50, TRUE, 1:50) -> out



# Compare
true.vals <- c(beta0 = unique(out$pars$true.beta0), beta1 = unique(out$pars$true.beta1))
               
out$pars %>% 
  group_by(model, par) %>%
  reframe(
    mean_est = mean(estimate),
    bias = mean_est - true.vals[par],
    mse =  mean((estimate - true.vals[par])^2)
  ) %>%
  distinct()

out$preds %>%
  pivot_longer(cols = 2:3, names_to = "weights", values_to = "pmolt") %>%
  group_by(cpue, weights) %>%
  mutate(mean.pmolt = mean(pmolt)) %>%
  ungroup() -> plot.dat

ggplot()+
  geom_line(plot.dat, mapping = aes(cpue, pmolt, group = interaction(weights, iteration)), color = "darkgrey", alpha = 0.5)+
  facet_wrap(~weights, nrow = 2)+
  geom_line(plot.dat, mapping = aes(cpue, true.pmolt), color = "blue", inherit.aes = FALSE)+
  geom_line(plot.dat, mapping = aes(cpue, mean.pmolt), linetype = "dashed")+
  theme_bw()


legacy.dat <- read.csv("./Maturity data processing/Data/legacy_matmale_ratio.csv")

legacy.mods <- list()

years <- unique(legacy.dat$YEAR)
  
for(ii in 1:length(years)){
  dat <- legacy.dat %>% filter(YEAR == years[ii])
  
  mod <- nls(PROP_MATURE ~ (1/(1 + exp(-a*(SIZE_BIN - b)))),
             data = dat,
             start = list(a = 0.10, b = 60.0),
             na.action = na.omit, nls.control(maxiter = 5000))
  
  saveRDS(mod, paste0("./Maturity research/Output/Legacy models/legacyNLS_", years[ii], ".rda"))
}
  
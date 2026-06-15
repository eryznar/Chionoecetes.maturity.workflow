# PURPOSE: to run model DHARMa residual diagnostics for maturity models (they should already have passed sanity())

# Author: Emily Ryznar

# NOTES:

# LOAD LIBS/PARAMS ------
source("./Scripts/Sourced scripts/load_libs_params.R")

# DIAGOSTIC FUNCTION ----
plot.resids <- function(model, species){
  resids <- simulate(model, nsim = 300, type= "mle-mvn")|>
    dharma_residuals(model, return_DHARMa = TRUE)
  
  dat <- cbind(model$data, DHARMa_resid = resids$scaledResiduals) %>%
          mutate(
            BIN_5MM = cut_width(SIZE, width = 5, center = 2.5,
                                closed = "left", dig.lab = 4),
            BIN2 = BIN_5MM) %>%
          separate(BIN2, sep = ",", into = c("LOWER", "UPPER")) %>%
          mutate(LOWER = as.numeric(sub('.', '', LOWER)),
                 UPPER = as.numeric(gsub('.$', '', UPPER)),
                 SIZE_5MM = (UPPER + LOWER) / 2) %>%
          dplyr::select(!c(BIN_5MM, LOWER, UPPER))
  
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
    ungroup() 
  
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
    ggtitle(paste0(species)) -> by_yr
  
  #ggsave(paste0("./Figures/", species, "_QQ_year.png", width = 10, height =10))
  
 
  rr_size <- dat %>%
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
    ggtitle(ggtitle(paste0(species))) -> by_size
  
  #ggsave(paste0("./Figures/", species, "_QQ_size.png", width = 10, height =10))
  
  
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
    ggtitle(ggtitle(paste0(species))) -> by_yr_sp
  
  #ggsave(paste0("./Figures/", species, "_spatialDHARMa_byYEAR.png"), width = 10, height = 9)
  
  
  ggplot(dat2, aes(LONGITUDE, LATITUDE, fill = DHARMa_resid))+
    geom_point(shape = 21, size = 1.75, stroke = NA)+
    facet_wrap(~SIZE_5MM)+
    scale_fill_gradient2(midpoint = 0.5)+
    theme_bw() +
    theme(legend.position = "bottom",
      legend.direction = "horizontal",
      strip.text = element_text(size = 10)) +
    ggtitle(ggtitle(paste0(species)))-> by_size_sp
  
  #ggsave(paste0("./Figures/", species, "_spatialDHARMa_bySIZE.png"), width = 10, height = 10)
  
  return(list(by_yr = by_yr, by_size = by_size, by_yr_sp = by_yr_sp, by_size_sp = by_size_sp,
         rr_yr = rr_yr, rr_size = rr_size))
}

# SNOW RESIDUALS ----
# Load model
snow_mod <- readRDS("./Models/snow_maturity_sdmTMB.rda")

# Run function
plot.resids(snow_mod, "Snow") -> snow.dharma 

# Check residuals
snow.dharma$by_yr # do residuals follow 1-1 line by year?
snow.dharma$by_size # do residuals follow 1-1 line by size?

snow.dharma$by_yr_sp # do residuals look uniformly distributed in space and across years (no clustering)?
snow.dharma$by_yr_size # do residuals look uniformly distributed in space and across sizes (no clustering)?

# Save residual plots
ggsave(plot = snow.dharma$by_yr, "./Figures/snow_QQ_year.png", width = 10, height = 10)
ggsave(plot = snow.dharma$by_size, "./Figures/snow_QQ_size.png", width = 10, height = 10)
ggsave(plot = snow.dharma$by_yr_sp, "./Figures/snow_spatialDHARMa_year.png", width = 10, height = 10)
ggsave(plot = snow.dharma$by_size_sp, "./Figures/snow_spatialDHARMa_size.png", width = 10, height = 10)

# TANNER RESIDUALS ----
# Load model
tanner_mod <- readRDS("./Models/tanner_maturity_sdmTMB.rda")

# Run function
plot.resids(tanner_mod, "Tanner") -> tanner.dharma 

# Check residuuals
tanner.dharma$by_yr # do residuals follow 1-1 line by year?
tanner.dharma$by_size # do residuals follow 1-1 line by size?

tanner.dharma$by_yr_sp # do residuals look uniformly distributed in space and across years (no clustering)?
tanner.dharma$by_yr_size # do residuals look uniformly distributed in space and across sizes (no clustering)?

# Save residual plots
ggsave(plot = tanner.dharma$by_yr, "./Figures/tanner_QQ_year.png", width = 10, height = 10)
ggsave(plot = tanner.dharma$by_size, "./Figures/tanner_QQ_size.png", width = 10, height = 10)
ggsave(plot = tanner.dharma$by_yr_sp, "./Figures/tanner_spatialDHARMa_year.png", width = 10, height = 10)
ggsave(plot = tanner.dharma$by_size_sp, "./Figures/tanner_spatialDHARMa_size.png", width = 10, height = 10)

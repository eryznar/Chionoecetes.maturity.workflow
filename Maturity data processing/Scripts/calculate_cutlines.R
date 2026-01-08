# PURPOSE: to calculate minima between dominant distribution modes of opilio chela height. Minima will serve as 
# the distribution cutline to define morphometric maturity. 

# Author: Emily Ryznar, Jon Richar, Shannon Hennessey

# LOAD LIBS/PARAMS/DATA ----
source("./Maturity data processing/Scripts/1) load_libs_params.R")

#  BIN DATA ----
read.csv("./Maturity data processing/Data/snow_tanner_cheladatabase.csv") %>% 
  mutate(BIN = cut_width(LN_CW, width = 0.025, center = 0.0125, closed = "left", dig.lab = 4), # Subset data into size intervals at ln(CW) of 0.025
         BIN2 = BIN) %>%
  separate(BIN2, sep = ",", into = c("LOWER", "UPPER")) %>%
  mutate(LOWER = as.numeric(sub('.', '', LOWER)),
         UPPER = as.numeric(gsub('.$', '', UPPER)),
         MIDPOINT = (UPPER + LOWER)/2) -> bin.dat


# CALCULATE MINIMA ----
# Sequentially apply KDE to the ln-chela height data for each interval, and identify minima of the resulting
# density distributions to define maturity classes in a given interval, looping over species

min.dat <- data.frame()
#bandwidth <- 0.04 # can adjust this for kernal smoothing (lower = less smooth)
species <- c("TANNER", "SNOW")

for(s in 1:length(species)){
  
  # Filter species and truncate sizes based on Jon's code
  if(species[s] == "SNOW"){
    bin.dat.temp <- bin.dat %>% 
      filter(SPECIES == species[s],
             LN_CW >= 3.9 & LN_CW <= 4.6)
  }
  
  if(species[s] == "TANNER"){
    bin.dat.temp <- bin.dat %>% 
      filter(SPECIES == species[s],
             LN_CW >= 4.3 & LN_CW <= 4.95)
  }
  
  
  # Set bins
  bins <- sort(unique(bin.dat.temp$BIN))
  
  for(i in 1:length(bins)){
    # Filter data by bin of interest
    opt.dat <- bin.dat.temp %>%
      filter(BIN == bins[i])
    
    # Calculate density for chela heights within that bin
    if(species[s] == "TANNER"){
      d <- density(opt.dat$LN_CH, kernel = "gaussian")
    }
    
    if(species[s] == "SNOW"){
      d <- density(opt.dat$LN_CH, kernel = "gaussian", bw = 0.04)
    }
    
    
    # Identify local minima and maxima:
    local_min_max <- function(x){
      signs <- sign(diff(x))
      list(minima = which(diff(signs) > 0) + 1,
           maxima = which(diff(signs) < 0) + 1)
    }
    
    extrema <- local_min_max(d$y) # both local minima and maxima
    
    maxima <- d$x[extrema$maxima] # isolate maxima index locations
    d.max <- d$y[extrema$maxima] # identify density at maxima index locations
    
    
    # Pull out maxima locations with top two greatest densities
    max.df <- cbind(maxima, d.max) %>%
      as.data.frame() %>%
      slice_max(d.max, n = 2) 
    
    
    # Identify the two dominant modes based on maxima to calculate 
    # minimum outside distribution tails
    ints <- round(max.df$maxima, 1)
    ints <- ints[order(ints)]
    
    # Add buffer to identifying minimum around maxima
    if(species[s] == "TANNER"){
      ints[2] <- ints[1] + 0.3
    }
    
    if(species[s] == "SNOW"){
      ints[2] <- ints[1] + 0.2
    }
    
    
    # find minimum in between dominant modes
    min <- optimize(approxfun(d$x, d$y), interval = ints)$minimum
    
    # create dataframe for output
    out <- data.frame(SPECIES = species[s],
                      LOWER = ints[1], 
                      UPPER = ints[2],
                      MINIMUM = min,
                      MIDPOINT = opt.dat$MIDPOINT,
                      BIN = bins[i]) %>%
      distinct()
    
    min.dat <- rbind(min.dat, out)
  } # end bin loop
} # end species loop


# Write output for use in subsequent scripts
min.dat %>%
  dplyr::select(SPECIES, BIN, LOWER, UPPER, MIDPOINT, MINIMUM) %>%
  write.csv("./Maturity data processing/Output/chela_cutline_minima.csv", row.names = FALSE) 


# Plot to make sure minima were calculated correctly
plot.dat <- right_join(bin.dat %>% dplyr::select(-UPPER, -LOWER), 
                       min.dat)

ggplot() +
  geom_density(plot.dat %>% filter(SPECIES == "SNOW"), mapping = aes(LN_CH), linewidth = 1, bw = 0.04) +
  geom_vline(plot.dat %>% filter(SPECIES == "SNOW"), mapping = aes(xintercept = MINIMUM), 
             color = "blue", linetype = "dashed", linewidth = 1) +
  facet_wrap(~BIN) +
  ylab("Density") +
  xlab("ln(chela height)") +
  ggtitle("Snow minima")+
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))
ggsave("./Maturity data processing/Figures/snow_density_plots.png", width = 20, height = 20)

ggplot() +
  geom_density(plot.dat %>% filter(SPECIES == "TANNER"), mapping = aes(LN_CH), linewidth = 1) +
  geom_vline(plot.dat %>% filter(SPECIES == "TANNER"), mapping = aes(xintercept = MINIMUM), 
             color = "blue", linetype = "dashed", linewidth = 1) +
  facet_wrap(~BIN) +
  ylab("Density") +
  xlab("ln(chela height)") +
  ggtitle("Tanner minima")+
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12))

ggsave("./Maturity data processing/Figures/tanner_density_plots.png", width = 20, height = 20)


# Fit linear model to evaluate ln(CH) minima (imm/mat division) and ln(CW) bin midpoint
for(s in 1:length(species)){
  mod <- lm(MINIMUM ~ MIDPOINT, min.dat %>% filter(SPECIES == species[s]))
  prd <- predict(mod, data = min.dat %>% filter(SPECIES == species[s]), se.fit = TRUE, interval = "confidence", level = 0.95)
  min.dat2 <- cbind(min.dat %>% filter(SPECIES == species[s]) %>% dplyr::select(!c(LOWER, UPPER)), prd$fit)
  
  labs <- data.frame(x = 4, y = 3, lab = paste0("R-squared = ", round(summary(mod)$r.squared, 2), "\np < 0.001"))
  
  cut.plot <- ggplot() +
    geom_ribbon(min.dat2, mapping = aes(x = MIDPOINT, ymin = lwr, ymax = upr), alpha = 0.4, fill = "cadetblue") +
    geom_line(min.dat2, mapping = aes(MIDPOINT, fit), linewidth = 1, color = "cadetblue") +
    geom_point(min.dat2, mapping = aes(MIDPOINT, MINIMUM)) +
    # geom_density(plot.dat, aes(LN_CH), color = "red", linetype = "dashed", linewidth = 1) +
    annotate("text", x = 4.4, y = 3, label = paste0("\nR-squared = ", round(summary(mod)$r.squared, 2), "\np < 0.001")) +
    ylab("Cutline (ln(chela height))") +
    xlab("Bin (ln(carapace width))") +
    ggtitle(paste0(tolower(species[s]), " cutline"))+
    theme_bw() +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 12))
  print(cut.plot)
  
  ggsave(paste0("./Maturity data processing/Figures/cutline_lm_plot_", species[s], ".png"),
         cut.plot, width = 6, height = 5)
}


# Pull out cutline parameters by species
cutline.params <- min.dat %>%
                    group_by(SPECIES) %>%
                    mutate(BETA0 = coef(lm(MINIMUM ~ MIDPOINT))[1],
                           BETA1 = coef(lm(MINIMUM ~ MIDPOINT))[2]) %>%
                  dplyr::select(SPECIES, BETA0, BETA1) %>%
                  distinct()

write.csv(cutline.params, "./Maturity data processing/Output/cutline_parameters.csv")

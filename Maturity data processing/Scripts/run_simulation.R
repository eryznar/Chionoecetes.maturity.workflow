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

# Exploratory plots ----
plot.dat <- snow_chela %>%
  mutate(SIZE_1MM = floor(SIZE),
         MAT_TEXT = case_when((MATURE == 1) ~ "Mature",
                              TRUE ~ "Immature")) %>% # could change this to 10mm or 5mm bins instead of 1mm
  right_join(snow_cpue, .) %>% # adding in haul-level CPUE for corresponding 1mm size bin
  mutate(log.CPUE = as.integer(round(log(CPUE+10))),
         CPUE = as.integer(round(CPUE))) %>%
  na.omit() 

plot.dat %>%
  filter(SIZE_1MM > 40 & SIZE_1MM <=95) %>%
  group_by(YEAR, LATITUDE, LONGITUDE, STATION_ID) %>%
  reframe(total = n(),
          mat = sum(MATURE),
          prop_mat = mat/total,
          CPUE = sum(CPUE)) -> prop.dat


ggplot() +
  geom_point(prop.dat, mapping = aes(LONGITUDE, LATITUDE, fill = prop_mat, size = CPUE), shape = 21, color = "black")+
  scale_fill_gradient2(high = scales::muted("red"), low = scales::muted("blue"),
                       midpoint=  0.5, 
                       name = "Proportion mature")+
  facet_wrap(~YEAR)+
  theme_bw()+
  ggtitle("Snow crab (40-95mm CW)")+
  theme(strip.text = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12))

ggsave("./Figures/propmat_CPUE_maps.png", width = 20, height = 10, units = "in")

ggplot(prop.dat, aes(CPUE, prop_mat))+
  geom_point()+
  geom_smooth()


ggplot() +
  geom_point(prop.dat, mapping = aes(LONGITUDE, LATITUDE, fill = log(total+10), size = CPUE), shape = 21, color = "black")+
  scale_fill_gradient2(high = scales::muted("red"), low = scales::muted("blue"),
                       midpoint=  3.5,
                       name = "log(N+10)")+
  facet_wrap(~YEAR)+
  theme_bw()+
  theme(strip.text = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12))

ggsave("./Figures/samplesize_CPUE_maps.png", width = 20, height = 10, units = "in")


plot.dat %>%
  group_by(STATION_ID) %>%
  mutate(N = n()) %>%
  ungroup() %>%
  filter(N >10, YEAR == 2019) -> tt

ggplot(tt, aes(LN_CW, LN_CH, color = as.factor(MATURE)))+
  geom_point()+
  facet_wrap(~STATION_ID)



# SIMULATION ----
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
    mutate(mature = case_when((size <40) ~ 0,
                              (size >= 115) ~ 1,
                              TRUE ~ mature)) %>%
    group_by(station) %>%
    mutate(N = n()) %>%
    ungroup() 

 ggplot(chela.data, aes(cpue, ))

  # Format population (specimen data) ----
  snow.spec <- readRDS("./Data/snow_survey_specimenEBS.rda")$specimen %>%
    filter(HAUL_TYPE !=17, SEX == 1, SHELL_CONDITION == 2)
           #YEAR ==yr)  # make sure filter for males, sh2, not HT17
  
  population.data <- snow.spec %>%
    dplyr::select(YEAR, STATION_ID, SIZE_1MM) %>%
    rename(station = STATION_ID, size = SIZE_1MM) %>%
    right_join(snow_cpue %>% 
                 #filter(YEAR == yr) %>% 
                 dplyr::select(YEAR, STATION_ID, SIZE_1MM, CPUE) %>% 
                 rename(station = STATION_ID, size=  SIZE_1MM, cpue = CPUE), .) %>%
    filter(station %in% mod.dat$station) %>% # have to filter to stations with chela crab
    mutate(bin = case_when((size >=40 & size <=59) ~ "40-59mm",
                           (size >= 60 & size<=99) ~ "60-99mm",
                           (size >99) ~ ">99mm"),
           log.cpue = as.integer(round((log(cpue + 10)))),
           cpue = as.integer(round(cpue)))
  
  ggplot()

  # Set up simulation function ----
  sim_fun <- function(yr, chela.dat, pop.dat, n_samples, weights, by_block, n_blocks, n_reps){
    
    # 1) Fit "true" maturity curves using chela data either by block or universal ----
    if(by_block == FALSE){
      chela.dat <- chela.dat
      
      true.mod <- gam(mature ~ s(size), chela.dat, family = "binomial")
      true.curve <- data.frame(size= 0:150, true = predict(true.mod, data.frame(size= 0:150), type = "response"))
      
      ggplot(true.curve, aes(size, true))+
        geom_line(linewidth = 1)+
        theme_bw()+
        ylab("p(molt)")+
        ggtitle(paste0(yr, " 'true' ogives"))+
        theme(axis.text = element_text(size = 12),
              axis.title = element_text(size = 12),
              legend.text = element_text(size = 12),
              legend.title = element_text(size = 12))
      
      ggsave(paste0("./Figures/Simulation/", yr, "_true_ogives.png"), width = 10, height = 9)
     
    }else{
      # set up spatial blocking 
      chela.dat <- mod.dat
      
      coords <- st_coordinates(chela.dat)
      
      clust <- kmeans(coords, centers = n_blocks, nstart = 25)
      
      chela.dat$block.id <- clust$cluster
      
      chela.dat <- cbind(st_coordinates(chela.dat), st_drop_geometry(chela.dat))
      
      station.block <- chela.dat %>% dplyr::select(station, block.id) %>% group_by(block.id) %>% mutate(N_stations = length(unique(station)),
                                                                                                        N_chela = n()) %>% ungroup()
      
      # fit gams by spatial block
      gam.dat <-  chela.dat %>%
                  group_by(block.id) %>%
                  mutate(N.stations = n()) %>%
                  nest() %>%
                  mutate(true.mod = map(data, ~gam(mature ~ s(size), family = binomial, data = .x)), gam.tidy = map(true.mod, broom::tidy))
      
      # Generate predictions for "true" maturity ogives 
      pred.dat <- gam.dat %>%
        mutate(
          predictions = map2(
            true.mod,
            data,
            ~{newd = data.frame(size = 0:150)
            preds = predict(.x, newdata = newd, type = "response")
            tibble(size = 0:150, p.molt = preds)}))
      
      true.curve <- pred.dat %>%
        dplyr::select(block.id, predictions) %>%
        unnest(predictions) %>%
        right_join(., station.block, relationship = "many-to-many")
      
      # Plot and save predictions for "true" maturity ogives
      labs <- unique(paste0("Block ", true.curve$block.id, " (# stations=", true.curve$N_stations, ";# chela=", true.curve$N_chela,")"))
      names(labs) <- unique(true.curve$block.id)
      
      ggplot(true.curve, aes(size, p.molt))+
        geom_line(linewidth = 1)+
        facet_wrap(~block.id, labeller = as_labeller(labs))+
        theme_bw()+
        ylab("p(molt)")+
        ggtitle(paste0(yr, " 'true' ogives by spatial block (# blocks=", n_blocks, ")"))+
        theme(axis.text = element_text(size = 12),
              axis.title = element_text(size = 12),
              legend.text = element_text(size = 12),
              legend.title = element_text(size = 12))-> p.1
      
      ggplot(chela.dat, aes(X, Y, color = as.factor(block.id)))+
        geom_point(size = 2)+
        theme_bw()+
        theme(axis.text = element_text(size = 12),
              axis.title = element_text(size = 12),
              legend.text = element_text(size = 12),
              legend.title = element_text(size = 12)) -> p.2
      
      p.1/p.2
      
      ggsave(paste0("./Figures/Simulation/", yr, "_true_ogives_by_", n_blocks, "block.png"), width = 11, height = 9)
    
    }
    
    # 2) Add spatial blocks to population (specimen) data (should be = between chela and pop data) ----
    if(by_block == FALSE){
      pop.dat <- pop.dat
    }else{
      pop.dat <- pop.dat %>%
        right_join(., station.block, relationship = "many-to-many")
    }
   
    
    # 3) Apply maturity curves to population data to estimate pmolt and mature for all crab, not just chela measured ----
    if(by_block == "FALSE"){ 
      p.molt <- predict(true.mod, pop.dat, type = "response")
      
      pop.dat2 <- cbind(pop.dat, p.molt) %>%
        mutate(mature = rbinom(n = nrow(.), size = 1, prob = p.molt),
               mature = case_when((size <40) ~ 0,
                                       (size >= 115) ~ 1,
                                       TRUE ~ mature))
    }else{ # by_block ==TRUE
      pop.dat2 <- gam.dat %>%
        group_by(block.id) %>%
        mutate(predictions = map2(true.mod, block.id,
                                  ~{# Subset pop.dat by block
                                    newd <- pop.dat %>% filter(block.id == .y)
                                    
                                    # Predict pmolt
                                    preds <- predict(.x, newdata = newd, type = "response")
                                    
                                    # Simulate mature status from pmolt
                                    mature_sim <- rbinom(n = length(preds), size = 1, prob = preds)
                                    
                                    # Modify mature smaller and larger than certain cutoffs
                                    mature_sim <- case_when(
                                      newd$size < 40     ~ 0,
                                      newd$size >= 115   ~ 1,
                                      TRUE               ~ mature_sim)

                                    
                                    # Combine
                                    bind_cols(
                                      newd %>% dplyr::select(!block.id),
                                      tibble(
                                        p.molt = preds, mature = mature_sim))})) %>%
        dplyr::select(block.id, predictions) %>%
        unnest(predictions) 
    }
    
    # 4) Run simulation ----
    sim.out <- data.frame()

    for (ii in 1:n_reps){

      # Simulate sampling from the population n_samples number of times
      if(weights == "cpue" & by_block == FALSE){
        # Simulate cpue-dependent sampling
        samples <- pop.dat[sample(nrow(pop.dat), n_samples, prob = pop.dat$cpue),]
        
      } else if(weights == "log.cpue" & by_block == FALSE){
        # Simulate cpue-dependent sampling
        samples <- pop.dat2 %>%
          slice_sample(n = n_samples, weight_by = log.cpue)
        
      } else if(weights == "unweighted" & by_block == FALSE){
        # Simulate unweighted sampling
        samples <-  pop.dat2 %>%
          slice_sample(n = n_samples)
        
      } else if(weights == "cpue" & by_block == TRUE){
        samples <- pop.dat2 %>%
          group_by(block.id) %>%
          slice_sample(n = n_samples, weight_by = cpue)
        
      } else if(weights == "log.cpue" & by_block == TRUE){
        samples <- pop.dat2 %>%
          group_by(block.id) %>%
          slice_sample(n = n_samples, weight_by = log.cpue)
        
      } else{
        samples <- pop.dat2 %>%
          group_by(block.id) %>%
          slice_sample(n = n_samples)
      }

      # Fit models based on population samples
      if(by_block == FALSE){
        cpue.mod <- gam(mature ~ s(size), data = samples, weights = cpue, family = binomial)
        log.cpue.mod <- gam(mature ~ s(size), data = samples, weights = log.cpue, family = binomial)
        unwtd.mod <- gam(mature ~ s(size), data = samples, family = binomial)

        # Extract model parameters
        sample.params <- cbind(size = 0:150, cpue = predict(cpue.mod, data.frame(size = 0:150), type = "response"),
                               log.cpue = predict(log.cpue.mod, data.frame(size = 0:150), type = "response"),
                               unwtd = predict(unwtd.mod, data.frame(size = 0:150), type = "response")) %>%
          as.data.frame() %>%
          pivot_longer(., cols = 2:4, names_to = "weights", values_to = "preds") %>%
          mutate(iter = (1:n_reps)[ii])

      } else{ # by_block == TRUE
        # Prepare a weighting methods table
        weighting_methods <- tibble::tibble(
          method = c("cpue", "log.cpue", "unweighted"),
          weight_var = c("cpue", "log.cpue", NA)
        )
        
        # Expand data frame to have one row per block × method
        expand_blocks <- samples %>%
          group_by(block.id) %>%
          nest() %>%
          crossing(weighting_methods)
        
        # Fit models by block and weighting method
        sim.gam <- expand_blocks %>%
          mutate(
            model = pmap(list(data, weight_var), ~{
              if (is.na(..2)) {
                gam(mature ~ s(size), family = binomial, data = ..1)
              } else {
                gam(mature ~ s(size), family = binomial, data = ..1, weights = ..1[[..2]])
              }
            }),
            model_tidy = map(model, broom::tidy)
          )
        
        # Predict models by block and weighting method
        sim.gam2 <- sim.gam %>%
          mutate(
            pred = map(model, ~{
              p.molt = predict(.x, newdata = data.frame(size= 0:150), type = "response")
              bind_cols(data.frame(size= 0:150), tibble(p.molt = p.molt))
            })
          )
        
        # Extract model parameters
        sample.params <- sim.gam2 %>%
          dplyr::select(block.id, method, pred) %>%
          unnest(pred) %>%
          right_join(., station.block %>% dplyr::select(!station), relationship = "many-to-many") %>%
          mutate(iter = (1:n_reps)[ii]) %>%
          rename(weights= method, preds = p.molt)

      }

      # Bind across iterations
      sim.out <- rbind(sim.out, sample.params)
     }
    
    return(list(true.curve = true.curve, sim.out = sim.out))
  }
  
  # Run simulation by_block----
  set.seed(123)
  chela.dat = chela.data
  pop.dat = population.data
  n_samples= 200
  weights= "unweighted"
  by_block = TRUE
  n_blocks = 10
  n_reps = 10
  
  sim_fun(yr, chela.dat, pop.dat, n_samples, weights, by_block, n_blocks, n_reps) -> out
  
  # Compare simulation output to "true" values
  out$sim.out %>%
    distinct() %>%
    right_join(., out$true.curve %>% rename(true = p.molt) %>% dplyr::select(!station), 
               relationship = "many-to-many", by = c("size", "block.id", "N_chela", "N_stations")) %>%
    group_by(size, weights, block.id) %>%
    mutate(mean.preds = mean(preds)) %>%
    ungroup() -> sum.dat
  
  ggplot()+
    geom_line(sum.dat, mapping = aes(size, preds, group = interaction(weights, iter)), color = "grey", alpha = 0.75) +
    geom_line(sum.dat, mapping = aes(size, mean.preds), color = "black", alpha = 0.75, linetype = "dashed") +
    geom_line(sum.dat, mapping = aes(size, true),  color = "blue") +
    facet_grid(block.id~weights)+
    ggtitle(paste0("P(molt) by block (n=", n_samples, ")"))+
    theme_bw()+
    theme(strip.text = element_text(size = 12),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12))
  
  ggsave(paste0("./Figures/Simulation/", yr, "_", n_samples, "N_sim_curves_spatial_blocks.png"), width = 10, height= 9)
  
  sum.dat %>%
    group_by(block.id, weights) %>%
    reframe(mse = mean(preds-true)^2)
  
  # Run simulation ebs-wide----
  set.seed(123)
  chela.dat = chela.data
  pop.dat = population.data
  n_samples= 200
  weights= "unweighted"
  by_block = FALSE
  n_blocks = 10
  n_reps = 10
  
  sim_fun(yr, chela.dat, pop.dat, n_samples, weights, by_block, n_blocks, n_reps) -> out
  
  # Compare simulation output to "true" values 
  out$sim.out %>%
    distinct() %>%
    right_join(., out$true.curve) %>%
    group_by(size, weights) %>%
    mutate(mean.preds = mean(preds)) %>%
    ungroup() -> sum.dat
  
  ggplot()+
    geom_line(sum.dat, mapping = aes(size, preds, group = interaction(weights, iter)), linewidth = 1, color = "grey", alpha = 0.75) +
    geom_line(sum.dat, mapping = aes(size, mean.preds), linewidth = 1, color = "black", alpha = 0.75, linetype = "dashed") +
    geom_line(sum.dat, mapping = aes(size, true), linewidth = 1, color = "blue") +
    facet_wrap(~weights, nrow = 3)+
    theme_bw()+
    ggtitle(paste0("P(molt) (n=", n_samples, ")"))+
    theme(strip.text = element_text(size = 12),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12))
  
  ggsave(paste0("./Figures/Simulation/", yr, "_", n_samples, "_sim_curves_no_blocks.png"), width = 10, height= 9)
  
  
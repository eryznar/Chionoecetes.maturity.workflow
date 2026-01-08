# 5mm Size only (LOWER AIC) ----
# spacetime_term: link (->), lag, param_name
sizes <- sort(unique(mod.dat$size_5mm))
spacetime_term <- "pmolt <-> pmolt, 0, sd_spacetime"

# mesh
mat.mesh = fmesher::fm_mesh_2d(mod.dat[, c("lon", "lat")],cutoff = 15)

# dist list to assign distribution to each size
dist <- unique(mod.dat$mature) %>%
  set_names() %>%
  map(~binomial())

# add year-size interaction
mod.dat <- mod.dat %>%
  mutate(size_f = factor(paste0("size", size_5mm)),
         year_size = interaction(year, size_f))

# specify controls
control <- tinyVASTcontrol(getsd = FALSE,
                           profile = c("alpha_j", "alpha2_j"), # Spatial, temporal, or spatio-temporal path coefficient for dynamic SEM effects, representing interactions or dependencies among variables.
                           trace = 1, # do you want messaging?
                           calculate_deviance_explained = FALSE) 


# 2) Fit model 
mod.dat2 <- droplevels(mod.dat) %>%
  mutate(year = as.numeric(as.character(year)),
         size_category = as.factor(size_category),
         #size_5mm = as.factor(size_5mm),
         var= "pmolt")

mat.mod <- tinyVAST(
  data = mod.dat2,
  formula = mature ~ 0 + year_size,
  time_term = NULL,   # estimates constant variance across time    
  space_term = NULL,  # estimates constant variance across space
  spacetime_term = spacetime_term,  # does variance change in both space and time?
  family = binomial(),
  #delta_options = list(delta_formula = ~ 0 + year_size), # can only use this in when a delta family is used
  space_columns = c("lon", "lat"),
  variable_column = "var",
  time_column = "year",
  #distribution_column = "mature",
  spatial_domain = mat.mesh,
  control = control
)

# saveRDS(mat.mod, paste0(dir, "/SNOW/Model-based indices/Models/snow_maturity_tinyVAST_5mm_sizeonly.rda"))
# 
# mat.mod <- readRDS(paste0(dir, "/SNOW/Model-based indices/Models/snow_maturity_tinyVAST_5mm_sizeonly.rda"))

# 3) Predict model
ebs_grid2 <- ebs_grid %>%
  dplyr::select(area_km2, X, Y) %>%
  replicate_df(., "year", unique(mod.dat2$year)) %>%
  rename(lon = X, lat = Y) %>%
  replicate_df(., "size_5mm", unique(mod.dat2$size_5mm)) %>%
  mutate(year_size = paste0(year, ".size", size_5mm)) %>%
  filter(year_size %in% mod.dat2$year_size)

prds <- predict(mat.mod, ebs_grid2, type = "response", se_fit = TRUE)

pred.dat <- cbind(prds, ebs_grid2)

ogive.dat <- pred.dat %>%
  group_by(year, size_5mm) %>%
  reframe(pmolt = mean(prds)) %>%
  mutate(year = as.factor(year),
         size_5mm = as.numeric(as.character(size_5mm))) %>%
  distinct()

ggplot(ogive.dat, aes(size_5mm, pmolt, color = year))+
  geom_line(linewidth = 1.25)+
  facet_wrap(~year)+
  ggtitle("Snow ogives (5mm bins; size only)")+
  xlab("Carapace width (mm)")+
  scale_color_viridis_d(option = "mako", name = "Year")+
  theme_bw()+
  scale_x_continuous(breaks = seq(0, 150, by = 25))

# space + spacetime ----
# spacetime_term: link (->), lag, param_name
sizes <- sort(unique(mod.dat$size_5mm))
time_term <- "\n  "
for(ii in 1:length(sizes)) {
  time_term  <- c(paste0(time_term ,  sizes[ii], " -> ", sizes[ii], ", 1, space_sd_", sizes[ii], "\n"))
}

spacetime_term <- "\n  "
for(ii in 1:length(sizes)) {
  spacetime_term  <- c(paste0(spacetime_term ,  sizes[ii], " -> ", sizes[ii], ", 1, spacetime_sd_", sizes[ii], "\n"))
}

st_term <- c(time_term, spacetime_term)

# mesh
mat.mesh = fmesher::fm_mesh_2d(mod.dat[, c("lon", "lat")],cutoff = 15)

# dist list to assign distribution to each size
dist <- unique(mod.dat$mature) %>%
  set_names() %>%
  map(~binomial())

# add year-size interaction
mod.dat <- mod.dat %>%
  mutate(size_f = factor(paste0("size", size_5mm)),
         year_size = interaction(year, size_f))

# specify controls
control <- tinyVASTcontrol(getsd = FALSE,
                           profile = c("alpha_j", "alpha2_j"), # Spatial, temporal, or spatio-temporal path coefficient for dynamic SEM effects, representing interactions or dependencies among variables.
                           trace = 1, # do you want messaging?
                           calculate_deviance_explained = FALSE) 


# 2) Fit model 
mod.dat2 <- droplevels(mod.dat) %>%
  mutate(year = as.numeric(as.character(year)),
         size_category = as.factor(size_category),
         size_5mm = as.factor(size_5mm),
         var= "pmolt")

mat.mod2 <- tinyVAST(
  data = mod.dat2,
  formula = mature ~ 0 + size_5mm+year,
  time_term = NULL,   # estimates constant variance across time    
  space_term = NULL,  # estimates constant variance across space
  spacetime_term = st_term,  # does variance change in both space and time?
  family = binomial(),
  #delta_options = list(delta_formula = ~ 0 + year_size), # can only use this in when a delta family is used
  space_columns = c("lon", "lat"),
  variable_column = "size_5mm",
  time_column = "year",
  #distribution_column = "mature",
  spatial_domain = mat.mesh,
  control = control
)

saveRDS(mat.mod2, paste0(dir, "/SNOW/Model-based indices/Models/snow_maturity_tinyVAST_space-spacetime.rda"))
# 
# mat.mod <- readRDS(paste0(dir, "/SNOW/Model-based indices/Models/snow_maturity_tinyVAST_5mm_sizeonly.rda"))

# 3) Predict model
ebs_grid2 <- ebs_grid %>%
  dplyr::select(area_km2, X, Y) %>%
  replicate_df(., "year", unique(mod.dat2$year)) %>%
  rename(lon = X, lat = Y) %>%
  replicate_df(., "size_5mm", unique(mod.dat2$size_5mm)) %>%
  mutate(year_size = paste0(year, ".size", size_5mm)) %>%
  filter(year_size %in% mod.dat2$year_size)

prds <- predict(mat.mod, ebs_grid2, type = "response", se_fit = TRUE)

pred.dat <- cbind(prds, ebs_grid2)

ogive.dat <- pred.dat %>%
  group_by(year, size_5mm) %>%
  reframe(pmolt = mean(prds)) %>%
  mutate(year = as.factor(year),
         size_5mm = as.numeric(as.character(size_5mm))) %>%
  distinct()

ggplot(ogive.dat, aes(size_5mm, pmolt, color = year))+
  geom_line(linewidth = 1.25)+
  facet_wrap(~year)+
  ggtitle("Snow ogives (5mm bins; size only)")+
  xlab("Carapace width (mm)")+
  scale_color_viridis_d(option = "mako", name = "Year")+
  theme_bw()+
  scale_x_continuous(breaks = seq(0, 150, by = 25))




# spacetime ----
# spacetime_term: link (->), lag, param_name
sizes <- sort(unique(mod.dat$size_5mm))
spacetime_term <- "\n  "
for(ii in 1:length(sizes)) {
  spacetime_term  <- c(paste0(spacetime_term ,  sizes[ii], " -> ", sizes[ii], ", 1, spacetime_sd_", sizes[ii], "\n"))
}

st_term <- c(spacetime_term)

# mesh
mat.mesh = fmesher::fm_mesh_2d(mod.dat[, c("lon", "lat")],cutoff = 15)

# dist list to assign distribution to each size
dist <- unique(mod.dat$mature) %>%
  set_names() %>%
  map(~binomial())

# add year-size interaction
mod.dat <- mod.dat %>%
  mutate(size_f = factor(paste0("size", size_5mm)),
         year_size = interaction(year, size_f))

# specify controls
control <- tinyVASTcontrol(getsd = FALSE,
                           profile = c("alpha_j", "alpha2_j"), # Spatial, temporal, or spatio-temporal path coefficient for dynamic SEM effects, representing interactions or dependencies among variables.
                           trace = 1, # do you want messaging?
                           calculate_deviance_explained = FALSE) 


# 2) Fit model 
mod.dat2 <- droplevels(mod.dat) %>%
  mutate(year = as.numeric(as.character(year)),
         size_category = as.factor(size_category),
         size_5mm = as.factor(size_5mm),
         var= "pmolt")
































































































































# mat.mod <- readRDS(paste0(dir, "/SNOW/Model-based indices/Models/snow_maturity_tinyVAST_5mm_sizeonly.rda"))

# 3) Predict model
ebs_grid2 <- ebs_grid %>%
  dplyr::select(area_km2, X, Y) %>%
  replicate_df(., "year", unique(mod.dat2$year)) %>%
  rename(lon = X, lat = Y) %>%
  replicate_df(., "size_5mm", unique(mod.dat2$size_5mm)) %>%
  mutate(year_size = paste0(year, ".size", size_5mm)) %>%
  filter(year_size %in% mod.dat2$year_size)

prds <- predict(mat.mod, ebs_grid2, type = "response", se_fit = TRUE)

pred.dat <- cbind(prds, ebs_grid2)

ogive.dat <- pred.dat %>%
  group_by(year, size_5mm) %>%
  reframe(pmolt = mean(prds)) %>%
  mutate(year = as.factor(year),
         size_5mm = as.numeric(as.character(size_5mm))) %>%
  distinct()

ggplot(ogive.dat, aes(size_5mm, pmolt, color = year))+
  geom_line(linewidth = 1.25)+
  facet_wrap(~year)+
  ggtitle("Snow ogives (5mm bins; size only)")+
  xlab("Carapace width (mm)")+
  scale_color_viridis_d(option = "mako", name = "Year")+
  theme_bw()+
  scale_x_continuous(breaks = seq(0, 150, by = 25))




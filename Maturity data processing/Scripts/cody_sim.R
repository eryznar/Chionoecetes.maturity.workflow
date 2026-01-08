source("./Scripts/load_libs_params.R")

set.seed(123)
true_b0 <- -4
true_b1 <- 0.5
size_bins <- seq(50, 150, by=10)
n_pop <- 10000

# True maturity curve
population <- data.frame(
  size = sample(size_bins, n_pop, replace=TRUE),
  cpue = NA
)
population$mat_prob <- plogis(true_b0 + true_b1*population$size/10)
population$mature <- rbinom(n_pop, 1, population$mat_prob)

# Simulate CPUE weighting (size-dependent sampling)
population$cpue <- 2*(population$size/100)^1.5  # Gamma=1.5 bias

# Biased sampling based on CPUE
sample_size <- 500
sampled_indices <- sample(n_pop, sample_size, prob=population$cpue)
cpue_sample <- population[sampled_indices, ]

# Unweighted comparison sample
unweighted_sample <- population[sample(n_pop, sample_size), ]

# Fit models ----
# Weighted model
cpue_model <- glm(mature ~ I(size/10), 
                  data = cpue_sample,
                  family = binomial,
                  weights = as.integer(round(cpue)))  # CPUE as weights

# Unweighted model
unweighted_model <- glm(mature ~ I(size/10), 
                        data = unweighted_sample,
                        family = binomial)



# Simulation wrapper function
run_simulation <- function(n_sims=500) {
  results <- list()
  for(i in 1:n_sims) {
    # Resample data
    cpue_sample <- population[sample(n_pop, sample_size, prob=population$cpue), ]
    
    # Fit models
    mod_weighted <- glm(mature ~ I(size/10), data=cpue_sample, 
                        family=binomial, weights=as.integer(round(cpue)))
    mod_unweighted <- glm(mature ~ I(size/10), data=cpue_sample, # should this be unweighted???????
                          family=binomial)
    
    # Store results
    results[[i]] <- data.frame(
      sim = i,
      model = rep(c("weighted","unweighted"), each=2),
      term = rep(c("intercept","slope"), 2),
      estimate = c(coef(mod_weighted), coef(mod_unweighted))
    )
  }
  return(bind_rows(results))
}

# Analyze simulation results
sim_results <- run_simulation() %>% 
  group_by(model, term) %>% 
  reframe(
    mean_est = mean(estimate),
    bias = mean_est - c(true_b0, true_b1)[as.numeric(factor(term))],
    mse = mean((estimate - c(true_b0, true_b1)[as.numeric(factor(term))])^2)
  )

distinct(sim_results)

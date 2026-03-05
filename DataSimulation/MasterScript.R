# Data generation - Master Script
nsim <- 100
N_sample_bootstrap <- 100
treatment_effect <- 2
n_covariates <- 8

nstudies <- 10 
tau_effect <- c(0.2, 0.5) # treatment heterogeneity 


set.seed(42)
# Generate beta from a Uniform(-2, 2)
beta <- round(runif(n_covariates, -2, 2), 2)
# main effects of covariates
effect_modification <- c("SMALL", "LARGE")
# From a uniform (-1,1)
gamma1 <- round(runif(n_covariates, -1.5, 1.5), 2)
# From a uniform(-2, 2)
gamma2 <- round(runif(n_covariates, -3, 3), 2)

sample_sizes <- c("SMALL", "LARGE")
sample_size_min <- c(50, 300)
sample_size_max <- c(200, 600)

complexity <- c(T, F)

parameters <- expand.grid(heterogeneity = tau_effect, 
                          effect_modification = effect_modification, 
                          sample_size = sample_sizes, 
                          complex_relationship = complexity)

nscenarios <- nrow(parameters)

data_training <- vector("list", length = nscenarios)
data_testing <- vector("list", length = nscenarios)

set.seed(42)
for(i in 1:nscenarios){
  train_temp <- vector("list", length = nsim); test_temp <- vector("list", length = nsim)
  if(parameters[i, "effect_modification"] == "SMALL"){
    gammas <- gamma1
  } else {
    gammas <- gamma2
  }
  
  if(parameters[i, "sample_size"] == "SMALL"){
    smin <- sample_size_min[1]; smax <- sample_size_max[1]
  } else {
    smin <- sample_size_min[2]; smax <- sample_size_max[2]
  }
  
  for(j in 1:nsim){
    obj <- generate_cont_data(nstudies = nstudies, 
                              tau_effect = parameters[i, "heterogeneity"], # treatment heterogeneity 
                              beta = beta, # main effects of covariates
                              gamma = gammas, # effects of interactions
                              treatment_effect = treatment_effect, # treatment effect
                              ss_min = smin,
                              ss_max = smax,
                              complexity = parameters[i, "complex_relationship"])
    
    train_temp[[j]] <- split(obj$data, obj$data$ID_study)
    test_temp[[j]] <- obj$true_ITE
  }
  
  data_training[[i]] <- train_temp; data_testing[[i]] <- test_temp
}






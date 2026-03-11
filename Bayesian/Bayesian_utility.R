# Bayesian utility

# Function for Bayesian Methods

BayesianLM <- function(studies_data, n.covariates,
                       nstudies){
  Model.ITE <- "
    model{
      for(i in 1:n_patients){
		    y[i] ~ dnorm(mu[i], prec)
		    mu[i] = alpha + inprod(beta[], x[i, ]) + delta * trt[i] +
		    inprod(gamma[], x_int[i, ]) 	     
      }
	    
	    prec <- pow(sigma, -2)
    	sigma ~ dunif(0, 10)

	    alpha ~ dnorm(0, 0.001)

	    for(j in 1:N.covariates){ 
		    beta[j] ~ dnorm(0, 0.001)
	    }
	    
	    for(l in 1:N.covariates){
	      gamma[l] ~ dnorm(0, 0.001)
	    }

      delta ~ dnorm(0, 0.001)
    }
  "
  
  # First approach
  samps.all <- list()
  k_star <- list()
  sigma_star <- list()
  samps.sd <- matrix(NA, nrow = nstudies, ncol = 1 + n.covariates)
  samps.mean <- matrix(NA, nrow = nstudies, ncol = 1 + n.covariates)
  
  for(i in 1:nstudies){
    data_ith_study <- studies_data[[i]]
    
    BayesLASSO.spec <- textConnection(Model.ITE)
    
    input_data <- list(y = data_ith_study$y,
                       x = as.matrix(select(data_ith_study, starts_with("x"))),
                       x_int = as.matrix(select(data_ith_study, starts_with("x")) *
                                           data_ith_study$treatment),
                       trt = data_ith_study$treatment, 
                       n_patients = nrow(data_ith_study),
                       N.covariates = n.covariates)
    
    jags.m <- jags.model(file = BayesLASSO.spec, data = input_data,
                         n.chains = 2, n.adapt = 1000, quiet = T)
    
    update(jags.m, n.iter = 2000, progress.bar = "none")
    
    params <- c("delta", "gamma")
    
    samps <- coda.samples(model = jags.m, params, n.iter = 1e4, 
                          progress.bar = "none")
    samps.all[[i]] <- do.call(rbind, samps)
    
    # Bayesian meta-analysis
    param_labels <- c("delta", glue("gamma{1:n.covariates}"))
    k_star[[i]] <- colMeans(samps.all[[i]]); names(k_star[[i]]) <- param_labels
    sigma_star[[i]] <- cov(samps.all[[i]])
    colnames(sigma_star[[i]]) <- rownames(sigma_star[[i]]) <- param_labels
  }
  
  closeAllConnections()
  
  # to return from first approach
  inv_sigma <- lapply(sigma_star, solve)
  
  k_numeric <- bind_rows(k_star)
  inv_sigma_array3D <- array(
    unlist(inv_sigma),
    dim = c(nrow(sigma_star[[1]]), ncol(sigma_star[[1]]), nstudies)
  )
  
  to_return.1 <- list(study_specific_effect_sizes = k_numeric,
                      inv_var_cov_matrix = inv_sigma_array3D,
                      nstudies = nstudies,
                      n.covariates = n.covariates)
  
  
  # Final return
  return(list(to_return.1, samps.all))
  
}

get_ITE <- function(covariates_test_data, parameters_posterior_sample){
  x <- cbind(1, covariates_test_data) %*% t(parameters_posterior_sample)
  to_return <- list(point_estimates = rowMeans(x),
                    weights = apply(X = x, MARGIN = 1, FUN = var))
  return(to_return)
}

second_stage_ITE <- function(covariates_test_data, posterior_samples){
  samps.mean <- matrix(ncol = nstudies, nrow = 3e3)
  samps.sd <- matrix(ncol = nstudies, nrow = 3e3)
  for(l in 1:nstudies){
    z <- get_ITE(covariates_test_data = covariates_test_data,
                 parameters_posterior_sample = posterior_samples[[l]])
    samps.mean[, l] <- z$point_estimates
    samps.sd[, l] <- z$weights
  }
  
  ITE.to.return <- rowSums(samps.mean * 1 / samps.sd) / rowSums(1 / samps.sd)
  
  return(ITE.to.return)
}



# Second stage
second_stage <- function(obj_first_stage){
  model.meta_analysis <- "
  model {
    for(i in 1:n_studies){
      y[i, ] ~ dmnorm(Mu[i, 1:total_number_covariates], inv_var_cov_matrix[, , i])
      
      Mu[i, 1] <- delta[i]
      
      for(j in 1:N.covariates){
        Mu[i, j + 1] <- gamma[j]
      }
    }
    
    # Fixed effects for gamma
    for(j in 1:N.covariates){
      # for(i in 1:n_studies){
        # gamma[i, j] ~ dnorm(theta_gamma, 0.001)
      #
      gamma[j] ~ dnorm(0, 0.001)
      # theta_gamma[j] ~ dnorm(0, 0.001)
      # prior for precisions
      # tau_gamma[j] ~ dgamma(0.001, 0.001)
    }

    # Random effects for delta
    for(i in 1:n_studies){
      delta[i] ~ dnorm(theta_delta, tau_delta)
    }
    theta_delta ~ dnorm(0, 0.001)
    tau ~ dnorm(0, 1) T(0, )
    tau_sq <- pow(tau, 2)
    tau_delta <- 1 / tau_sq
  }"
  
  
  model.MA.spec <- textConnection(model.meta_analysis)
  data_to_pass <- list(n_studies = obj_first_stage[["nstudies"]], 
                       N.covariates = obj_first_stage[["n.covariates"]],
                       y = obj_first_stage[["study_specific_effect_sizes"]], 
                       total_number_covariates = 1 + obj_first_stage[["n.covariates"]], # effect modifiers + treatment variable
                       inv_var_cov_matrix = obj_first_stage[["inv_var_cov_matrix"]])
  
  MA.jags.model <- jags.model(file = model.MA.spec, data = data_to_pass,
                              n.chains = 2, quiet = T)
  closeAllConnections()
  
  params <- c("theta_delta", "gamma")
  samps <- coda.samples(MA.jags.model, variable.names = params, n.iter = 2000,
                        progress.bar = "none")
  samps_dataset <- do.call(rbind, samps)
  to_return <- data.frame(parameter = colnames(samps_dataset),
                          Estimate = colMeans(samps_dataset))
  
  return(to_return)
}


BayesLASSO.I_stage <- function(studies_data, n.covariates,
                               nstudies){
  samps.all <- list()
  BayesLASSO <- "
   model{
  for(i in 1:n_patients){
    y[i] ~ dnorm(mu[i], prec)
    mu[i] <- alpha + inprod(beta[], x[i, ]) + inprod(gamma[], x_trt[i, ]) +
             delta * trt[i]
  }

  # Likelihood
	prec <- pow(sigma, -2)
  sigma ~ dunif(0, 10)

  # Intercept
  alpha ~ dnorm(0, 0.001)

  # Main effects (weakly informative Gaussian)
  for(j in 1:N.covariates){
    beta[j] ~ dnorm(0, 0.001)
  }

  # LASSO on interactions
  
  for(j in 1:N.covariates){
    gamma[j] ~ dnorm(0, prec / tau2[j])
    tau2[j] ~ dexp(lambda2 / 2)
  }
  lambda2 ~ dgamma(0.01, 0.01)
  
  # Treatment main effect
  delta ~ dnorm(0, 0.001)
 }
  "
  
  params <- c("delta", "gamma")

  for(i in 1:nstudies){
    data_ith_study <- studies_data %>% filter(ID_study == i)
    
    BayesLASSO.spec <- textConnection(BayesLASSO)
    cov <- data_ith_study %>% select(starts_with("x_")) %>% as.matrix()
    
    input_data <- list(y = data_ith_study$y,
                       x = cov,
                       x_trt = cov * data_ith_study$treatment,
                       trt = data_ith_study$treatment, 
                       n_patients = nrow(data_ith_study),
                       N.covariates = n.covariates)
    
    jags.m <- jags.model(file = BayesLASSO.spec, data = input_data,
                         n.chains = 2, n.adapt = 100, quiet = T)
    
    update(jags.m, n.iter = 2000, progress.bar = "none")
    
    samps <- coda.samples(model = jags.m, params, n.iter = 5e3, 
                          progress.bar = "none")
    
    samps.all[[i]] <- do.call(rbind, samps)
  }
  close(BayesLASSO.spec)
  return(samps.all)
}    


second_stage_LASSO <- function(covariates_test_data, posterior_samples,
                               s.means, s.sds){
  samps.mean <- matrix(ncol = nstudies, nrow = 3e3)
  samps.sd <- matrix(ncol = nstudies, nrow = 3e3)
  for(l in 1:nstudies){
    covariates_test_data_std <- scale(covariates_test_data,
                                      center = s.means[l, -1],
                                      scale = s.sds[l, -1])
    
    z <- get_ITE(covariates_test_data = covariates_test_data_std,
                 parameters_posterior_sample = posterior_samples[[l]])
    
    samps.mean[, l] <- z$point_estimates
    samps.sd[, l] <- z$weights
  }
  
  ITE.to.return <- rowSums(samps.mean * 1 / samps.sd) / rowSums(1 / samps.sd)
  return(ITE.to.return)
}











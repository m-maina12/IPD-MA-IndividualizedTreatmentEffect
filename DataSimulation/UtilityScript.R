# Utility script
# Function for Data Generation

# Required libraries: dplyr, glue
library(dplyr)
library(glue)
library(mvtnorm)

generate_cont_data <- function(
    nstudies,
    tau_effect,
    beta,
    gamma,
    treatment_effect,
    ss_min,
    ss_max,
    complexity,
    N.covariates = 8
){
  covariates_names <- glue::glue("x_{1:N.covariates}")
  N <- round(runif(nstudies, ss_min, ss_max))
  
  x_mean <- rep(0, N.covariates)
  x_cov  <- toeplitz(0.35^(0:(N.covariates - 1)))

  # thresholds per categorical transformation
  threshold_x2 <- 0
  threshold_x4 <- runif(1, -1, 1)
  threshold_x6 <- runif(1, 0, 1)
  threshold_x8 <- runif(1, -1, 0)
  
  # random study effects
  study_re <- rnorm(nstudies, 0, sqrt(0.5))
  HTE_study <- rnorm(nstudies, 0, sqrt(tau_effect))
  
  intercept <- 2
  
  # ---- ALLOC OUTPUT ----
  data_list <- vector("list", nstudies)
  
  # ---- TRAINING DATA ----
  for (i in 1:nstudies) {
    n_i <- N[i]
    
    # treatment + patient noise
    treatment <- rbinom(n_i, 1, 0.5)
    patient_noise <- rnorm(n_i, 0, sqrt(0.5))
    
    # PATIENT-LEVEL MEAN-SHIFT 
    mean_shift <- matrix(rnorm(n_i * N.covariates, 0, sqrt(1.5)),
                         nrow = n_i, ncol = N.covariates)
    
    # COVARIATES GENERATION MVN 
    cov <- rmvnorm(n_i, mean = x_mean, sigma = x_cov) + mean_shift
    cov <- as.data.frame(cov)
    colnames(cov) <- covariates_names
    
    # binary transformations
    cov$x_2 <- as.numeric(cov$x_2 < threshold_x2)
    cov$x_4 <- as.numeric(cov$x_4 < threshold_x4)
    cov$x_6 <- as.numeric(cov$x_6 < threshold_x6)
    cov$x_8 <- as.numeric(cov$x_8 < threshold_x8)
    
    # nonlinear parts
    if (complexity) {
      int <- with(cov, 1.5 * log(1 + abs(x_1)) + 1.5 * x_2 * sqrt(x_3^2 + 1) -
                    1.4 * (abs(x_5) + 2) -  1.8 * (abs(x_3*x_5))^(2/3))
    } else {
      int  <- 0
    }
    
    # linear parts
    lin  <- as.matrix(cov) %*% beta
    int_lin <- as.matrix(cov) %*% gamma
    
    # outcome
    y <- intercept +
      study_re[i] +
      lin +
      (treatment_effect + HTE_study[i]) * treatment +
       (int + int_lin) * treatment +
      patient_noise
    
    data_list[[i]] <- data.frame(
      ID_study = i,
      ID_patient = 1:n_i,
      cov,
      treatment,
      y
    )
  }
  
  data <- dplyr::bind_rows(data_list)
  
  # ---- TESTING DATA ----
  N_study_testing <- 3
  N_patient_testing <- 1000
  
  study_re_test <- rnorm(N_study_testing, 0, sqrt(0.5))
  HTE_test <- rnorm(N_study_testing, 0, sqrt(tau_effect))
  
  ITE_list <- vector("list", N_study_testing)
  
  for (t in 1:N_study_testing) {
    mean_shift <- matrix(rnorm(N_patient_testing * N.covariates, 0, sqrt(1.5)),
                         nrow = N_patient_testing, ncol = N.covariates)
    
    cov_test <- rmvnorm(N_patient_testing,
                        mean = x_mean, sigma = x_cov) + mean_shift
    cov_test <- as.data.frame(cov_test)
    colnames(cov_test) <- covariates_names
    
    # binary transformations
    cov_test$x_2 <- as.numeric(cov_test$x_2 < threshold_x2)
    cov_test$x_4 <- as.numeric(cov_test$x_4 < threshold_x4)
    cov_test$x_6 <- as.numeric(cov_test$x_6 < threshold_x6)
    cov_test$x_8 <- as.numeric(cov_test$x_8 < threshold_x8)
    
    # nonlinear
    if (complexity) {
      int <- with(cov_test, 1.5 * log(1 + abs(x_1)) + 1.5 * x_2 * sqrt(x_3^2 + 1) -
                    1.4 * (abs(x_5) + 2) -  1.8 * (abs(x_3*x_5))^(2/3))
    } else {
      int <- 0
    }
    
    int_lin <- as.matrix(cov_test) %*% gamma
    
    ITE_list[[t]] <- data.frame(
      ID_study = t,
      ID_patient = 1:N_patient_testing,
      cov_test,
      true_ITE_cont = treatment_effect + HTE_test[t] + int + int_lin
    )
  }
  
  true_ITE <- dplyr::bind_rows(ITE_list)
  
  list(data = data, true_ITE = true_ITE)
}

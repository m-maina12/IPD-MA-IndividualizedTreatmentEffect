# Utility script
# Function for model performances
assess_performances <- function(ITE_actual, ITE_predicted){
  ok <- complete.cases(ITE_actual, ITE_predicted)
  ITE_actual <- ITE_actual[ok]
  ITE_predicted <- ITE_predicted[ok]
  
  # Median Absolute Error
  MAE <- median(abs(ITE_actual - ITE_predicted))
  
  # Mean Bias
  MB <- mean(ITE_actual - ITE_predicted)
  
  # R-squared
  R2 <- 1 - mean((ITE_actual - ITE_predicted)^2)/var(ITE_actual)
  
  to_return <- c(MAE = MAE, 
                 MSE = MSE, 
                 MB = MB,
                 R2 = R2)
  
  return(t(to_return))
}

# Pooling methods

# 2. Bootstrap

# We said 100 bootstrap sample

get_bootstrap_sample <- function(X){
  n <- nrow(X)
  index <- sample(1:n, size = n, replace = T)
  return(X[index, ])
}

# 3. Linear model variance
get_linear_model_variance <- function(data_train,
                                      new_patients_data,
                                      formula) {
  
  # Model fit
  ausiliar_model <- lm(formula, data = data_train)
  
  # Cov matrixes
  Vbeta <- vcov(ausiliar_model)
  
  # New patients data; under treatment and control
  data_1 <- new_patients_data %>% mutate(treatment = 1)
  data_0 <- new_patients_data %>% mutate(treatment = 0)
  
  formula_no_y <- delete.response(terms(formula))
  x1_test <- model.matrix(formula_no_y, data = data_1)
  x0_test <- model.matrix(formula_no_y, data = data_0)
  
  diff_x <- x1_test - x0_test
  
  # Linear Variances
  variances <- rowSums((diff_x %*% Vbeta) * diff_x)
  
  return(as.numeric(variances))
}

# LINEAR MODEL VARIANCES

# Names of the covariates and create the formula

x_vars <- grep("^x", names(data_training[[1]][[1]][[1]]), value = TRUE)
formula_str <- paste("y ~ treatment * (", paste(x_vars, collapse = " + "), ")", 
                     sep = "")
model_formula <- as.formula(formula_str)

linear_model_variances <- vector("list", length = nscenarios)
for(i in 1:nscenarios){
  temp_matrix <- vector("list", length = nsim)
  for(j in 1:nsim){
    temp_simulation <- matrix(nrow = nrow(data_testing[[i]][[j]]),
                              ncol = nstudies)
    for(l in 1:nstudies){
      temp_simulation[, l] <- get_linear_model_variance(data_train = data_training[[i]][[j]][[l]], 
                                                        new_patients_data = data_testing[[i]][[j]],
                                                        formula = model_formula)
    }
    temp_matrix[[j]] <- temp_simulation
  }
  linear_model_variances[[i]] <- temp_matrix
}


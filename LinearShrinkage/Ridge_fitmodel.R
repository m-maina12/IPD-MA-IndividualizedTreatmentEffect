# Ridge 
library(glmnet)
library(matrixStats)

ridge_build_design_matrix <- function(X, mod.formula){
  return(model.matrix(mod.formula, data = X)[, -1])
}


ridge_predict_ite <- function(model, # Model object
                              x1_test,
                              x0_test){
  ite <- as.numeric(predict(model, newx = x1_test, s = "lambda.min") -
                      predict(model, newx = x0_test, s = "lambda.min"))

  return(ite)
}

ridge_bootstrap_variance <- function(data_train, # Original training data
                                     newX, # New covariates
                                     model_formula, # Model formula
                                     penalty_factors, # Penalty factors for the regression 
                                     lambda_seq, # Grid of lambda parameters
                                     N_boot = 100, # Numbers of bootstrap iteration
                                     x1_test,
                                     x0_test){

  bootstrap_preds <- matrix(nrow = nrow(newX), ncol = N_boot)

  for(m in 1:N_boot){
      bs_data <- get_bootstrap_sample(X = data_train)
      x_train_bs <- ridge_build_design_matrix(X = bs_data, mod.formula = model_formula)
      cv_bs <- cv.glmnet(
          x = x_train_bs,
          y = bs_data$y,
          nfolds = 5,
          alpha = 0,
          penalty.factor = penalty_factors,
          lambda = lambda_seq
      )

      bootstrap_preds[, m] <- ridge_predict_ite(model = cv_bs, x1_test = x1_test, x0_test = x0_test)
  }

  matrixStats::rowVars(bootstrap_preds)
}


get_variance_wrapper <- function(second_stage, 
                                 newX,
                                 model_formula,
                                 penalty_factors = NULL,
                                 lambda_seq = NULL,
                                 N_boot = 100,
                                 x1_test,
                                 x0_test){

      second_stage <- toupper(second_stage)
      
      if(second_stage == "OLS+IV"){
        # ritorna una funzione che prende solo data_train
        return(function(data_train){
          get_linear_model_variance(
            data_train = data_train,
            newX = newX,
            model_formula = model_formula
          )
        })
        
      } else if(second_stage == "BS"){
        # ritorna una funzione che prende solo data_train
        return(function(data_train){
          ridge_bootstrap_variance(
            data_train = data_train,
            newX = newX,
            model_formula = model_formula,
            penalty_factors = penalty_factors,
            lambda_seq = lambda_seq,
            N_boot = N_boot,
            x1_test = x1_test,
            x0_test = x0_test
          )
        })
        
      } else {
        stop("second_stage must be either 'OLS+IV' or 'BS'")
      }
    } 


ridge_fit_model <- function(data, model_formula, lambda_seq, 
                            penalty_factors){

      x_train <- ridge_build_design_matrix(X = data, mod.formula = model_formula)
      
      # fit cv.glmnet su dataset originale
      cv.obj <- cv.glmnet(x_train, y = data$y, nfolds = 10, alpha = 0,
                            penalty.factor = penalty_factors, lambda = lambda_seq)

      return(cv.obj)
    
}


# 1. Model specification

ridge_ipd_spec <- function(covariate_names,
                           outcome = "y",
                           treatment = "treatment",
                           lambda_seq = NULL) {

  if(is.null(lambda_seq)){
    lambda_seq <- 10^seq(2, -3, by = -0.3)
  }

  formula_str <- paste(
  "y ~ treatment * (",
  paste(covariate_names, collapse = " + "),
  ")"
  )

  model_formula <- as.formula(formula_str)

  penalty_factors <- c(rep(0, 1 + length(covariate_names)), rep(1, length(covariate_names)))

  list(
    covariate_names = covariate_names,
    outcome = outcome,
    treatment = treatment,
    lambda_seq = lambda_seq,
    model_formula = model_formula,
    penalty_factors = penalty_factors
  )
}

# 2. Training
ridge_ipd_train <- function(data, spec){

  ridge.models <- vector("list", length = length(data))
  nstudies <- length(data)

  for(l in 1:nstudies){
     # fit cv.glmnet su dataset originale
     data[[l]]$y <- data[[l]][[spec$outcome]]
     data[[l]]$treatment <- data[[l]][[spec$treatment]]
     ridge.models[[l]] <- ridge_fit_model(data = data[[l]], model_formula = spec$model_formula, 
                                lambda_seq = spec$lambda_seq, penalty_factors = spec$penalty_factors)       
    }

  to_return <- list(
    data = data,
    spec = spec,
    nstudies = nstudies,
    models = ridge.models
  )

  class(to_return) <- "ipd_ridge"

  return(to_return)
}



# 3. Prediction
predict.ipd_ridge <- function(object,
                              newX,
                              second_stage = "BS",
                              N_boot = 100){

  predictions <- matrix(nrow = nrow(newX), ncol = object$nstudies)
  variance <- matrix(nrow = nrow(newX), ncol = object$nstudies)
  x1_test <- model.matrix(object$spec$model_formula, data = cbind(newX, treatment = 1))[,-1]
  x0_test <- model.matrix(object$spec$model_formula, data = cbind(newX, treatment = 0))[,-1] 

  variance_method <- get_variance_wrapper(
      second_stage = second_stage,
      newX = newX,
      model_formula = object$spec$model_formula,
      penalty_factors = object$spec$penalty_factors,
      lambda_seq = object$spec$lambda_seq,
      N_boot = N_boot,
      x1_test = x1_test,
      x0_test = x0_test
    )
  
  for(l in seq_len(object$nstudies)){
    predictions[, l] <- ridge_predict_ite(model = object$models[[l]], x1_test = x1_test, x0_test = x0_test)
    variance[, l] <- variance_method(object$data[[l]])
  }

  weights <- 1 / variance
  pooled_predictions <- rowSums(predictions * weights) / rowSums(weights)

  return(pooled_predictions)
}
















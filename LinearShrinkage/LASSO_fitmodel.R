# LASSO 
library(glmnet)
library(matrixStats)

lasso_build_design_matrix <- function(X, mod.formula){
  return(model.matrix(mod.formula, data = X)[, -1])
}


lasso_predict_ite <- function(model, # Model object
                              x1_test,
                              x0_test){
  ite <- as.numeric(predict(model, newx = x1_test, s = "lambda.min") -
                      predict(model, newx = x0_test, s = "lambda.min"))

  return(ite)
}

lasso_bootstrap_variance <- function(data_train, # Original training data
                                     newX, # New covariates
                                     model_formula, # Model formula
                                     penalty_factors, # Penalty factors for the regression 
                                     lambda_seq, # Grid of lambda parameters
                                     N_boot = 100 # Numbers of bootstrap iteration
                                     ){

  bootstrap_preds <- matrix(nrow = nrow(newX), ncol = N_boot)

  for(m in 1:N_boot){
      bs_data <- get_bootstrap_sample(X = data_train)
      x_train_bs <- lasso_build_design_matrix(X = bs_data, mod.formula = model_formula)
      cv_bs <- cv.glmnet(
          x = x_train_bs,
          y = bs_data$y,
          nfolds = 5,
          alpha = 1,
          penalty.factor = penalty_factors,
          lambda = lambda_seq
      )

      bootstrap_preds[, m] <- lasso_predict_ite(model = cv_bs, newX = newX)
  }

  matrixStats::rowVars(bootstrap_preds)
}


get_variance_wrapper <- function(second_stage, 
                                 newX,
                                 model_formula,
                                 penalty_factors = NULL,
                                 lambda_seq = NULL,
                                 N_boot = 100){

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
          lasso_bootstrap_variance(
            data_train = data_train,
            newX = newX,
            model_formula = model_formula,
            penalty_factors = penalty_factors,
            lambda_seq = lambda_seq,
            N_boot = N_boot
          )
        })
        
      } else {
        stop("second_stage must be either 'OLS+IV' or 'BS'")
      }
    } 


LASSO_fit_model <- function(data, model_formula, lambda_seq, 
                            penalty_factors){

      x_train <- lasso_build_design_matrix(X = data, mod.formula = model_formula)
      
      # fit cv.glmnet su dataset originale
      cv.obj <- cv.glmnet(x_train, y = data$y, nfolds = 10, alpha = 1,
                            penalty.factor = penalty_factors, lambda = lambda_seq)

      return(cv.obj)
    
}


lasso_ipd_ite <- function(data, 
                          newX, 
                          covariate_names, 
                          outcome = "response",
                          treatment = "trt",  
                          second_stage, 
                          N_sample_bootstrap = 100
                          ){
    n_covariates <- length(covariate_names)
    nstudies <- length(data)

    lambda_seq <- 10 ^ seq(2, -3, by = -0.3)
    penalty_factors <- c(rep(0, 1 + n_covariates), rep(1, n_covariates))

    formula_str <- paste("y ~ ", treatment, " * (", paste(covariate_names, collapse = " + "), ")", 
                         sep = "")
    model_formula <- as.formula(formula_str)

    predictions <- variance <- matrix(nrow = nrow(newX), ncol = nstudies)
    x1_test <- model.matrix(model_formula, data = cbind(newX, treatment = 1))[,-1]
    x0_test <- model.matrix(model_formula, data = cbind(newX, treatment = 0))[,-1] 
    
    variance_method <- get_variance_wrapper(
      second_stage = second_stage,
      newX = newX,
      model_formula = model_formula,
      penalty_factors = penalty_factors,
      lambda_seq = lambda_seq,
      N_boot = N_sample_bootstrap
    )

    for(l in 1:nstudies){
     # fit cv.glmnet su dataset originale
     data[[l]]$y <- data[[l]][[outcome]]
      cv.obj <- LASSO_fit_model(data = data[[l]], model_formula = model_formula, 
                                lambda_seq = lambda_seq, penalty_factors = penalty_factors,
                                cov_names = covariate_names)
      
      # predictions per ITE
      predictions[, l] <- lasso_predict_ite(model = cv.obj, x1_test = x1_test, x0_test = x0_test) 

      variance[, l] <- variance_method(data[[l]])

      
    }

    weights <- 1 / variance
    pooled_predictions <- rowSums(predictions * weights) / rowSums(weights)

    return(pooled_predictions)
}






















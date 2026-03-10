# LASSO 
library(glmnet)
library(matrixStats)

LASSO_fit <- function(y, X, t, newX, second_stage, N_sample_bootstrap = 100){
    # y = list of outcomes
    # X = list of datasets with covariates
    # t = list of treatment

    ols_iv_variance <- function(data_train, newX, model_formula){
      get_linear_model_variance(
          data_train = data_train,
          new_patients_data = newX,
          formula = model_formula
      )
    }

    bootstrap_variance <- function(data_train, newX, model_formula,
                               x1_test, x0_test,
                               penalty_factors, lambda_seq,
                               N_boot = 100){

  bootstrap_preds <- matrix(nrow = nrow(newX), ncol = N_boot)

  for(m in 1:N_boot){

      bs_data <- get_bootstrap_sample(X = data_train)

      x_train_bs <- model.matrix(model_formula, data = bs_data)[, -1]

      cv_bs <- cv.glmnet(
          x_train_bs,
          y = bs_data$y,
          nfolds = 5,
          alpha = 1,
          penalty.factor = penalty_factors,
          lambda = lambda_seq
      )

      bootstrap_preds[, m] <- as.numeric(
          predict(cv_bs, newx = x1_test, s = "lambda.min") -
          predict(cv_bs, newx = x0_test, s = "lambda.min")
      )
  }

  matrixStats::rowVars(bootstrap_preds)
}

    n_covariates <- ncol(X[[1]])
    nstudies <- length(y)

    lambda_seq <- 10 ^ seq(2, -3, by = -0.3)
    penalty_factors <- c(rep(0, 1 + n_covariates), rep(1, n_covariates))

    new_cov_names <- paste0("x", 1:n.covariates)
    formula_str <- paste("y ~ treatment * (", paste(new_cov_names, collapse = " + "), ")", 
                         sep = "")
    model_formula <- as.formula(formula_str)

    x1_test <- cbind(1, newX, newX) 
    x0_test <- cbind(0, newX, 0 * newX) 
    predictions <- variance <- matrix(nrow = nrow(newX), ncol = nstudies)
    
    second_stage <- toupper(second_stage)
    if(second_stage == "OLS+IV"){
      variance_method <- ols_iv_variance
    } else (second_stage == "BS"){
      variance_method <- bootstrap_variance
    } 

    for(l in 1:nstudies){
      data_lth_study <- cbind(X[[l]], t[[l]], y[[l]]); colnames(data_lth_study) <- c(new_cov_names, "treatment", "y")
      x_train <- model.matrix(model_formula, data = data_lth_study)[, -1]
      
      # fit cv.glmnet su dataset originale
      cv.obj <- cv.glmnet(x_train, y = y[[l]], nfolds = 10, alpha = 1,
                            penalty.factor = penalty_factors, lambda = lambda_seq)
      
      # predictions per ITE
      predictions[, l] <- as.numeric(predict(cv.obj, newx = x1_test, s = "lambda.min") -
                               predict(cv_lasso, newx = x0_test, s = "lambda.min"))

      variance[, l] <- variance_method(
      data_train = data_lth_study,
      newX = newX,
      model_formula = model_formula,
      x1_test = x1_test,
      x0_test = x0_test,
      penalty_factors = penalty_factors,
      lambda_seq = lambda_seq
  )

      
    }

    weights <- 1 / variance
    pooled_predictions <- rowSums(predictions * weights) / rowSums(weights)

    return(pooled_predictions)
}

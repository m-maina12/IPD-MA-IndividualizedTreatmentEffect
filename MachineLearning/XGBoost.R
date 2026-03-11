# XGBoost
library(dplyr)
library(matrixStats)
library(xgboost)


params_df <- expand.grid(
  max_depth = 2:6,
  eta = c(0.1, 0.01, 0.001),  
  nrounds = c(100, 500, 1000)
)

xgb_predict_ite <- function(model, x1_test, x0_test){
	xgb_test_1 <- xgb.DMatrix(data = x1_test)
    xgb_test_0 <- xgb.DMatrix(data = x0_test)
	return(predict(model, newdata = xgb_test_1) -
	          predict(model, newdata = xgb_test_0))
}

xgb_bootstrap_variance <- function(data_train, # Original training data
                                     newX, # New covariates
                                     model_formula, # Model formula
                                     penalty_factors, # Penalty factors for the regression 
                                     lambda_seq, # Grid of lambda parameters
                                     N_boot = 100, # Numbers of bootstrap iteration
                                     x1_test,
                                     x0_test){

  bootstrap_preds <- matrix(nrow = nrow(x1_test), ncol = N_boot)

  for(m in 1:N_boot){
        bs_data <- get_bootstrap_sample(X = data_train)
        x_bs <- bs_data %>% select(treatment, starts_with("x_")) %>% as.matrix()
        y_bs <- bs_data$y
        xgboost.data <- xgb.DMatrix(data = x_bs, label = y_bs)
        xgb_bs <- xgboost(data = xgboost.data, params = best_params,
                          nrounds = best_row$nrounds, verbose = 0)
        bootstrap_preds[, m] <- predict(xgb_bs, newdata = xgb_test_1) -
          predict(xgb_bs, newdata = xgb_test_0)
      }

  matrixStats::rowVars(bootstrap_preds)
}

xgb_get_variance_wrapper <- function(second_stage, 
	                                 newX,
	                                 N_boot = 100,
	                                 x1_test,
	                                 x0_test, 
	                                 hyperparameters){

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
            hyperparameters = hyperparameters,
            N_boot = N_boot,
            x1_test = x1_test,
            x0_test = x0_test
          )
        })
        
      } else {
        stop("second_stage must be either 'OLS+IV' or 'BS'")
      }
    } 


xgb_cross_validation <- function(xgb.data, parameters){
	performances <- numeric(length = nrow(parameters))

	# devo aggiugnere un check sui parametri
	for(i in 1:nrow(parameters)){
		params <- list(objective = "reg:squarederror", 
				   eval_metric = "rmse",
                   max_depth = parameters$max_depth[i], 
                   eta = parameters$eta[i],
                   colsample_bytree = 1, 
                   min_child_weight = 1, 
                   subsample = 1)

		cv <- xgb.cv(params = params, data = xgb.data, 
                     nrounds = parameters$nrounds[i], nfold = 10,
                     early_stopping_rounds = 20, verbose = 0
        )
        performances[[i]] <- min(cv$evaluation_log$test_rmse_mean)
	}

	best_row <- params_df[which.min(cv_results), ]

	to_return <- list(best_nrounds = best_row$nrounds,
		best_params = list(objective = "reg:squarederror",
                        eval_metric = "rmse", max_depth = best_row$max_depth,
                        eta = best_row$eta, colsample_bytree = 1,
                        min_child_weight = 1, subsample = 1))
      
    return(to_return)
}

xgb_fit_model <- function(){}

# Training
xgb_train <- function(data, params_df, covariate_names, 
                        outcome = "y",
                        treatment = "treatment"){
	nstudies <- length(data)
	ncovariates <- length(covariate_names)

	xgb.models <- vector("list", length = nstudies)
	hyperparameters <- vector("list", length = nstudies)

	for(l in 1:nstudies){
		data_lth_study <- data[[l]]
		x_train <- as.matrix(data[[l]][, covariate_names])
		y_train <- data[[l]][[outcome]]
		dtrain <- xgb.DMatrix(data  = x_train, label = y_train)

		# Cross-Validation for hyperparameters tuning 
		hyperparameters[[l]] <- xgb_cross_validation(xgb.data = dtrain, parameters = params_df)

		# Training model
      	xgb.models[[l]] <- xgb.train(data = dtrain, params = hyperparameters[[l]]$best_params,
                               nrounds = hyperparameters[[l]]$best_nrounds)
	}

	return(list(final_models = xgb.models,
		final_hyperparameters = hyperparameters))
}

# Prediction
predict.xgb <- function(newX, obj_xgb_train, second_stage, N_boot){
	x1_test <- as.matrix(cbind(treatment = 1, newX))
    x0_test <- as.matrix(cbind(treatment = 0, newX))  

    predictions <- matrix(nrow = nrow(newX), ncol = object$nstudies)
  	variance <- matrix(nrow = nrow(newX), ncol = object$nstudies) 


    variance_method <- xgb_get_variance_wrapper(
      second_stage = second_stage,
      newX = newX,
      model_formula = ,
      hyperparameters = obj_xgb_train$final_hyperparameters,
      N_boot = N_boot,
      x1_test = x1_test,
      x0_test = x0_test
    )
  
  for(l in seq_len(object$nstudies)){
    predictions[, l] <- xgb_predict_ite(model = obj_xgb_train$final_models[[l]], x1_test = x1_test, x0_test = x0_test)
    variance[, l] <- variance_method(object$data[[l]])
  }

  weights <- 1 / variance
  pooled_predictions <- rowSums(predictions * weights) / rowSums(weights)

  return(pooled_predictions)


}      
      
      
      
      pred_xgb <- predict(final_model, dtest1) - predict(final_model, dtest0)
      
      bootstrap_preds <- matrix(NA, nrow = nrow(X_test), ncol = N_sample_bootstrap)
      
      
      for(m in 1:N_sample_bootstrap){
        bs_data <- get_bootstrap_sample(X = data_l)
        x_bs <- bs_data %>% select(treatment, starts_with("x_")) %>% as.matrix()
        y_bs <- bs_data$y
        xgboost.data <- xgb.DMatrix(data = x_bs, label = y_bs)
        xgb_bs <- xgboost(data = xgboost.data, params = best_params,
                          nrounds = best_row$nrounds, verbose = 0)
        bootstrap_preds[, m] <- predict(xgb_bs, newdata = xgb_test_1) -
          predict(xgb_bs, newdata = xgb_test_0)
      }
      
      list(pred_xgb = pred_xgb,
           bootstrap_var = rowVars(bootstrap_preds))
    }
  
  




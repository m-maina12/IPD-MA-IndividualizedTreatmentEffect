# LASSO 
library(glmnet)
library(matrixStats)

prepare_lasso_training_data <- function(X, formula){
  x_train <- model.matrix(formula, data = X)
  return(x_train[, -1])
}

prepare_lasso_testing_data <- function(X, formula){
  X1 <- X %>% mutate(treatment = 1)
  X0 <- X %>% mutate(treatment = 0)
  formula_no_y <- delete.response(terms(formula))
  
  x1_test <- model.matrix(formula_no_y, data = X1)[, -1]
  x0_test <- model.matrix(formula_no_y, data = X0)[, -1]
  to_return <- list(x1_test = x1_test, x0_test = x0_test)
}


lambda_seq <- 10 ^ seq(2, -3, by = -0.3)
penalty_factors <- c(rep(0, 1 + n_covariates), rep(1, n_covariates))


obj_test <- prepare_lasso_testing_data(X = data_testing, formula = model_formula)
x1_test <- obj_test$x1_test; x0_test <- obj_test$x0_test

lasso_predictions <- bootstrap_variances <- matrix(nrow = nrow(x1_test), ncol = nstudies)

for(l %in% 1:nstudies){
  data_l <- data_training[[l]]
  x_train <- prepare_lasso_training_data(X = data_l, formula = model_formula)
  
  # fit cv.glmnet su dataset originale
  cv_lasso <- cv.glmnet(x_train, y = data_l$y, nfolds = 10, alpha = 1,
                        penalty.factor = penalty_factors, lambda = lambda_seq)
  
  # predictions per ITE
  lasso_predictions[, l] <- as.numeric(predict(cv_lasso, newx = x1_test, s = "lambda.min") -
                           predict(cv_lasso, newx = x0_test, s = "lambda.min"))
  
  # bootstrap
  bootstrap_preds <- matrix(nrow = nrow(x1_test), ncol = N_sample_bootstrap)
  for(m in 1:N_sample_bootstrap){
    bs_data <- get_bootstrap_sample(X = data_l)
    x_train_bs <- prepare_lasso_training_data(X = bs_data, formula = model_formula)
    cv_bs <- cv.glmnet(x_train_bs, y = bs_data$y, nfolds = 5, alpha = 1,
                       penalty.factor = penalty_factors, lambda = lambda_seq)
    bootstrap_preds[, m] <- as.numeric(predict(cv_bs, newx = x1_test, s = "lambda.min") -
                                         predict(cv_bs, newx = x0_test, s = "lambda.min"))
  }
  bootstrap_variances[, l] <- rowVars(bootstrap_preds)
}

# ponderazione bootstrap
weights_bs <- 1/bootstrap_variances
prediction_bootstrap <- rowSums(lasso_predictions * weights_bs) / rowSums(weights_bs)
LASSO_bootstrap <- cbind(j, assess_performances(ITE_actual = data_testing$true_ITE_cont,
                                                     ITE_predicted = prediction_bootstrap))

# ponderazione linear model variance
inv_var <- 1 / LM_IV
prediction_lm_variance <- rowSums(lasso_predictions * inv_var) / rowSums(inv_var)
LASSO_lm_variance <- cbind(j, assess_performances(ITE_actual = data_testing$true_ITE_cont,
                                                       ITE_predicted = prediction_lm_variance))





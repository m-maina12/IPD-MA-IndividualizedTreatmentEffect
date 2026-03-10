# Second stage: Bootstrap, OLS-IV, Bayesian MA, Weighted mixing of posteriors

bootstrapping <- function(M = 100){

}


OLS.IV <- function(predictions){
	LM_IV <- 
	inv_var <- 1 / LM_IV
	prediction_lm_variance <- rowSums(lasso_predictions * inv_var) / rowSums(inv_var)
	LASSO_lm_variance <- cbind(j, assess_performances(ITE_actual = data_testing$true_ITE_cont,
	                                                       ITE_predicted = prediction_lm_variance))

}
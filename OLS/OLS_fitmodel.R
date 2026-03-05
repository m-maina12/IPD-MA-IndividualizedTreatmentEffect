# Ordinary Least Square Regression - Inverse Variance Multivariate meta-analysis
n.parameters_to_meta_analyze <- n_covariates + 1

# Covariates of testing dataset
test_dataset <- data_testing %>% select(starts_with("x_")) %>%  as.matrix()
true_test_ITE <- data_testing[$true_ITE_cont

coefficients_to_meta_analyze <- matrix(nrow = nstudies, 
                               ncol = n.parameters_to_meta_analyze)

# Let's create a list that will contain the variance/covariance matrix of the coefficients
vcov_parameters_to_meta_analyze <- vector("list", length = nstudies)

for(l in 1:nstudies){
# Fit linear model
mod <- lm(model_formula,  data = data_training[[l]])

# Save the punctual estimates
coefficients_to_meta_analyze[l, ] <- mod$coefficients[startsWith(names(mod$coefficients), 
                                                           "treatment")]

# Save the variance/cov matrix
vcov_parameters_to_meta_analyze[[l]] <- vcov(mod)[startsWith(rownames(vcov(mod)), "treatment"), 
                                            startsWith(colnames(vcov(mod)), "treatment")]


# Meta-analysis of the parameters: 
meta_obj <- mvmeta::mvmeta(coefficients_to_meta_analyze,
                   vcov_parameters_to_meta_analyze, 
                   method = "reml") 
parameters_meta_analyzed <- coef(meta_obj)

# Predict ITE on new data
predictions <- cbind(1, test_dataset) %*% parameters_meta_analyzed

# Assess performances
performance <- cbind(j, assess_performances(ITE_actual = true_test_ITE,
                                                ITE_predicted = predictions))

}
 


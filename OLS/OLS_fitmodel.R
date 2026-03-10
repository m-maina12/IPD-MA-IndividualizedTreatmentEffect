# Ordinary Least Square Regression - Inverse Variance Multivariate meta-analysis

OLS.fit <- function(y, X, t, newX, second_stage){ # second stage can specify if the MVMA is at random or fixed effects so it can be either reml or fixed
    # y = list of outcomes
    # X = list of datasets with covariates
    # t = list of treatment
    n.covariates <- ncol(X[[1]])
    nstudies <- length(y)

    coefficients_to_meta_analyze <- matrix(nrow = nstudies, ncol = n.covariates + 1)

    # Let's create a list that will contain the variance/covariance matrix of the coefficients
    vcov_parameters_to_meta_analyze <- vector("list", length = nstudies)

    # Model formula
    new_cov_names <- paste("x", 1:n.covariates, sep = "")
    formula_str <- paste("y ~ treatment * (", paste(new_cov_names, collapse = " + "), ")", 
                         sep = "")
    model_formula <- as.formula(formula_str)

    for(l in 1:nstudies){
        # Fit linear model
        data_lth_study <- cbind(X[[l]], t[[l]]); colnames(data_lth_study) <- c(new_cov_names, "treatment")
        mod <- lm(model_formula,  data = data_lth_study)

        # Save the punctual estimates
        coefficients_to_meta_analyze[l, ] <- mod$coefficients[startsWith(names(mod$coefficients), 
                                                                   "treatment")]

        # Save the variance/cov matrix
        vcov_parameters_to_meta_analyze[[l]] <- vcov(mod)[startsWith(rownames(vcov(mod)), "treatment"), 
                                                    startsWith(colnames(vcov(mod)), "treatment")]
    }

        # Meta-analysis of the parameters: 
        meta_obj <- mvmeta::mvmeta(coefficients_to_meta_analyze,
                           vcov_parameters_to_meta_analyze, 
                           method = second_stage) 

        parameters_meta_analyzed <- coef(meta_obj)

        # Predict ITE on new data
        predictions <- cbind(1, newX) %*% parameters_meta_analyzed

  return(predictions)  
}




 


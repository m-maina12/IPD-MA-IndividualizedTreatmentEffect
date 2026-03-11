# Ordinary Least Square Regression - Inverse Variance Multivariate meta-analysis

OLS_ipd_train <- function(data, second_stage, # can be either reml or fixed
                    covariate_names, outcome, treatment){ 
    n.covariates <- ncol(X[[1]])
    nstudies <- length(y)

    coefficients_to_meta_analyze <- matrix(nrow = nstudies, ncol = n.covariates + 1)

    # Let's create a list that will contain the variance/covariance matrix of the coefficients
    vcov_parameters_to_meta_analyze <- vector("list", length = nstudies)

    # Model formula
    formula_str <- paste("y ~ treatment * (", paste(covariate_names, collapse = " + "), ")", 
                         sep = "")
    model_formula <- as.formula(formula_str)

    for(l in 1:nstudies){
        # Fit linear model
        data[[l]]$y <- data[[l]][[outcome]]
        data[[l]]$treatment <- data[[l]][[treatment]]
        mod <- lm(model_formula,  data = data[[l]])

        # Save the punctual estimates
        coef_idx <- grep(paste0("^", treatment), names(mod$coefficients))
        coefficients_to_meta_analyze[l, ] <- mod$coefficients[coef_idx]

        # Save the variance/cov matrix
        vcov_parameters_to_meta_analyze[[l]] <- vcov(mod)[coef_idx, coef_idx]
    }

        # Meta-analysis of the parameters: 
        meta_obj <- mvmeta::mvmeta(coefficients_to_meta_analyze,
                           vcov_parameters_to_meta_analyze, 
                           method = second_stage) 

        parameters_meta_analyzed <- coef(meta_obj)

        class(parameters_meta_analyzed) <- "OLS_MVMA"

  return(parameters_meta_analyzed)  
}


predict.OLS_MVMA <- function(obj.OLS, newX){
    cbind(1, newX) %*% obj.OLS
}

 


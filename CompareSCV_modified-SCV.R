
## Add range of alphas

nested_cv_diff <- function(X_m , X, Y_m, Y, Treat, range.tau, funcs_m, funcs_s, reps , n_folds,  alpha_values = c(0.01, 0.05, 0.1, 0.25, 0.5), bias_reps = NA,
                           funcs_params, funcs_params_m, n_cores = 1, verbose = F) {
  min_diff <- Inf
  best_alpha <- NULL
  result <- NULL
  all_results <- list()
  for (alpha0 in alpha_values) {

    res <- nested_cv(X, Y, funcs_s, reps, n_folds,  alpha = alpha0, bias_reps = NA,
                     funcs_params, n_cores = 1, verbose = F)
    res_m <- nested_cv_m(X_m, Y_m, Treat, range.tau, 
                         funcs_m, reps , n_folds,  alpha = alpha0, bias_reps = NA,
                         funcs_params_m, n_cores = 1, verbose = F)
    all_results[[as.character(alpha0)]] <- list(res, res_m)
    
    diff <- abs(res$ci_hi - res_m$ci_lo)
    if (diff < min_diff) {
      min_diff <- diff
      best_alpha <- alpha0
    }
  }
  
  result$best_alpha <- best_alpha
  result$all_results <- all_results
  
  return(result)
}
## linear
nested_cv_diff(X_m = data.frame(DATA_m), X = data.frame(DATA), Y_m = as.vector(Y_m), Y = as.vector(Y), 
               Treat = as.vector(Treat_m), range.tau = tau.range, funcs_m = linear_regression_funs_m, funcs_s = linear_regression_funs, reps = nested_cv_reps , n_folds = n_folds,  
               alpha_values = c(0.01, 0.05, 0.1, 0.25, 0.5), bias_reps = NA,
               funcs_params_m = NULL, funcs_params = NULL, n_cores = 1, verbose = F)
## glmboost

nested_cv_diff(X_m = data.frame(DATA_m), X = data.frame(DATA), Y_m = as.vector(Y_m), Y = as.vector(Y), 
               Treat = as.vector(Treat_m), range.tau = tau.range, funcs_m = glmboost_funs_m, funcs_s = glmboost_funs, reps = nested_cv_reps , n_folds = n_folds,  
               alpha_values = c(0.01, 0.05, 0.1, 0.25, 0.5), bias_reps = NA,
               funcs_params_m = NULL, funcs_params = NULL, n_cores = 1, verbose = F)
## glmnet
nested_cv_diff(X_m = data.frame(DATA_m), X = data.frame(DATA), Y_m = as.vector(Y_m), Y = as.vector(Y), 
               Treat = as.vector(Treat_m), range.tau = tau.range, funcs_m = gaussian_lasso_funs_m, funcs_s = gaussian_lasso_funs, reps = nested_cv_reps , n_folds = n_folds,  
               alpha_values = c(0.01, 0.05, 0.1, 0.25, 0.5), bias_reps = NA,
               funcs_params_m = list("lambdas" = lambdas_m, "best_lam" = best_lam_m), funcs_params = list("lambdas" = lambdas, "best_lam" = best_lam), n_cores = 1, verbose = F)


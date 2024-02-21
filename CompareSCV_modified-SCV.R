
## Add range of alphas

## linear
nested_cv_diff <- function(X_m , X, Y_m, Y, Treat, range.tau, funcs_m, funcs_s, reps , n_folds,  alpha_values = c(0.01, 0.05, 0.1, 0.25, 0.5), bias_reps = NA,
                           funcs_params = NULL, n_cores = 1, verbose = F) {
  min_diff <- Inf
  best_alpha <- NULL
  result <- NULL
  all_results <- list()
  for (alpha0 in alpha_values) {
    res <- nested_cv(X = data.frame(train.set), Y = as.vector(Y.train), funcs = funcs_s, reps, n_folds,  alpha = alpha0, bias_reps = NA,
                     funcs_params = NULL, n_cores = 1, verbose = F)
    res_m <- nested_cv_m(X_m = data.frame(train.set_m), Y_m = as.vector(Y.train_m), Treat = as.vector(Treat.train_m), range.tau = tau.range, 
                         funcs = funcs_m, reps = nested_cv_reps , n_folds = n_folds,  alpha = alpha0, bias_reps = NA,
                         funcs_params = NULL, n_cores = 1, verbose = F)
    all_results[[as.character(alpha0)]] <- list(res, res_m)
    
    diff_res <- res$ci_hi - res$ci_lo
    diff_res_m <- res_m$ci_hi - res_m$ci_lo
    diff <- abs(diff_res - diff_res_m)
    
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
nested_cv_diff(X_m = data.frame(train.set_m), X = data.frame(train.set), Y_m = as.vector(Y.train_m), Y = as.vector(Y.train), 
               Treat = as.vector(Treat.train_m), range.tau = tau.range, funcs_m = linear_regression_funs_m, funcs_s = linear_regression_funs, reps = nested_cv_reps , n_folds = n_folds,  
               alpha_values = c(0.01, 0.05, 0.1, 0.25, 0.5), bias_reps = NA,
                           funcs_params = NULL, n_cores = 1, verbose = F)

## glmboost

nested_cv_diff(X_m = data.frame(train.set_m), X = data.frame(train.set), Y_m = as.vector(Y.train_m), Y = as.vector(Y.train), 
               Treat = as.vector(Treat.train_m), range.tau = tau.range, funcs_m = glmboost_funs_m, funcs_s = glmboost_funs, reps = nested_cv_reps , n_folds = n_folds,  
               alpha_values = c(0.01, 0.05, 0.1, 0.25, 0.5), bias_reps = NA,
               funcs_params = NULL, n_cores = 1, verbose = F)



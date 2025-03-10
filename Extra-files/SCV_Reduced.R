

naive_cv_m <- function(X_m, Y_m, Treat, range.tau, funcs, n_folds, alpha = c(0.01, 0.05, 0.1, 0.25, 0.5),
                     trans = list(identity), funcs_params = NULL, fold_id = NULL) {
  min_gp_errors <- NULL
  if(is.null(fold_id)) {
    fold_id <- (1:nrow(X_m)) %% n_folds + 1
    fold_id <- sample(fold_id)
  }
  
  errors <- c()
  gp_errors <- c()
  for(tau in range.tau){
    for(k in 1:n_folds) {
      fit <- funcs$fitter(X_m[fold_id !=k, ], Y_m[fold_id != k], Treat[fold_id != k], tau, funcs_params = funcs_params)
      Y_m_hat <- funcs$predictor(fit, X_m[fold_id == k, ], funcs_params = funcs_params)
      error_k <- funcs$loss(Y_m_hat, Y_m[fold_id == k], Treat[fold_id == k], tau, funcs_params = funcs_params)
      errors <- c(errors, error_k)
      
      temp_vec <- c()
      for(tran in trans) {
        temp_vec <- c(temp_vec, tran(mean(error_k)))
      }
      gp_errors <- rbind(gp_errors, temp_vec)
    }
    if (is.null(min_gp_errors) || sum(gp_errors) < sum(min_gp_errors)) {
      min_gp_errors <- gp_errors
      best_tau <- tau
      best_errors <- errors
    }
  }
  
  results <- list()
  for (a in alpha) {
    ci_lo <- mean(errors) - qnorm(1 - a / 2) * sd(errors) / sqrt(length(Y))
    ci_hi <- mean(errors) + qnorm(1 - a / 2) * sd(errors) / sqrt(length(Y))
    
    results[[paste0("ci_lo_alpha_", a)]] <- ci_lo
    results[[paste0("ci_hi_alpha_", a)]] <- ci_hi
  }
  
  results[["err_hat"]] <- mean(errors)
  results[["raw_mean"]] <- mean(errors)
  results[["sd"]] <- sd(errors)
  results[["group_err_hat"]] <- apply(gp_errors, 2, mean)
  results[["group_sd"]] <- apply(gp_errors, 2, sd)
  results[["raw_errors"]] <- errors
  results[["fold_id"]] <- fold_id
  return(results)
}


# nested CV

nested_cv_m <- function(X_m, Y_m, Treat, range.tau, funcs, reps, n_folds,  alpha = c(0.01, 0.05, 0.1, 0.25, 0.5), bias_reps = NA,
                      funcs_params = NULL, n_cores = 1, verbose = F) {
  #estimate time required
  if(verbose) {
    t1 <- Sys.time()
    temp <- nested_cv_helper_m(X_m, Y_m, Treat, range.tau, funcs, n_folds,
                             funcs_params = funcs_params)
    t2 <- Sys.time()
    print(paste0("Estimated time required: ", (t2 - t1) * reps))
  }
  
  #compute out-of-fold errors on SE scale
  var_pivots <- c()
  gp_errs <- c()
  ho_errs <- c()
  tau <- c()
  if(n_cores == 1){
    raw <- lapply(1:reps, function(i){nested_cv_helper_m(X_m, Y_m, Treat, range.tau,funcs, n_folds,
                                                       funcs_params = funcs_params)})
  } else {
    raw <- parallel::mclapply(1:reps, function(i){nested_cv_helper_m(X_m, Y_m, Treat, range.tau, funcs, n_folds,
                                                                   funcs_params = funcs_params)},
                              mc.cores = n_cores)
  }
  for(i in 1:reps) {
    temp <- raw[[i]]
    var_pivots <- rbind(var_pivots, temp$pivots)
    ho_errs <- c(ho_errs, temp$errs)
    tau <- c(tau, temp$tau)
  }
  
  ########
  n_sub <- floor(length(Y_m) * (n_folds - 1) / n_folds)
  # look at the estimate of inflation after each repetition
  ugp_infl <- sqrt(max(0, mean(var_pivots[, 1]^2 - var_pivots[, 2]))) / (sd(ho_errs) / sqrt(n_sub))
  ugp_infl <- max(1, min(ugp_infl, sqrt(n_folds)))
  
  #estimate of inflation at each time step
  infl_est2 <- sqrt(pmax(0, sapply(1:reps, function(i){mean(var_pivots[1:(i*n_folds), 1]^2 - var_pivots[1:(i*n_folds), 2])}))) /
    (sd(ho_errs) / sqrt(n_sub))
  
  #bias correction
  cv_means <- c() #estimated pred error from normal CV
  bias_est <- 0
  if(is.na(bias_reps)) {
    bias_reps <- ceiling(reps / 5) #fewer reps for bias estimation
  }
  if(bias_reps == 0) {
    bias_est <- 0
  }
  else {
    for(i in 1:bias_reps) {
      temp <- naive_cv_m(X_m, Y_m, Treat, range.tau, funcs, n_folds, funcs_params = funcs_params, alpha = alpha)
      cv_means <- c(cv_means, temp$err_hat)
    }
    
    bias_est <- (mean(ho_errs) - mean(cv_means)) * (1 + ((n_folds - 2) / (n_folds ))^(1.5))
  }
  pred_est <- mean(ho_errs) - bias_est #debiased estimate
  
  results <- list()
  
  for (a in alpha) {
    ci_lo <- pred_est - qnorm(1 - a / 2) * sd(ho_errs) / sqrt(length(Y)) * ugp_infl
    ci_hi <- pred_est + qnorm(1 - a / 2) * sd(ho_errs) / sqrt(length(Y)) * ugp_infl
    
    results[[paste0("ci_lo_alpha_", a)]] <- ci_lo
    results[[paste0("ci_hi_alpha_", a)]] <- ci_hi
  }
  
  results[["sd_infl"]] <- ugp_infl
  results[["err_hat"]] <- pred_est
  results[["raw_mean"]] <- mean(ho_errs)
  results[["bias_est"]] <- bias_est
  results[["sd"]] <- sd(ho_errs)
  results[["running_sd_infl"]] <- infl_est2
  
  return(results)
}

nested_cv_helper_m <- function(X_m, Y_m, Treat, range.tau, funcs, n_folds, funcs_params = NULL) {
  min_ho_errors <- NULL
  fold_id <- 1:nrow(X_m) %% n_folds + 1
  fold_id <- sample(fold_id[1:(nrow(X_m) %/% n_folds * n_folds)]) #handle case where n doesn't divide n_folds
  fold_id <- c(fold_id, rep(0, nrow(X_m) %% n_folds))
  
  #nested CV model fitting
  ho_errors <- array(0, dim = c(n_folds, n_folds, nrow(X_m) %/% n_folds))
  #^ entry i, j is error on fold i,
  # when folds i & j are not used for model fitting.
  for (tau in range.tau) {
  for(f1 in 1:(n_folds - 1)) {
    for(f2 in (f1+1):n_folds) {
        test_idx <- c(which(fold_id == f1), which(fold_id == f2))
        fit <- funcs$fitter(as.matrix(X_m[-test_idx, ]), Y_m[-test_idx], Treat[-test_idx], tau, funcs_params = funcs_params)
        preds <- funcs$predictor(fit, X_m, funcs_params = funcs_params)
        ho_errors[f1, f2, ] <- funcs$loss(preds[fold_id == f1], Y_m[fold_id == f1], Treat[fold_id == f1], tau, funcs_params = funcs_params)
        ho_errors[f2, f1, ] <- funcs$loss(preds[fold_id == f2], Y_m[fold_id == f2], Treat[fold_id == f2], tau, funcs_params = funcs_params)
      }
      
  }
    if (is.null(min_ho_errors) || sum(ho_errors) < sum(min_ho_errors)) {
      min_ho_errors <- ho_errors
      best_tau <- tau
    }
  }
  
  #e_bar - f_bar in the notation of the paper. n_folds X_m 1 matrix
  out_mat <- matrix(0, n_folds, 2)
  for(f1 in 1:(n_folds)) {
    test_idx <- which(fold_id == f1)
    fit <- funcs$fitter(as.matrix(X_m[-test_idx, ]), Y_m[-test_idx], Treat[-test_idx], best_tau, funcs_params = funcs_params)
    preds <- funcs$predictor(fit, X_m[test_idx, ], funcs_params = funcs_params)
    e_out <- funcs$loss(preds, Y_m[test_idx], Treat[test_idx], best_tau)
    
    #loop over other folds
    e_bar_t <- c() # errors from internal CV
    for(f2 in 1:n_folds) {
      if(f2 == f1) {next}
      e_bar_t <- c(e_bar_t, min_ho_errors[f2, f1, ])
    }
    out_mat[f1, 1] <- mean(e_bar_t) - mean(e_out) # (a) terms
    out_mat[f1, 2] <- var(e_out) / length(test_idx) # (b) terms
  }
  
  #errors on points not used for fitting, combined across all runs (used for point estimate)
  all_ho_errs <- c()
  for(f1 in 1:(n_folds - 1)) {
    for(f2 in (f1+1):n_folds) {
      all_ho_errs <- c(all_ho_errs, min_ho_errors[f1, f2, ], min_ho_errors[f2, f1, ])
    }
  }
  
  return(list("pivots" = out_mat,
              "errs" = all_ho_errs,
              "tau" = best_tau))
}



naive_cv_BART <- function(X, Y, Treat,funcs, n_folds, alpha = alpha,
                       trans = list(identity), funcs_params = NULL, fold_id = NULL) {
  min_gp_errors <- NULL
  if(is.null(fold_id)) {
    fold_id <- (1:nrow(X)) %% n_folds + 1
    fold_id <- sample(fold_id)
  }
  
  errors <- c()
  gp_errors <- c()
    for(k in 1:n_folds) {
      fit <- funcs$fitter(X[fold_id !=k, ], Y[fold_id != k], Treat[fold_id != k], funcs_params = funcs_params)
      Y_hat <- funcs$predictor(fit, X[fold_id == k, ], Treat[fold_id == k], funcs_params = funcs_params)
      error_k <- funcs$loss(Y_hat, Y[fold_id == k], funcs_params = funcs_params)
      errors <- c(errors, error_k)
      
      temp_vec <- c()
      for(tran in trans) {
        temp_vec <- c(temp_vec, tran(mean(error_k)))
      }
      gp_errors <- rbind(gp_errors, temp_vec)
    }
    if (is.null(min_gp_errors) || sum(gp_errors) < sum(min_gp_errors)) {
      min_gp_errors <- gp_errors
      best_errors <- errors
    }

  return(list("err_hat" = mean(best_errors),
              "ci_lo" = mean(best_errors) - qnorm(1-alpha/2) * sd(best_errors) / sqrt(length(Y)),
              "ci_hi" = mean(best_errors) + qnorm(1-alpha/2) * sd(best_errors) / sqrt(length(Y)),
              "raw_mean" = mean(best_errors),
              "sd" = sd(best_errors),
              "group_err_hat" = apply(min_gp_errors, 2, mean),
              "group_sd" = apply(min_gp_errors, 2, sd),
              "raw_errors" = best_errors,
              "fold_id" = fold_id))
}


# nested CV

nested_cv_BART <- function(X, Y, Treat,  funcs, reps, n_folds,  alpha = alpha, bias_reps = NA,
                        funcs_params = NULL, n_cores = 1, verbose = F) {
  #estimate time required
  if(verbose) {
    t1 <- Sys.time()
    temp <- nested_cv_helper_BART(X, Y, Treat,  funcs, n_folds,
                               funcs_params = funcs_params)
    t2 <- Sys.time()
    print(paste0("Estimated time required: ", (t2 - t1) * reps))
  }
  
  #compute out-of-fold errors on SE scale
  var_pivots <- c()
  gp_errs <- c()
  ho_errs <- c()
  if(n_cores == 1){
    raw <- lapply(1:reps, function(i){nested_cv_helper_BART(X, Y, Treat,funcs, n_folds,
                                                         funcs_params = funcs_params)})
  } else {
    raw <- parallel::mclapply(1:reps, function(i){nested_cv_helper_BART(X, Y, Treat, funcs, n_folds,
                                                                     funcs_params = funcs_params)},
                              mc.cores = n_cores)
  }
  for(i in 1:reps) {
    temp <- raw[[i]]
    var_pivots <- rbind(var_pivots, temp$pivots)
    ho_errs <- c(ho_errs, temp$errs)
  }
  
  ########
  n_sub <- floor(length(Y) * (n_folds - 1) / n_folds)
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
      temp <- naive_cv_BART(X, Y, Treat, funcs, n_folds, funcs_params = funcs_params, alpha = alpha)
      cv_means <- c(cv_means, temp$err_hat)
    }
    
    bias_est <- (mean(ho_errs) - mean(cv_means)) * (1 + ((n_folds - 2) / (n_folds ))^(1.5))
  }
  pred_est <- mean(ho_errs) - bias_est #debiased estimate
  
  list("sd_infl" = ugp_infl,
       "err_hat" = pred_est,
       "ci_lo" = pred_est - qnorm(1-alpha/2) * sd(ho_errs) / sqrt(length(Y)) * ugp_infl,
       "ci_hi" = pred_est + qnorm(1-alpha/2) * sd(ho_errs) / sqrt(length(Y)) * ugp_infl,
       "raw_mean" = mean(ho_errs),
       "bias_est" = bias_est,
       "sd" = sd(ho_errs),
       "running_sd_infl" = infl_est2)
}

nested_cv_helper_BART <- function(X, Y, Treat, funcs, n_folds, funcs_params = NULL) {
  min_ho_errors <- NULL
  fold_id <- 1:nrow(X) %% n_folds + 1
  fold_id <- sample(fold_id[1:(nrow(X) %/% n_folds * n_folds)]) #handle case where n doesn't divide n_folds
  fold_id <- c(fold_id, rep(0, nrow(X) %% n_folds))
  
  #nested CV model fitting
  ho_errors <- array(0, dim = c(n_folds, n_folds, nrow(X) %/% n_folds))
  #^ entry i, j is error on fold i,
  # when folds i & j are not used for model fitting.
    for(f1 in 1:(n_folds - 1)) {
      for(f2 in (f1+1):n_folds) {
        test_idx <- c(which(fold_id == f1), which(fold_id == f2))
        fit <- funcs$fitter(as.matrix(X[-test_idx, ]), Y[-test_idx], Treat[-test_idx], funcs_params = funcs_params)
        preds <- funcs$predictor(fit, X, Treat, funcs_params = funcs_params)
        ho_errors[f1, f2, ] <- funcs$loss(preds[fold_id == f1], Y[fold_id == f1], funcs_params = funcs_params)
        ho_errors[f2, f1, ] <- funcs$loss(preds[fold_id == f2], Y[fold_id == f2], funcs_params = funcs_params)
      }
      
    }
    if (is.null(min_ho_errors) || sum(ho_errors) < sum(min_ho_errors)) {
      min_ho_errors <- ho_errors
    }

  
  #e_bar - f_bar in the notation of the paper. n_folds X 1 matrix
  out_mat <- matrix(0, n_folds, 2)
  for(f1 in 1:(n_folds)) {
    test_idx <- which(fold_id == f1)
    fit <- funcs$fitter(as.matrix(X[-test_idx, ]), Y[-test_idx], Treat[-test_idx], funcs_params = funcs_params)
    preds <- funcs$predictor(fit, X[test_idx, ], Treat[test_idx], funcs_params = funcs_params)
    e_out <- funcs$loss(preds, Y[test_idx])
    
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
              "errs" = all_ho_errs
              ))
}

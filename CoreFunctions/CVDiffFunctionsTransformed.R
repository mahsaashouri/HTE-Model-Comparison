

naive_cv <- function(X, W, funcs, n_folds, alpha = c(0.01, 0.05, 0.1, 0.25, 0.5),
                     trans = list(identity), fold_id = NULL) {
  min_gp_errors <- NULL
  if(is.null(fold_id)) {
    fold_id <- (1:nrow(X)) %% n_folds + 1
    fold_id <- sample(fold_id)
  }

  errors <- c()
  mse_full <- c()
  mse_reduced <- c()
  gp_errors <- c()
  for(k in 1:n_folds) {

    fit <- funcs$fitter(as.matrix(X[fold_id !=k, ]), W[fold_id != k])
    preds <- funcs$predictor(fit, X[fold_id == k, ])
    error_k <- funcs$loss(preds$pred_full, preds$tau_star, W[fold_id == k])
    errors <- c(errors, error_k)
    temp_vec <- c()
    for(tran in trans) {
      temp_vec <- c(temp_vec, tran(mean(error_k)))
    }
    gp_errors <- rbind(gp_errors, temp_vec)
  }

  return(list("err_hat" = mean(errors),
              "ci_lo" = mean(errors) - qnorm(1-alpha/2) * sd(errors) / sqrt(length(W)),
              "ci_hi" = mean(errors) + qnorm(1-alpha/2) * sd(errors) / sqrt(length(W)),
              "raw_mean" = mean(errors),
              "sd" = sd(errors),
              "group_err_hat" = apply(gp_errors, 2, mean),
              "group_sd" = apply(gp_errors, 2, sd),
              "raw_errors" = errors,
              "fold_id" = fold_id))

}


# nested CV

nested_cv <- function(X, W, funcs, reps, n_folds,  alpha = c(0.01, 0.05, 0.1, 0.25, 0.5),
                      bias_reps = NA, n_cores = 1) {

  #compute out-of-fold errors on SE scale
  var_pivots <- c()
  gp_errs <- c()
  ho_errs <- c()
  mse_full <- c()
  mse_reduced <- c()

  #if(n_cores == 1){
  raw <- lapply(1:reps, function(i){nested_cv_helper(X, W, funcs, n_folds)})
  #}
  for(i in 1:reps) {
    temp <- raw[[i]]
    var_pivots <- rbind(var_pivots, temp$pivots)
    ho_errs <- c(ho_errs, temp$errs)
  }
  n_sub <- floor(length(W) * (n_folds - 1) / n_folds)
  # look at the estimate of inflation after each repetition
  ugp_infl <- sqrt(max(0, mean(var_pivots[, 1]^2 - var_pivots[, 2]))) / (sd(ho_errs) / sqrt(n_sub))
  ugp_infl <- max(1, min(ugp_infl, sqrt(n_folds)))

  #estimate of inflation at each time step
  infl_est2 <- sqrt(pmax(0, sapply(1:reps, function(i){mean(var_pivots[1:(i*n_folds), 1]^2 - var_pivots[1:(i*n_folds), 2])}))) /
    (sd(ho_errs) / sqrt(n_sub))

  #bias correction
  cv_means <- c() #estimated pred error from normal CV
  bias_est <- 0
  print(n_folds)
  if(is.na(bias_reps)) {
    bias_reps <- ceiling(reps / 5) #fewer reps for bias estimation
  }
  if(bias_reps == 0) {
    bias_est <- 0
  }
  else {
    for(i in 1:bias_reps) {
      temp <- naive_cv(X, W, funcs, n_folds)
      cv_means <- c(cv_means, temp$err_hat)
    }

    bias_est <- (mean(ho_errs) - mean(cv_means)) * (1 + ((n_folds - 2) / (n_folds ))^(1.5))
  }
  pred_est <- mean(ho_errs) - bias_est #debiased estimate

  results <- list()
  zero_between_bounds <- rep(NA, length(alpha))
  std.err <- sd(ho_errs) / sqrt(length(W)) * ugp_infl
  results[["ci_lo"]] <- pred_est - std.err*qnorm(1 - 0.05/2)
  results[["ci_hi"]] <- pred_est + std.err*qnorm(1 - 0.05/2)
  results[["hvalue"]] <- 2*pnorm(-abs(pred_est)/std.err)
  results[["hvalue1sided"]] <- pnorm(pred_est/std.err, lower.tail=FALSE)
  results[["sd_infl"]] <- ugp_infl
  results[["err_hat"]] <- pred_est
  results[["raw_mean"]] <- mean(ho_errs)
  results[["bias_est"]] <- bias_est
  results[["sd"]] <- sd(ho_errs)
  results[["running_sd_infl"]] <- infl_est2
  #results[["prop_zero_CI"]] <- mean(zero_between_bounds)
  #results[["mse_full"]] <- mse_full
  # results[["mse_reduced"]] <- mse_reduced

  return(results)
}

nested_cv_helper <- function(X, W, funcs, n_folds = 10, funcs_params = NULL) {
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
      fit <- funcs$fitter(X[-test_idx, ], W[-test_idx])
      preds <- funcs$predictor(fit, X)
      ho_errors[f1, f2, ] <- funcs$loss(preds$pred_full[fold_id == f1], preds$tau_star, W[fold_id == f1])
      ho_errors[f2, f1, ] <- funcs$loss(preds$pred_full[fold_id == f2], preds$tau_star, W[fold_id == f2])
    }
  }

  #e_bar - f_bar in the notation of the paper. n_folds x 1 matrix
  out_mat <- matrix(0, n_folds, 2)
  for(f1 in 1:(n_folds)) {
    test_idx <- which(fold_id == f1)
    fit <- funcs$fitter(X[-test_idx, ], W[-test_idx])
    preds <- funcs$predictor(fit, X[test_idx, ])
    e_out <- funcs$loss(preds$pred_full, preds$tau_star, W[test_idx]) ## Change this

    #loop over other folds
    e_bar_t <- c() # errors from internal CV
    for(f2 in 1:n_folds) {
      if(f2 == f1) {next}
      e_bar_t <- c(e_bar_t, ho_errors[f2, f1, ])
    }
    out_mat[f1, 1] <- mean(e_bar_t) - mean(e_out) # (a) terms
    out_mat[f1, 2] <- var(e_out) / length(test_idx) # (b) terms
  }

  #errors on points not used for fitting, combined across all runs (used for point estimate)
  all_ho_errs <- c()
  for(f1 in 1:(n_folds - 1)) {
    for(f2 in (f1+1):n_folds) {
      all_ho_errs <- c(all_ho_errs, ho_errors[f1, f2, ], ho_errors[f2, f1, ])
    }
  }

  return(list("pivots" = out_mat,
              "errs" = all_ho_errs))
}



naive_cv <- function(X, X0, Y, Trt, tau.seq, funcs, n_folds, alpha = c(0.01, 0.05, 0.1, 0.25, 0.5),
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
      fit <- funcs$fitter(as.matrix(X[fold_id !=k, ]), as.matrix(X0[fold_id !=k, ]), Y[fold_id != k], Trt[fold_id != k], tau.seq)
      Y_hat <- funcs$predictor(fit$full, X[fold_id == k, ])
      Y_hat_reduced <- funcs$predictor(fit$reduced, X0[fold_id == k, ])
      error_k <- funcs$loss(Y_hat, Y_hat_reduced, Y[fold_id == k], Trt[fold_id == k], fit$tau)
      mse_k <- funcs$mse(Y_hat, Y_hat_reduced, Y[fold_id == k], Trt[fold_id == k], fit$tau)
      errors <- c(errors, error_k)
      mse_full <- c(mse_full, mse_k$full)
      mse_reduced <- c(mse_reduced, mse_k$reduced)
      
      temp_vec <- c()
      for(tran in trans) {
        temp_vec <- c(temp_vec, tran(mean(error_k)))
      }
      gp_errors <- rbind(gp_errors, temp_vec)
    }
  
  results <- list()
  for (a in alpha) {
    ci_lo <- mean(errors) - qnorm(1 - a / 2) * sd(errors) / sqrt(length(Y))
    ci_hi <- mean(errors) + qnorm(1 - a / 2) * sd(errors) / sqrt(length(Y))
    
    results[[paste0("ci_lo_alpha_", a)]] <- ci_lo
    results[[paste0("ci_hi_alpha_", a)]] <- ci_hi
  }
  
  results[["mse_full"]] <- mean(mse_full)
  results[["mse_reduced"]] <- mean(mse_reduced)
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

nested_cv <- function(X, X0, Y, Trt, tau.seq, funcs, reps, n_folds,  alpha = c(0.01, 0.05, 0.1, 0.25, 0.5), bias_reps = NA, n_cores = 1, verbose = F) {
  #estimate time required
  if(verbose) {
    t1 <- Sys.time()
    temp <- nested_cv_helper(X, X0, Y, Trt, tau.seq, funcs, n_folds)
    t2 <- Sys.time()
    print(paste0("Estimated time required: ", (t2 - t1) * reps))
  }
  
  #compute out-of-fold errors on SE scale
  var_pivots <- c()
  gp_errs <- c()
  ho_errs <- c()
  mse_full <- c()
  mse_reduced <- c()
  
  if(n_cores == 1){
    raw <- lapply(1:reps, function(i){nested_cv_helper(X, X0, Y, Trt, tau.seq, funcs, n_folds)})
  } else {
    raw <- parallel::mclapply(1:reps, function(i){nested_cv_helper(X, X0, Y, Trt, tau.seq, funcs, n_folds)},
                              mc.cores = n_cores)
  }
  for(i in 1:reps) {
    temp <- raw[[i]]
    var_pivots <- rbind(var_pivots, temp$pivots)
    ho_errs <- c(ho_errs, temp$errs)
    mse_full <- c(mse_full, temp$mse_full)
    mse_reduced <- c(mse_reduced, temp$mse_reduced)
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
      temp <- naive_cv(X, X0, Y, Trt, tau.seq, funcs, n_folds, alpha = alpha)
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
  results[["mse_full"]] <- mse_full
  results[["mse_reduced"]] <- mse_reduced
  
  return(results)
}

nested_cv_helper <- function(X, X0, Y, Trt, tau.seq, funcs, n_folds) {
  fold_id <- 1:nrow(X) %% n_folds + 1
  fold_id <- sample(fold_id[1:(nrow(X) %/% n_folds * n_folds)]) #handle case where n doesn't divide n_folds
  fold_id <- c(fold_id, rep(0, nrow(X) %% n_folds))
  
  #nested CV model fitting
  ho_errors <- array(0, dim = c(n_folds, n_folds, nrow(X) %/% n_folds))
  ho_mse_full <- matrix(0, ncol = n_folds, nrow = n_folds)
  ho_mse_reduced <- matrix(0, ncol = n_folds, nrow = n_folds)
  #^ entry i, j is error on fold i,
  # when folds i & j are not used for model fitting.
    for(f1 in 1:(n_folds - 1)) {
      for(f2 in (f1+1):n_folds) {
        test_idx <- c(which(fold_id == f1), which(fold_id == f2))
        fit <- funcs$fitter(as.matrix(X[-test_idx, ]), as.matrix(X0[-test_idx, ]), Y[-test_idx], Trt[-test_idx], tau.seq)
        preds <- funcs$predictor(fit$full, X)
        preds_reduced <- funcs$predictor(fit$reduced, X0)
        ho_errors[f1, f2, ] <- funcs$loss(preds[fold_id == f1], preds_reduced[fold_id == f1], Y[fold_id == f1], Trt[fold_id == f1], fit$tau)
        ho_errors[f2, f1, ] <- funcs$loss(preds[fold_id == f2], preds_reduced[fold_id == f1], Y[fold_id == f2], Trt[fold_id == f2], fit$tau)
        mse_full_reduced12 <- funcs$mse(preds[fold_id == f1], preds_reduced[fold_id == f1], Y[fold_id == f1], Trt[fold_id == f1], fit$tau)
        mse_full_reduced21 <- funcs$mse(preds[fold_id == f2], preds_reduced[fold_id == f1], Y[fold_id == f2], Trt[fold_id == f2], fit$tau)
        ho_mse_full[f1, f2] <- mse_full_reduced12$full; ho_mse_full[f2, f1] <- mse_full_reduced21$full
        ho_mse_reduced[f1, f2] <- mse_full_reduced12$reduced;  ho_mse_reduced[f2, f1] <- mse_full_reduced21$reduced
      }
      
    }
  
  #e_bar - f_bar in the notation of the paper. n_folds X 1 matrix
  out_mat <- matrix(0, n_folds, 2)
  for(f1 in 1:(n_folds)) {
    test_idx <- which(fold_id == f1)
    fit <- funcs$fitter(as.matrix(X[-test_idx, ]), as.matrix(X0[-test_idx, ]), Y[-test_idx], Trt[-test_idx], tau.seq)
    preds <- funcs$predictor(fit$full, X[test_idx, ])
    preds_reduced <- funcs$predictor(fit$reduced, X0[test_idx, ])
    e_out <- funcs$loss(preds, preds_reduced, Y[test_idx], Trt[test_idx], fit$tau)
    
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
              "errs" = all_ho_errs,
              "mse_full" = mean(as.vector(ho_mse_full)),
              "mse_reduced" = mean(as.vector(ho_mse_reduced))
              ))
}

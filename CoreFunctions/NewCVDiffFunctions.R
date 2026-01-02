# Load required library for parallel processing
library(parallel)

ncv_nb <- function(X, X0, Y, Trt, tau.range, funcs, n_folds, alpha = 0.05,
                   trans = list(identity)) {
  n <- length(Y)
  fold_id <- sample(1:n_folds, size=n, replace=TRUE)
  
  ehat <- rep(NA, n_folds)
  errors <- c()
  for(k in 1:n_folds) {
    
    fit <- funcs$fitter(as.matrix(X[fold_id !=k, ]), as.matrix(X0[fold_id !=k, ]), Y[fold_id != k],
                        Trt[fold_id != k], tau.range)
    preds <- funcs$predictor(fit, X[fold_id == k, ], X0[fold_id==k,], Trt[fold_id==k])
    
    error_k <- funcs$loss(preds$pred_full, preds$pred_reduced, Y[fold_id == k])
    ehat[k] <- mean(error_k)
    errors <- c(errors, error_k)
  }
  err_hat <- mean(ehat)
  se_naive <- sd(errors)/sqrt(n)
  
  split_id <- sample(0:1, size=n, replace=TRUE)
  ehat0 <- ehat1 <- rep(NA, n_folds)
  for(k in 1:n_folds) {
    
    fit <- funcs$fitter(as.matrix(X[fold_id !=k & split_id==0, ]), as.matrix(X0[fold_id !=k & split_id==0, ]), 
                        Y[fold_id != k & split_id==0], Trt[fold_id != k & split_id==0], tau.range)
    preds <- funcs$predictor(fit, X[fold_id == k & split_id==0, ], X0[fold_id==k & split_id==0,], Trt[fold_id==k & split_id==0])
    
    error_k <- funcs$loss(preds$pred_full, preds$pred_reduced, Y[fold_id == k & split_id==0])
    ehat0[k] <- mean(error_k)
  }
  for(k in 1:n_folds) {
    
    fit <- funcs$fitter(as.matrix(X[fold_id !=k & split_id==1, ]), as.matrix(X0[fold_id !=k & split_id==1, ]), 
                        Y[fold_id != k & split_id==1], Trt[fold_id != k & split_id==1], tau.range)
    preds <- funcs$predictor(fit, X[fold_id == k & split_id==1, ], X0[fold_id==k & split_id==1,], Trt[fold_id==k & split_id==1])
    
    error_k <- funcs$loss(preds$pred_full, preds$pred_reduced, Y[fold_id == k & split_id==1])
    ehat1[k] <- mean(error_k)
  }
  err_hat <- mean(ehat)
  se_nb <- sqrt(0.5*((mean(ehat0) - mean(ehat1))^2))
  
  print(c(se_nb, se_naive))
  std_err <- max(c(se_naive, se_nb))
  return(list("err_hat" = err_hat,
              "std_err" = std_err,
              "ci_lo" = err_hat - qnorm(1-alpha/2) * std_err,
              "ci_hi" = err_hat + qnorm(1-alpha/2) * std_err))
}



data_splitting <- function(X, X0, Y, Trt, tau.range, funcs, train_prob, alpha = 0.05,
                           trans = list(identity)) {
  min_gp_errors <- NULL
  n <- length(Y)
  train_ind <- sample(0:1, size=n, replace=TRUE, prob=c(train_prob, 1 - train_prob))
  
  
  fit <- funcs$fitter(as.matrix(X[train_ind ==1, ]), as.matrix(X0[train_ind ==1, ]), Y[train_ind == 1],
                      Trt[train_ind == 1], tau.range)
  preds <- funcs$predictor(fit, X[train_ind == 0, ], X0[train_ind==0,], Trt[train_ind==0])
  
  errors <- funcs$loss(preds$pred_full, preds$pred_reduced, Y[train_ind==0])
  
  err_hat <- mean(errors)
  std_err <- sd(errors)/sqrt(sum(train_ind==0))
  
  return(list("err_hat" = err_hat,
              "std_err" = std_err,
              "ci_lo" = err_hat - qnorm(1-alpha/2) * std_err,
              "ci_hi" = err_hat + qnorm(1-alpha/2) * std_err))
}


naive_cv <- function(X, X0, Y, Trt, tau.range, funcs, n_folds, alpha = 0.05,
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
    
    fit <- funcs$fitter(as.matrix(X[fold_id !=k, ]), as.matrix(X0[fold_id !=k, ]), Y[fold_id != k],
                        Trt[fold_id != k], tau.range)
    preds <- funcs$predictor(fit, X[fold_id == k, ], X0[fold_id==k,], Trt[fold_id==k])
    
    error_k <- funcs$loss(preds$pred_full, preds$pred_reduced, Y[fold_id == k])
    
    errors <- c(errors, error_k)
    temp_vec <- c()
    for(tran in trans) {
      temp_vec <- c(temp_vec, tran(mean(error_k)))
    }
    gp_errors <- rbind(gp_errors, temp_vec)
  }
  
  return(list("err_hat" = mean(errors),
              "ci_lo" = mean(errors) - qnorm(1-alpha/2) * sd(errors) / sqrt(length(Y)),
              "ci_hi" = mean(errors) + qnorm(1-alpha/2) * sd(errors) / sqrt(length(Y)),
              "raw_mean" = mean(errors),
              "sd" = sd(errors),
              "group_err_hat" = apply(gp_errors, 2, mean),
              "group_sd" = apply(gp_errors, 2, sd),
              "raw_errors" = errors,
              "fold_id" = fold_id))
  
}


# nested CV with parallel processing support

nested_cv <- function(X, X0, Y, Trt, tau.range, funcs, reps, n_folds,  
                      alpha = c(0.01, 0.05, 0.1, 0.25, 0.5),
                      bias_reps = NA, n_cores = 4) {
  
  # Compute out-of-fold errors on SE scale
  var_pivots <- c()
  gp_errs <- c()
  ho_errs <- c()
  mse_full <- c()
  mse_reduced <- c()
  
  # Parallel or sequential processing
  if(n_cores > 1) {
    # Setup parallel cluster
    cl <- makeCluster(n_cores)
    on.exit(stopCluster(cl), add = TRUE)  # Ensure cleanup
    
    # Export necessary objects to cluster
    clusterExport(cl, c("X", "X0", "Y", "Trt", "tau.range", "funcs", "n_folds",
                        "nested_cv_helper"), 
                  envir = environment())
    
    # Load required libraries on each worker
    clusterEvalQ(cl, {
      library(glmnet)
      library(mboost)
      library(dbarts)
      library(randomForest)
      library(bcf)
      library(grf)
    })
    
    # Run parallel computation
    raw <- parLapply(cl, 1:reps, function(i) {
      nested_cv_helper(X, X0, Y, Trt, tau.range, funcs, n_folds)
    })
    
  } else {
    # Sequential processing
    raw <- lapply(1:reps, function(i) {
      nested_cv_helper(X, X0, Y, Trt, tau.range, funcs, n_folds)
    })
  }
  
  # Process results
  for(i in 1:reps) {
    temp <- raw[[i]]
    var_pivots <- rbind(var_pivots, temp$pivots)
    ho_errs <- c(ho_errs, temp$errs)
  }
  
  n_sub <- floor(length(Y) * (n_folds - 1) / n_folds)
  
  # Look at the estimate of inflation after each repetition
  ugp_infl <- sqrt(max(0, mean(var_pivots[, 1]^2 - var_pivots[, 2]))) / (sd(ho_errs) / sqrt(n_sub))
  ugp_infl <- max(1, min(ugp_infl, sqrt(n_folds)))
  
  # Estimate of inflation at each time step
  infl_est2 <- sqrt(pmax(0, sapply(1:reps, function(i){
    mean(var_pivots[1:(i*n_folds), 1]^2 - var_pivots[1:(i*n_folds), 2])
  }))) / (sd(ho_errs) / sqrt(n_sub))
  
  # Bias correction
  cv_means <- c() # Estimated pred error from normal CV
  bias_est <- 0
  if(is.na(bias_reps)) {
    bias_reps <- ceiling(reps / 5) # Fewer reps for bias estimation
  }
  if(bias_reps == 0) {
    bias_est <- 0
  } else {
    for(i in 1:bias_reps) {
      temp <- naive_cv(X, X0, Y, Trt, tau.range, funcs, n_folds)
      cv_means <- c(cv_means, temp$err_hat)
    }
    
    bias_est <- (mean(ho_errs) - mean(cv_means)) * (1 + ((n_folds - 2) / (n_folds ))^(1.5))
  }
  pred_est <- mean(ho_errs) - bias_est # Debiased estimate
  
  results <- list()
  zero_between_bounds <- rep(NA, length(alpha))
  std.err <- sd(ho_errs) / sqrt(length(Y)) * ugp_infl
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
  
  return(results)
}

nested_cv_helper <- function(X, X0, Y, Trt, tau.range, funcs, n_folds = 10, funcs_params = NULL) {
  fold_id <- 1:nrow(X) %% n_folds + 1
  fold_id <- sample(fold_id[1:(nrow(X) %/% n_folds * n_folds)]) # Handle case where n doesn't divide n_folds
  fold_id <- c(fold_id, rep(0, nrow(X) %% n_folds))
  
  # Nested CV model fitting
  ho_errors <- array(0, dim = c(n_folds, n_folds, nrow(X) %/% n_folds))
  # ^ entry i, j is error on fold i,
  # when folds i & j are not used for model fitting. 
  for(f1 in 1:(n_folds - 1)) {
    for(f2 in (f1+1):n_folds) {
      test_idx <- c(which(fold_id == f1), which(fold_id == f2))
      fit <- funcs$fitter(X[-test_idx, ], X0[-test_idx,], Y[-test_idx], Trt[-test_idx],
                          tau.range)
      preds <- funcs$predictor(fit, X, X0, Trt)
      ho_errors[f1, f2, ] <- funcs$loss(preds$pred_full[fold_id == f1], preds$pred_reduced[fold_id == f1], Y[fold_id == f1])
      ho_errors[f2, f1, ] <- funcs$loss(preds$pred_full[fold_id == f2], preds$pred_reduced[fold_id == f2], Y[fold_id == f2])
    }
  }
  
  # e_bar - f_bar in the notation of the paper.  n_folds x 1 matrix
  out_mat <- matrix(0, n_folds, 2)
  for(f1 in 1:(n_folds)) {
    test_idx <- which(fold_id == f1)
    fit <- funcs$fitter(X[-test_idx, ], X0[-test_idx,], Y[-test_idx], Trt[-test_idx],
                        tau.range)
    preds <- funcs$predictor(fit, X[test_idx, ], X0[test_idx, ], Trt[test_idx])
    e_out <- funcs$loss(preds$pred_full, preds$pred_reduced, Y[test_idx])
    
    # Loop over other folds
    e_bar_t <- c() # Errors from internal CV
    for(f2 in 1:n_folds) {
      if(f2 == f1) {next}
      e_bar_t <- c(e_bar_t, ho_errors[f2, f1, ])
    }
    out_mat[f1, 1] <- mean(e_bar_t) - mean(e_out) # (a) terms
    out_mat[f1, 2] <- var(e_out) / length(test_idx) # (b) terms
  }
  
  # Errors on points not used for fitting, combined across all runs (used for point estimate)
  all_ho_errs <- c()
  for(f1 in 1:(n_folds - 1)) {
    for(f2 in (f1+1):n_folds) {
      all_ho_errs <- c(all_ho_errs, ho_errors[f1, f2, ], ho_errors[f2, f1, ])
    }
  }
  
  return(list("pivots" = out_mat,
              "errs" = all_ho_errs))
}
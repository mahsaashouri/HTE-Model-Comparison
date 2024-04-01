
squared_loss <- function(y1, y2, y3, tau, Trt, funcs_params = NA) {
  ## y1 - full model
  ## y2 - reduced model 
  ## y3 - outcome
  (y1 - y3)^2 - (y2 - (y3 - tau*Trt))^2
}

fitter_glmnet <- function(X, X0, Y, Trt, tau.seq, idx = NA, funcs_params = NA) {
  ## X0 does not have treatment column
  if(sum(is.na(idx)) > 0) {
    idx <- 1:nrow(X)
  }
  fit <- glmnet(X[idx, ], Y[idx], family = "gaussian", lambda = funcs_params$lambdas) 
  
  mse <- rep(NA, length(tau.seq))
  for(k in 1:length(tau.seq)) {
    glmnet_tmp <- glmnet(X0[idx, ], Y[idx]-(tau.seq[k]*Treat[idx]), family = "gaussian", lambda = funcs_params$lambdas) 
    mse[k] <- mean(((Y[idx]-(tau.seq[k]*Treat[idx])) - glmnet_tmp$fitted)^2)
  }
  tau.star <- tau.seq[which.min(mse)]
  
  fit_reduced <- glmnet(X0[idx, ], Y[idx]-(tau.star*Treat[idx]), family = "gaussian", lambda = funcs_params$lambdas) 
  
  return(list(full=fit, reduced=fit_reduced, tau = tau.star))
}

predictor_glmnet <- function(fit, X_new, funcs_params = NA) {
  beta_hat <- fit$beta[, funcs_params$best_lam] 
  a0_hat <- fit$a0[funcs_params$best_lam]
  preds <- (as.matrix(X_new) %*% beta_hat + a0_hat)
  preds
}

gaussian_lasso_funs <- list(fitter = fitter_glmnet,
                            predictor = predictor_glmnet,
                            loss = squared_loss,
                            name = "gaussian_lasso")


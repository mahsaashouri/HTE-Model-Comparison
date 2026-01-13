library(glmnet)     # Already loaded, but ridge uses this too
library(bcf)        # For Bayesian Causal Forest
library(grf)        # For Causal Forest
library(bartCause)  # For bartCause
#library(devtools) 
#install_github("xnie/rlearner")
library(rlearner)   # For rlearner

# Your existing functions remain the same
squared_loss <- function(fhat_full, fhat_reduced, y) {
  ## y1 - full model
  ## y2 - reduced model
  ## y3 - outcome
  (fhat_reduced - y)^2 - (fhat_full - y)^2
}

mse <- function(y1, y2, y3, Trt, tau) {
  ## y1 - full model
  ## y2 - reduced model
  ## y3 - outcome
  mse_full <- mean((y1 - y3)^2)
  mse_reduced <- mean((y2 - (y3 - tau*Trt))^2)
  return(list(full = mse_full, reduced = mse_reduced))
}

### Fitting Functions 
fitter_lm <- function(X, X0, Y, Trt, tau.range, idx = NA) {
  ## X0 does not have treatment column
  if(sum(is.na(idx)) > 0) {
    idx <- 1:nrow(X)
  }
  data.all <- cbind.data.frame(X[idx, ], Y = Y[idx])
  fit <- lm(Y ~., data = data.all)
  
  f0fn <- function(tau) {
    data.reduced <- cbind.data.frame(X0[idx, ], Y= Y[idx] - tau*Trt[idx])
    fit_reduced <- lm(Y ~., data = data.reduced)
    mse_local <- mean((data.reduced$Y - fit_reduced$fitted.values )^2)
    return(mse_local)
  }
  tau.star <- optimize(f0fn, interval=tau.range)$minimum
  data.reduced <- cbind.data.frame(X0[idx, ], Y= Y[idx] - tau.star*Trt[idx])
  fit_reduced <- lm(Y ~., data = data.reduced)
  return(list(full=fit, reduced=fit_reduced, tau = tau.star))
}


fitter_glmnet <- function(X, X0, Y, Trt, tau.range, idx = NA) {
  ## X0 does not have treatment column
  if(sum(is.na(idx)) > 0) {
    idx <- 1:nrow(X)
  }
  Xmat <- as.matrix(X[idx,])
  X0mat <- as.matrix(X0[idx,])
  fit <- cv.glmnet(Xmat, Y[idx], family = "gaussian", nfolds = 5)
  
  f0fn <- function(tau) {
    Wtau <- Y[idx] - tau*Trt[idx]
    fit_reduced <- cv.glmnet(X0mat, Wtau, family = "gaussian", nfolds = 5)
    mse_local <- mean((Wtau - predict(fit_reduced, newx=X0mat))^2)
    return(mse_local)
  }
  tau.star <- optimize(f0fn, interval=tau.range)$minimum
  
  Wtau.star <- Y[idx] - tau.star*Trt[idx]
  fit_reduced <- cv.glmnet(X0mat, Wtau.star, family = "gaussian", nfolds = 5)
  return(list(full=fit, reduced=fit_reduced, tau = tau.star))
}

fitter_glmboost <- function(X, X0, Y, Trt, tau.range, idx = NA) {
  ## X0 does not have treatment column
  if(sum(is.na(idx)) > 0) {
    idx <- 1:nrow(X)
  }
  Xmat <- as.matrix(X[idx,])
  X0mat <- as.matrix(X0[idx,])
  fit <- glmboost(x=Xmat, y=Y[idx], family = Gaussian())
  
  f0fn <- function(tau) {
    Wtau <- Y[idx] - tau*Trt[idx]
    fit_reduced <- glmboost(x=X0mat, y=Wtau,family = Gaussian())
    mse_local <- mean((Wtau - predict(fit_reduced, newdata=X0mat))^2)
    return(mse_local)
  }
  tau.star <- optimize(f0fn, interval=tau.range)$minimum
  Wtau.star <- Y[idx] - tau.star*Trt[idx]
  fit_reduced <- glmboost(x=X0mat, y=Wtau.star,family = Gaussian())
  
  return(list(full=fit, reduced=fit_reduced, tau = tau.star))
}

fitter_bart <- function(X, X0, Y, Trt, tau.range, idx = NA) {
  ## X0 does not have treatment column
  if(sum(is.na(idx)) > 0) {
    idx <- 1:nrow(X)
  }
  Xmat <- as.matrix(X[idx,])
  X0mat <- as.matrix(X0[idx,])
  fit <- bart(x.train=Xmat, y.train=Y[idx],
              ntree=50L, numcut=10L, nskip=100L, ndpost=500L, keeptrees=TRUE, printevery=1000, verbose=FALSE)
  
  f0fn <- function(tau) {
    Wtau <- Y[idx] - tau*Trt[idx]
    fit_reduced <- bart(x.train=X0mat, y.train=Wtau,
                        ntree=50L, numcut=10L, nskip=100L, ndpost=500L, keeptrees=TRUE, printevery=1000, verbose=FALSE)
    pp <- fit_reduced$yhat.train.mean
    mse_local <- mean((Wtau - pp)^2)
    return(mse_local)
  }
  tau.star <- optimize(f0fn, interval=tau.range)$minimum
  Wtau.star <- Y[idx] - tau.star*Trt[idx]
  fit_reduced <- bart(x.train=X0mat, y.train=Wtau.star,
                      ntree=50L, numcut=10L, nskip=100L, ndpost=500L, keeptrees=TRUE, printevery=1000, verbose=FALSE)
  
  return(list(full=fit, reduced=fit_reduced, tau = tau.star))
}

fitter_rf <- function(X, X0, Y, Trt, tau.range, idx = NA) {
  if(sum(is.na(idx)) > 0) {
    idx <- 1:nrow(X)
  }
  data.all <- cbind.data.frame(X[idx, ], Y = Y[idx])
  fit <- randomForest(Y ~., data = data.all, maxnodes=16, ntree=100)
  
  
  data.reduced <- cbind.data.frame(X0[idx, ], Y= Y[idx])
  f0fn <- function(tau) {
    data.reduced$Y <- Y[idx] - tau*Trt[idx]
    #data.reduced <- cbind.data.frame(X0[idx, ], Y= Y[idx] - tau*Trt[idx])
    fit_reduced <- randomForest(Y ~., data = data.reduced, maxnodes=16, ntree=100)
    mse_local <- mean((data.reduced$Y - fit_reduced$predicted)^2)
    return(mse_local)
  }
  tau.star <- optimize(f0fn, interval=tau.range)$minimum
  data.reduced <- cbind.data.frame(X0[idx, ], Y= Y[idx] - tau.star*Trt[idx])
  fit_reduced <- randomForest(Y ~., data = data.reduced, maxnodes=16, ntree=100)
  return(list(full=fit, reduced=fit_reduced, tau = tau.star))
}


# Ridge Regression
fitter_ridge <- function(X, X0, Y, Trt, tau.range, idx = NA) {
  ## X0 does not have treatment column
  if(sum(is.na(idx)) > 0) {
    idx <- 1:nrow(X)
  }
  Xmat <- as.matrix(X[idx,])
  X0mat <- as.matrix(X0[idx,])
  fit <- cv.glmnet(Xmat, Y[idx], family = "gaussian", nfolds = 5, alpha = 0)  # Ridge:  alpha = 0
  
  f0fn <- function(tau) {
    Wtau <- Y[idx] - tau*Trt[idx]
    fit_reduced <- cv.glmnet(X0mat, Wtau, family = "gaussian", nfolds = 5, alpha = 0)
    mse_local <- mean((Wtau - predict(fit_reduced, newx=X0mat))^2)
    return(mse_local)
  }
  tau.star <- optimize(f0fn, interval=tau.range)$minimum
  
  Wtau.star <- Y[idx] - tau.star*Trt[idx]
  fit_reduced <- cv.glmnet(X0mat, Wtau.star, family = "gaussian", nfolds = 5, alpha = 0)
  return(list(full=fit, reduced=fit_reduced, tau = tau.star))
}

# Bayesian Causal Forest
# Simplified and more robust BCF fitting function
fitter_bcf <- function(X, X0, Y, Trt, tau.range, idx = NA) {
  if(sum(is.na(idx)) > 0) {
    idx <- 1:nrow(X)
  }
  
  X_bcf <- as.matrix(X0[idx,])
  Y_bcf <- Y[idx]
  Trt_bcf <- Trt[idx]
    fit <- bcf(y = Y_bcf, 
               z = Trt_bcf, 
               x_control = X_bcf, 
               x_moderate = X_bcf,
               pihat = rep(0.5, length(Y_bcf)),
               nburn = 100,
               nsim = 100,
               ntree_control = 25, ## 25 at least
               ntree_moderate = 10, ## 10 at least
               n_chains = 1,
               save_tree_directory = tempdir(), 
               #save_tree_directory = NULL, 
               verbose = FALSE)
  # Fit reduced model (always use linear model for consistency)
  f0fn <- function(tau) {
    Wtau <- Y_bcf - tau * Trt_bcf
    fit_reduced <- lm(Wtau ~ X_bcf)
    mse_local <- mean((Wtau - fitted(fit_reduced))^2)
    return(mse_local)
  }
  tau.star <- optimize(f0fn, interval = tau.range)$minimum
  Wtau.star <- Y_bcf - tau.star * Trt_bcf
  fit_reduced <- lm(Wtau.star ~ X_bcf)
  
  return(list(full = fit, reduced = fit_reduced, tau = tau.star, 
              X_bcf = X_bcf))
}

# Causal Forest
fitter_causal_forest <- function(X, X0, Y, Trt, tau.range, idx = NA) {
  if(sum(is.na(idx)) > 0) {
    idx <- 1:nrow(X)
  }
  
  # For Causal Forest, X and X0 should be the same (original covariates)
  # CF handles treatment internally, so we use X0 which has the original covariates
  X_cf <- as.matrix(X0[idx,])
  Y_cf <- Y[idx] 
  Trt_cf <- Trt[idx]
  
  # Remove any NA or infinite values
  valid_idx <- complete.cases(X_cf, Y_cf, Trt_cf)
  X_cf <- X_cf[valid_idx, , drop=FALSE]
  Y_cf <- Y_cf[valid_idx]
  Trt_cf <- Trt_cf[valid_idx]
  
  # Ensure X_cf is numeric matrix
  X_cf <- apply(X_cf, 2, as.numeric)
  # Fit Causal Forest (full model with heterogeneous treatment effects)
  fit <- causal_forest(X = X_cf,
                       Y = Y_cf,
                       W = Trt_cf,
                       num.trees = 100,
                       honesty = TRUE,
                       honesty.fraction = 0.5,
                       ci.group.size = 2)
  
  # For reduced model, fit without treatment effect heterogeneity
  f0fn <- function(tau) {
    Wtau <- Y_cf - tau * Trt_cf
    # Simple regression forest for reduced model
    fit_reduced <- regression_forest(X = X_cf, Y = Wtau, num.trees = 100)
    pred_reduced <- predict(fit_reduced, X_cf)$predictions
    mse_local <- mean((Wtau - pred_reduced)^2)
    return(mse_local)
  }
  tau.star <- optimize(f0fn, interval=tau.range)$minimum
  
  Wtau.star <- Y_cf - tau.star * Trt_cf
  fit_reduced <- regression_forest(X = X_cf, Y = Wtau.star, num.trees = 100)
  
  return(list(full=fit, reduced=fit_reduced, tau = tau.star, X_cf = X_cf))
}


# bartCause fitter
fitter_bartcause <- function(X, X0, Y, Trt, tau.range, idx = NA) {
  if(sum(is.na(idx)) > 0) {
    idx <- 1:nrow(X)
  }
  
  # Combine X and X0 as confounders (assuming X includes treatment interactions)
  confounders <- as.matrix(X0[idx,])
  Y_bc <- Y[idx]
  Trt_bc <- as.numeric(Trt[idx])
  
  # Fit bartCause model (this is the "full" model)
  fit_full <- bartc(response = Y_bc,
                    treatment = Trt_bc,
                    confounders = confounders,
                    estimand = "ate",
                    method.rsp = "bart",
                    method.trt = "bart",
                    verbose = FALSE)
  
  # For reduced model, find optimal tau and fit baseline model
  f0fn <- function(tau) {
    Wtau <- Y_bc - tau * Trt_bc
    fit_reduced <- lm(Wtau ~ confounders)
    mse_local <- mean((Wtau - fitted(fit_reduced))^2)
    return(mse_local)
  }
  tau.star <- optimize(f0fn, interval=tau.range)$minimum
  
  Wtau.star <- Y_bc - tau.star * Trt_bc
  fit_reduced <- lm(Wtau.star ~ confounders)
  
  return(list(full = fit_full, reduced = fit_reduced, tau = tau.star, 
              confounders = confounders))
}

# rlearner fitter (using rlasso as default)
fitter_rlearner <- function(X, X0, Y, Trt, tau.range, idx = NA, method = "rlasso") {
  if(sum(is.na(idx)) > 0) {
    idx <- 1:nrow(X)
  }
  
  # Use X0 as features for rlearner
  features <- as.matrix(X0[idx,])
  Y_rl <- Y[idx]
  Trt_rl <- as.numeric(Trt[idx])
  
  # Calculate baseline (mean of control outcomes)
  baseline <- ifelse(sum(Trt_rl == 0) > 0, mean(Y_rl[Trt_rl == 0]), mean(Y_rl))
  
  # Fit R-learner
  if (method == "rlasso") {
    fit_full <- rlasso(features, Trt_rl, Y_rl)
  } else if (method == "rboost") {
    fit_full <- rboost(features, Trt_rl, Y_rl)
  } else if (method == "rkern") {
    fit_full <- rkern(features, Trt_rl, Y_rl)
  } else {
    stop("Unsupported rlearner method.  Choose from: rlasso, rboost, rkern")
  }
  
  # For reduced model, use glmnet to avoid ALL column name issues
  f0fn <- function(tau) {
    Wtau <- Y_rl - tau * Trt_rl
    tryCatch({
      fit_reduced <- cv.glmnet(features, Wtau, family = "gaussian", nfolds = 5, alpha = 1)
      pred <- predict(fit_reduced, newx = features, s = "lambda.min")
      mse_local <- mean((Wtau - as.vector(pred))^2)
      return(mse_local)
    }, error = function(e) {
      # Fallback:  simple mean
      return(mean(Wtau^2))
    })
  }
  
  tau.star <- optimize(f0fn, interval = tau.range)$minimum
  
  Wtau.star <- Y_rl - tau.star * Trt_rl
  fit_reduced <- cv.glmnet(features, Wtau.star, family = "gaussian", nfolds = 5, alpha = 1)
  
  return(list(full = fit_full, reduced = fit_reduced, tau = tau.star,
              features = features, baseline = baseline))
}


### Predictor Functions 
predictor_lm <- function(fit, X_new, X0_new, Trt_new) {
  pfull <- predict(fit$full, newdata = X_new)
  preduced <- predict(fit$reduced, newdata=X0_new) + fit$tau*Trt_new
  return(list(pred_full=pfull, pred_reduced=preduced))
}

predictor_glmnet <- function(fit, X_new, X0_new, Trt_new) {
  pfull <- predict(fit$full, newx=as.matrix(X_new))
  preduced <- predict(fit$reduced, newx=as.matrix(X0_new)) + fit$tau*Trt_new
  return(list(pred_full=pfull, pred_reduced=preduced))
}

predictor_glmboost <- function(fit, X_new, X0_new, Trt_new) {
  pfull <- as.numeric(predict(fit$full, newdata = as.matrix(X_new)))
  preduced <- as.numeric(predict(fit$reduced, newdata = as.matrix(X0_new))) + fit$tau*Trt_new
  return(list(pred_full=pfull, pred_reduced=preduced))
}

predictor_bart <- function(fit, X_new, X0_new, Trt_new) {
  pfull <- colMeans(predict(fit$full, newdata=as.matrix(X_new), ndpost=500L, printevery=1000, verbose=FALSE))
  preduced <- colMeans(predict(fit$reduced, newdata = as.matrix(X0_new), ndpost=500L, printevery=1000, verbose=FALSE)) + fit$tau*Trt_new
  return(list(pred_full=pfull, pred_reduced=preduced))
}

predictor_rf <- function(fit, X_new, X0_new, Trt_new) {
  pfull <- predict(fit$full, newdata = X_new)
  preduced <- predict(fit$reduced, newdata=X0_new) + fit$tau*Trt_new
  return(list(pred_full=pfull, pred_reduced=preduced))
}

# Ridge Regression Predictor
predictor_ridge <- function(fit, X_new, X0_new, Trt_new) {
  pfull <- predict(fit$full, newx=as.matrix(X_new))
  preduced <- predict(fit$reduced, newx=as.matrix(X0_new)) + fit$tau*Trt_new
  return(list(pred_full=pfull, pred_reduced=preduced))
}

# Bayesian Causal Forest Predictor
predictor_bcf <- function(fit, X_new, X0_new, Trt_new) {
  X_new_bcf <- as.matrix(X0_new)
    # BCF prediction with proper output handling
      pred_result <- predict(fit$full, 
                             x_predict_control = X_new_bcf, 
                             x_predict_moderate = X_new_bcf, 
                             z_pred = Trt_new,
                             pi_pred = rep(0.5, nrow(X_new_bcf)),
                             save_tree_directory = tempdir()
                             #save_tree_directory = NULL
                             )
      
      if(is.list(pred_result) && "mu" %in% names(pred_result) && "tau" %in% names(pred_result)) {
        mu_pred <- colMeans(pred_result$mu)  # Average across MCMC samples
        tau_pred <- colMeans(pred_result$tau)  # Average across MCMC samples
        # Full model
        pfull <- as.numeric(mu_pred + tau_pred * Trt_new)
      } else {
        pfull <- as.numeric(pred_result)
      }
      
      # Reduced model
      if(inherits(fit$reduced, "lm")) {
        preduced_base <- predict(fit$reduced, newdata = data.frame(X_new_bcf))
      } else {
        preduced_base <- rep(0, length(Trt_new))
      }
      preduced <- as.numeric(preduced_base) + fit$tau * Trt_new  
  pfull <- as.numeric(pfull)
  preduced <- as.numeric(preduced)
  
  return(list(pred_full = pfull, pred_reduced = preduced))
}

# Causal Forest Predictor  
predictor_causal_forest <- function(fit, X_new, X0_new, Trt_new) {
  X_new_cf <- as.matrix(X0_new)
  X_new_cf <- apply(X_new_cf, 2, as.numeric)
  tau_pred <- predict(fit$full, X_new_cf)$predictions
  baseline <- mean(fit$full$Y.orig)  
  pfull <- baseline + tau_pred * Trt_new
  preduced <- predict(fit$reduced, X_new_cf)$predictions + fit$tau * Trt_new
  
  return(list(pred_full=pfull, pred_reduced=preduced))
}


# bartCause predictor
# Updated bartCause predictor
predictor_bartcause <- function(fit, X_new, X0_new, Trt_new) {
  X_new_bc <- as.matrix(X0_new)
  
  # Get the posterior mean of the control outcome (mu.0) and treatment outcome (mu.1)
  mu_0 <- extract(fit$full, "mu.0")  # Control outcomes
  mu_1 <- extract(fit$full, "mu.1")  # Treatment outcomes
  
  # For new predictions, we use the posterior means
  mu_0_mean <- colMeans(mu_0)
  mu_1_mean <- colMeans(mu_1)
  
  # Full model prediction
  pfull <- ifelse(Trt_new == 1, mean(mu_1_mean), mean(mu_0_mean))
  
  # Alternatively, you can use the individual treatment effect
  # icate <- extract(fit$full, "icate")
  # baseline <- mean(mu_0_mean)
  # pfull <- baseline + colMeans(icate) * Trt_new
  
  # Reduced model prediction
  if(inherits(fit$reduced, "lm")) {
    preduced_base <- predict(fit$reduced, newdata = data.frame(fit$confounders))
    preduced <- mean(preduced_base) + fit$tau * Trt_new
  } else {
    preduced <- mean(mu_0_mean) + fit$tau * Trt_new
  }
  
  return(list(pred_full = as.numeric(pfull), pred_reduced = as.numeric(preduced)))
}
# rlearner predictor
predictor_rlearner <- function(fit, X_new, X0_new, Trt_new) {
  X_new_rl <- as.matrix(X0_new)
  
  # Get treatment effect predictions from rlearner
  tau_pred <- predict(fit$full, X_new_rl)
  
  # Use stored baseline
  baseline <- fit$baseline
  
  pfull <- baseline + tau_pred * Trt_new
  
  # Reduced model prediction using glmnet - pure matrix operations
  preduced_base <- predict(fit$reduced, newx = X_new_rl, s = "lambda.min")
  preduced <- as.vector(preduced_base) + fit$tau * Trt_new
  
  return(list(pred_full = as.numeric(pfull), pred_reduced = as.numeric(preduced)))
}


### Function Lists for New Methods ###

ridge_funs <- list(fitter = fitter_ridge,
                   predictor = predictor_ridge,
                   mse = mse,
                   loss = squared_loss,
                   name = "ridge")

bcf_funs <- list(fitter = fitter_bcf,
                 predictor = predictor_bcf,
                 mse = mse,
                 loss = squared_loss,
                 name = "bcf")

causal_forest_funs <- list(fitter = fitter_causal_forest,
                           predictor = predictor_causal_forest,
                           mse = mse,
                           loss = squared_loss,
                           name = "causal_forest")


linear_regression_funs <- list(fitter = fitter_lm,
                               predictor = predictor_lm,
                               mse = mse,
                               loss = squared_loss,
                               name = "linear_regression")

glmnet_funs <- list(fitter = fitter_glmnet,
                    predictor = predictor_glmnet,
                    mse = mse,
                    loss = squared_loss,
                    name = "glmnet")

glmboost_funs <- list(fitter = fitter_glmboost,
                      predictor = predictor_glmboost,
                      mse = mse,
                      loss = squared_loss,
                      name = "glmboost")

rf_funs <- list(fitter = fitter_rf,
                predictor = predictor_rf,
                mse = mse,
                loss = squared_loss,
                name = "rf")

bart_funs <- list(fitter = fitter_bart,
                  predictor = predictor_bart,
                  mse = mse,
                  loss = squared_loss,
                  name = "bart")

bartcause_funs <- list(fitter = fitter_bartcause,
                       predictor = predictor_bartcause,
                       mse = mse,
                       loss = squared_loss,
                       name = "bartcause")

rlearner_lasso_funs <- list(fitter = function(X, X0, Y, Trt, tau.range, idx = NA) {
  fitter_rlearner(X, X0, Y, Trt, tau.range, idx, method = "rlasso")
},
predictor = predictor_rlearner,
mse = mse,
loss = squared_loss,
name = "rlearner_lasso")

rlearner_boost_funs <- list(fitter = function(X, X0, Y, Trt, tau.range, idx = NA) {
  fitter_rlearner(X, X0, Y, Trt, tau.range, idx, method = "rboost")
},
predictor = predictor_rlearner,
mse = mse,
loss = squared_loss,
name = "rlearner_boost")

rlearner_kernel_funs <- list(fitter = function(X, X0, Y, Trt, tau.range, idx = NA) {
  fitter_rlearner(X, X0, Y, Trt, tau.range, idx, method = "rkern")
},
predictor = predictor_rlearner,
mse = mse,
loss = squared_loss,
name = "rlearner_kernel")


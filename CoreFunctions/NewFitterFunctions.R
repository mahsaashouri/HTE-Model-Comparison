# Load additional libraries for the new methods
library(glmnet)  # Already loaded, but ridge uses this too
library(bcf)     # For Bayesian Causal Forest
library(grf)     # For Causal Forest

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
               ntree_control = 10,
               ntree_moderate = 5,
               save_tree_directory = tempdir(), 
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
                       num.trees = 2000,
                       honesty = TRUE,
                       honesty.fraction = 0.5,
                       ci.group.size = 2)
  
  # For reduced model, fit without treatment effect heterogeneity
  f0fn <- function(tau) {
    Wtau <- Y_cf - tau * Trt_cf
    # Simple regression forest for reduced model
    fit_reduced <- regression_forest(X = X_cf, Y = Wtau, num.trees = 2000)
    pred_reduced <- predict(fit_reduced, X_cf)$predictions
    mse_local <- mean((Wtau - pred_reduced)^2)
    return(mse_local)
  }
  tau.star <- optimize(f0fn, interval=tau.range)$minimum
  
  Wtau.star <- Y_cf - tau.star * Trt_cf
  fit_reduced <- regression_forest(X = X_cf, Y = Wtau.star, num.trees = 2000)
  
  return(list(full=fit, reduced=fit_reduced, tau = tau.star, X_cf = X_cf))
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
                             save_tree_directory = tempdir())
      
      if(is.list(pred_result) && "mu" %in% names(pred_result) && "tau" %in% names(pred_result)) {
        mu_pred <- colMeans(pred_result$mu)  # Average across MCMC samples
        tau_pred <- colMeans(pred_result$tau)  # Average across MCMC samples
        
        # Full model: mu(x) + tau(x) * treatment
        pfull <- as.numeric(mu_pred + tau_pred * Trt_new)
      } else {
        # Fallback if unexpected format
        pfull <- as.numeric(pred_result)
      }
      
      # Reduced model:  homogeneous treatment effect
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
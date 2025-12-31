wsquared_loss <- function(thetahat_full, tau_star, w) {
  
  (thetahat_full - w)^2 - (tau_star - w)^2
}


### Fitting Functions

fitter_lm <- function(X, W, idx = NA) {
  ## X0 does not have treatment column
  if(sum(is.na(idx)) > 0) {
    idx <- 1:nrow(X)
  }
  data.all <- cbind.data.frame(X[idx, ], W = W[idx])
  fit <- lm(W ~., data = data.all)
  tau.star <- mean(W)

  return(list(full=fit, tau = tau.star))
}


fitter_glmnet <- function(X, W, idx = NA) {
  if(sum(is.na(idx)) > 0) {
    idx <- 1:nrow(X)
  }
  Xmat <- as.matrix(X[idx,])
  fit <- cv.glmnet(Xmat, W[idx], family = "gaussian", nfolds = 5)

  tau.star <- mean(W)
  return(list(full=fit, tau = tau.star))
}

fitter_ridge <- function(X, W, idx = NA) {
  if(sum(is.na(idx)) > 0) {
    idx <- 1:nrow(X)
  }
  Xmat <- as.matrix(X[idx,])
  fit <- cv.glmnet(Xmat, W[idx], family = "gaussian", nfolds = 5, alpha=0)
  
  tau.star <- mean(W)
  return(list(full=fit, tau = tau.star))
}


fitter_glmboost <- function(X, W, tau.range, idx = NA) {
  ## X0 does not have treatment column
  if(sum(is.na(idx)) > 0) {
    idx <- 1:nrow(X)
  }
  Xmat <- as.matrix(X[idx,])
  fit <- glmboost(x=Xmat, y=W[idx], family = Gaussian())

  tau.star <- mean(W)
  return(list(full=fit, tau = tau.star))
}



### Predictor Functions.

predictor_lm <- function(fit, X_new) {
  pfull <- predict(fit$full, newdata = X_new)
  return(list(pred_full=pfull, tau_star=fit$tau))
}

predictor_glmnet <- function(fit, X_new) {
  pfull <- predict(fit$full, newx = as.matrix(X_new))
  return(list(pred_full=pfull, tau_star=fit$tau))
}

predictor_ridge <- function(fit, X_new) {
  pfull <- predict(fit$full, newx = as.matrix(X_new))
  return(list(pred_full=pfull, tau_star=fit$tau))
}

predictor_glmboost <- function(fit, X_new, X0_new, Trt_new) {
  pfull <- predict(fit$full, newdata = as.matrix(X_new))
  return(list(pred_full=pfull, tau_star=fit$tau))
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


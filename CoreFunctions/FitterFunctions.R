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


### Predictor Functions.

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

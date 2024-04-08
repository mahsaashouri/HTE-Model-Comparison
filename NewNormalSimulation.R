
library(glmnet)
library(mboost)
library(BART)
library(dplyr)

# Set seed for reproducibility
set.seed(123)

# Define the number of observations
n <- 500

# Generate random values for x1 and x2 from a normal distribution
x1 <- rnorm(n, mean = 0, sd = 1)
x2 <- rnorm(n, mean = 0, sd = 1)

# Generate treatment variable A
A <- rbinom(n, size = 1, prob = 0.5)

# Define the coefficients
beta0 <- 2
beta1 <- 3
beta2 <- -1
beta3 <- 1.5
beta4 <- 0.5
beta5 <- -2

# Generate random error terms
epsilon <- rnorm(n, mean = 0, sd = 1)

Y <- beta0 + beta1 * x1 + beta2 * x2 + beta3 * A + beta4 * A * x1 + beta5 * A * x2 + epsilon

## linear

squared_loss <- function(y1, y2, y3, tau, Trt) {
  ## y1 - full model
  ## y2 - reduced model 
  ## y3 - outcome
  (y1 - y3)^2 - (y2 - (y3 - tau*Trt))^2
}

fitter_lm <- function(X, X0, Y, Trt, tau.seq, idx = NA) {
  ## X0 does not have treatment column
  if(sum(is.na(idx)) > 0) {
    idx <- 1:nrow(X)
  }
  data.all <- cbind.data.frame(X[idx, ], Y = Y[idx])
  fit <- lm(Y ~., data = data.all)
  
  mse <- rep(NA, length(tau.seq))
  for(k in 1:length(tau.seq)) {
    data.reduced <- cbind.data.frame(X0[idx, ], Y= Y[idx] - tau.seq[k]*Trt)
    lm_tmp <- lm(Y ~., data = data.reduced)
    mse[k] <- mean((data.reduced$Y - lm_tmp$fitted)^2)
  }
  tau.star <- tau.seq[which.min(mse)]
  
  data.reduced <- cbind.data.frame(X0[idx, ], Y= Y[idx] - tau.star*Trt)
  fit_reduced <- lm(Y ~., data = data.reduced)
  
  return(list(full=fit, reduced=fit_reduced, tau = tau.star))
}

predictor_lm <- function(fit, X_new) {
  predict(fit, newdata = X_new)
}

linear_regression_funs <- list(fitter = fitter_lm,
                               predictor = predictor_lm,
                               loss = squared_loss,
                               name = "linear_regression")
## glmnet

squared_loss <- function(y1, y2, y3, tau, Trt) {
  ## y1 - full model
  ## y2 - reduced model 
  ## y3 - outcome
  (y1 - y3)^2 - (y2 - (y3 - tau*Trt))^2
}

fitter_glmnet <- function(X, X0, Y, Trt, tau.seq, idx = NA) {
  ## X0 does not have treatment column
  if(sum(is.na(idx)) > 0) {
    idx <- 1:nrow(X)
  }
  fit <- cv.glmnet(X[idx, ], Y[idx], family = "gaussian", nfolds = 5)
  
  mse <- rep(NA, length(tau.seq))
  for(k in 1:length(tau.seq)) {
    glmnet_tmp <- cv.glmnet(X0[idx, ], Y[idx]-(tau.seq[k]*Treat[idx]), family = "gaussian", nfolds = 5) 
    mse[k] <- mean(((Y[idx]-(tau.seq[k]*Treat[idx])) - predict(glmnet_tmp, X0[idx, ]))^2)
  }
  tau.star <- tau.seq[which.min(mse)]
  
  
  
  fit_reduced <- cv.glmnet(X0[idx, ], Y[idx]-(tau.star*Treat[idx]), family = "gaussian", nfolds = 5) 
  
  return(list(full=fit, reduced=fit_reduced, tau = tau.star))
}

predictor_glmnet <- function(fit, X_new) {
  preds <- predict(fit, as.matrix(X_new))
  preds
}

glmnet_funs <- list(fitter = fitter_glmnet,
                    predictor = predictor_glmnet,
                    loss = squared_loss,
                    name = "glmnet")

## glmboost

squared_loss <- function(y1, y2, y3, tau, Trt) {
  ## y1 - full model
  ## y2 - reduced model 
  ## y3 - outcome
  (y1 - y3)^2 - (y2 - (y3 - tau*Trt))^2
}

fitter_glmboost <- function(X, X0, Y, Trt, tau.seq, idx = NA) {
  ## X0 does not have treatment column
  if(sum(is.na(idx)) > 0) {
    idx <- 1:nrow(X)
  }
  fit <- glmboost(x= as.matrix(X[idx, ]), y=Y[idx],family = Gaussian())
  
  mse <- rep(NA, length(tau.seq))
  for(k in 1:length(tau.seq)) {
    glmboost_tmp <- glmboost(x= as.matrix(X0[idx, ]), y=Y[idx]-(tau.seq[k]*Trt[idx]),family = Gaussian())
    mse[k] <- mean((Y[idx]-(tau.seq[k]*Trt[idx]) - glmboost_tmp$fitted())^2)
  }
  tau.star <- tau.seq[which.min(mse)]
  
  fit_reduced <- glmboost(x= as.matrix(X0[idx, ]), y=Y[idx]-(tau.star*Trt[idx]),family = Gaussian())
  
  return(list(full=fit, reduced=fit_reduced, tau = tau.star))
}

predictor_glmboost <- function(fit, X_new) {
  as.numeric(predict(fit, newdata = as.matrix(X_new)))
}

glmboost_funs <- list(fitter = fitter_glmboost,
                      predictor = predictor_glmboost,
                      loss = squared_loss,
                      name = "glmboost")

## bart

squared_loss <- function(y1, y2, y3, tau, Trt) {
  ## y1 - full model
  ## y2 - reduced model 
  ## y3 - outcome
  (y1 - y3)^2 - (y2 - (y3 - tau*Trt))^2
}

fitter_bart <- function(X, X0, Y, Trt, tau.seq, idx = NA) {
  ## X0 does not have treatment column
  if(sum(is.na(idx)) > 0) {
    idx <- 1:nrow(X)
  }
  
  fit <- wbart(X[idx, ], Y[idx], ntree=50L, numcut=10L, ndpost=50L, nskip=10L)
  
  mse <- rep(NA, length(tau.seq))
  for(k in 1:length(tau.seq)) {
    lm_tmp <- wbart(X0[idx, ], Y[idx] - tau.seq[k]*Trt, ntree=50L, numcut=10L, ndpost=50L, nskip=10L)
    mse[k] <- mean(((Y[idx] - tau.seq[k]*Trt) - lm_tmp$yhat.train)^2)
  }
  tau.star <- tau.seq[which.min(mse)]
  
  fit_reduced <- wbart(X0[idx, ], Y[idx] - tau.star*Trt, ntree=50L, numcut=10L, ndpost=50L, nskip=10L)
  
  return(list(full=fit, reduced=fit_reduced, tau = tau.star))
}

predictor_bart <- function(fit, X_new) {
  colMeans(predict(fit, newdata = X_new))
}

bart_funs <- list(fitter = fitter_bart,
                  predictor = predictor_bart,
                  loss = squared_loss,
                  name = "bart")


DATA_full <- data.frame('Y' = Y, 'x1' = x1, 'x2' = x2, 'A' = A, 'x1.t' = A*x1, 'x2.t' = A*x2)
DATA_reduced <- data.frame('Y' = Y, 'x1' = x1, 'x2' = x2, 'A' = A)

n_folds <- 10
nested_cv_reps <- 500

Y <- DATA_full$Y
DATA_full <- DATA_full[ , !(names(DATA_full) %in% c('Y'))]
DATA_full <- model.matrix(Y~.-1, data = DATA_full)


Treat <- DATA_reduced$A
DATA_reduced <- DATA_reduced[ , !(names(DATA_reduced) %in% c('Y', 'A'))]
DATA_reduced <- model.matrix(Y~.-1, data = DATA_reduced)
tau <- mean(Treat == 1) - mean(Treat == 0)
tau.range <- seq(-4*tau, 4*tau, length.out = 9)
#tau.range = seq(1,10, by =1)

## linear
nested_cv(data.frame(DATA_full), data.frame(DATA_reduced), as.vector(Y), as.vector(Treat),tau.seq = tau.range, linear_regression_funs, 
          n_folds = n_folds, reps  = nested_cv_reps, verbose = T)
## glmnet
nested_cv(data.frame(DATA_full), data.frame(DATA_reduced), as.vector(Y), as.vector(Treat), tau.seq = tau.range, glmnet_funs, 
          n_folds = n_folds, reps  = nested_cv_reps, verbose = T)
## glmboost
nested_cv(data.frame(DATA_full), data.frame(DATA_reduced), as.vector(Y), as.vector(Treat),tau.seq = tau.range, glmboost_funs, 
          n_folds = n_folds, reps  = nested_cv_reps, verbose = T)
## bart
nested_cv(data.frame(DATA_full), data.frame(DATA_reduced), as.vector(Y), as.vector(Treat),tau.seq = tau.range, bart_funs, 
          n_folds = n_folds, reps  = nested_cv_reps, verbose = T)



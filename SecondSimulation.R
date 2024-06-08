

library(glmnet)
library(mboost)
library(BART)
library(dplyr)

set.seed(123)
# Define mu and tau 
f1 <- function(x) {
  return((x[, 2]*x[, 4]*x[, 6]) + (2*x[, 2]*x[, 4]*(1-x[, 6])) + (3*x[, 2]*(1-x[, 4])*x[, 6]) + (4*x[, 2]*(1-x[, 4])*(1-x[, 4])) + (5*(1-x[, 2])*x[, 4]*x[, 6]) + 
    (6*(1-x[, 2])*x[, 4]*(1-x[, 6])) + (7*(1-x[, 2])*(1-x[, 4])*x[, 6]) + (8*(1-x[, 2])*(1-x[, 4])*(1-x[, 6])))
}

f2 <- function(x) {
  return((4*ifelse(x[, 1] > 1 & x[, 3] > 0, 1, 0)) + (4*ifelse(x[, 5] > 1 & x[, 7] > 1, 1, 0)) + (x[, 8]*x[, 9]))
}

f3 <- function(x) {
  return(0.5*(x[, 1]^2 + x[, 2] + x[, 3]^2 + x[, 4] + x[, 5]^2 + x[, 6] + x[, 7]^2 + x[, 8] + x[, 9]^2 -11))
}

f4 <- function(x) {
  (1/sqrt(2))*(f1(x) + (x[, 1] + x[, 3] + x[, 5] + x[, 7] + x[, 8] + x[, 9] - 2))
  return(x[, 7])
}

mu <- function(choice, x){
  if (choice == 1) 
    return(f1(x))
  else if (choice == 2) 
    return(f2(x))
  else if (choice == 3) 
    return(f3(x))
  #else if (choice == 4) 
  #  return(f4(x))
  else
    stop("Invalid choice for mu")
}

theta <- function(choice, x) {
  if (choice == 1) {
    return(rep(0, n))
  } else if (choice == 2) {
    return(rep(1, n))
  } else if (choice == 3) {
    return(2 + 0.1/(1 + exp(-x[, 2])))
  } else if (choice == 4) {
    return(f1(x))
  } #else if (choice == 5) {
    #return(f2(x))
  #} 
  else {
    stop("Invalid choice for tau")
  }
}

n <- 500  # Number of observations
p <- 9   # Number of features

# Generate x
x <- matrix(0, nrow = n, ncol = p)
for (j in 1:p) {
  if (j %% 2 == 0) {
    x[, j] <- rnorm(n, 0, 1)  # Normal(0, 1) if j is even
  } else {
    x[, j] <- rbinom(n, 1, 0.5)  # Bernoulli(0.5) if j is odd
  }
}
colnames(x) <-paste0("X", 1:p) 
# Generate treatment indicator A 
A <- rbinom(n, 1, 0.5)

# Generate outcome variable Y
mu_new <- mu(1, x) 
theta_new <- theta(2, x) 

Y <- numeric(n)
for (i in 1:n) {
  Y[i] <- mu_new[i] + A[i] * theta_new[i] + rnorm(1)  
}

x.t <- matrix(NA, nrow = nrow(x), ncol = ncol(x))
for(k in 1:ncol(x)){
  x.t[,k] <- x[,k]*A
}
colnames(x.t) <- paste0(colnames(x), ".t")

DATA_full <- data.frame(Y = Y, A = A, x, x.t)
DATA_reduced <- data.frame(Y = Y, A = A, x)

n_folds <- 10
nested_cv_reps <- 300

Y <- DATA_full$Y
DATA_full <- DATA_full[ , !(names(DATA_full) %in% c('Y'))]
DATA_full <- model.matrix(Y~.-1, data = DATA_full)


Treat <- DATA_reduced$A
DATA_reduced <- DATA_reduced[ , !(names(DATA_reduced) %in% c('Y', 'A'))]
DATA_reduced <- model.matrix(Y~.-1, data = DATA_reduced)
tau <- mean(Y[Treat == 1]) - mean(Y[Treat == 0])
tau.range <- seq(-4*tau, 4*tau, length.out = 20)

## linear

squared_loss <- function(y1, y2, y3, Trt, tau) {
  ## y1 - full model
  ## y2 - reduced model 
  ## y3 - outcome
  (y1 - y3)^2 - (y2 - (y3 - tau*Trt))^2
}


mse <- function(y1, y2, y3, Trt, tau) {
  ## y1 - full model
  ## y2 - reduced model 
  ## y3 - outcome
  mse_full <- mean((y1 - y3)^2)
  mse_reduced <- mean((y2 - (y3 - tau*Trt))^2)
  return(list(full = mse_full, reduced = mse_reduced))
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
                               mse = mse,
                               name = "linear_regression")
## glmnet

squared_loss <- function(y1, y2, y3, Trt, tau) {
  ## y1 - full model
  ## y2 - reduced model 
  ## y3 - outcome
  (y1 - y3)^2 - (y2 - (y3 - tau*Trt))^2
}

mse <- function(y1, y2, y3, Trt, tau) {
  ## y1 - full model
  ## y2 - reduced model 
  ## y3 - outcome
  mse_full <- mean((y1 - y3)^2)
  mse_reduced <- mean((y2 - (y3 - tau*Trt))^2)
  return(list(full = mse_full, reduced = mse_reduced))
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
                    mse = mse,
                    name = "glmnet")

## glmboost

squared_loss <- function(y1, y2, y3, Trt, tau) {
  ## y1 - full model
  ## y2 - reduced model 
  ## y3 - outcome
  (y1 - y3)^2 - (y2 - (y3 - tau*Trt))^2
}

mse <- function(y1, y2, y3, Trt, tau) {
  ## y1 - full model
  ## y2 - reduced model 
  ## y3 - outcome
  mse_full <- mean((y1 - y3)^2)
  mse_reduced <- mean((y2 - (y3 - tau*Trt))^2)
  return(list(full = mse_full, reduced = mse_reduced))
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
                      mse = mse,
                      name = "glmboost")

## bart

squared_loss <- function(y1, y2, y3, Trt, tau) {
  ## y1 - full model
  ## y2 - reduced model 
  ## y3 - outcome
  (y1 - y3)^2 - (y2 - (y3 - tau*Trt))^2
}

mse <- function(y1, y2, y3, Trt, tau) {
  ## y1 - full model
  ## y2 - reduced model 
  ## y3 - outcome
  mse_full <- mean((y1 - y3)^2)
  mse_reduced <- mean((y2 - (y3 - tau*Trt))^2)
  return(list(full = mse_full, reduced = mse_reduced))
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
                  mse = mse,
                  name = "bart")

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
#options(error = utils::dump.frames)

nested_cv(data.frame(DATA_full), data.frame(DATA_reduced), as.vector(Y), as.vector(Treat),tau.seq = tau.range, bart_funs, 
          n_folds = n_folds, reps  = nested_cv_reps, verbose = T)

#tryCatch(nested_cv(data.frame(DATA_full), data.frame(DATA_reduced), as.vector(Y), as.vector(Treat),tau.seq = tau.range, bart_funs, 
#                   n_folds = n_folds, reps  = nested_cv_reps, verbose = T), error = function(e) {
#                     cat("An error occurred:", conditionMessage(e), "\n")
#                     return(NA)  
#                   })


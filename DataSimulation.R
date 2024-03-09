# Set seed for reproducibility
set.seed(123)

# Define the number of observations
n <- 100

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
######################
## full model
######################

squared_loss <- function(y1, y2, funcs_params = NA) {
  (y1 - y2)^2
}

fitter_lm <- function(X, Y, idx = NA, funcs_params = NA) {
  if (sum(is.na(idx)) > 0) {
    idx <- 1:nrow(X)
  }
  data.all <- cbind.data.frame(X[idx, ], Y = Y[idx])
  fit <- lm(Y ~., data = data.all)  
  
  fit
}

predictor_lm <- function(fit, X_new, funcs_params = NA) {
  predict(fit, newdata = X_new)
}

linear_regression_funs <- list(fitter = fitter_lm,
                               predictor = predictor_lm,
                               loss = squared_loss,
                               name = "linear_regression")

## glmboost
library(mboost)

squared_loss <- function(y1, y2, funcs_params = NA) {
  (y1 - y2)^2
}

fitter_glmboost <- function(X, Y, idx = NA, funcs_params = NA) {
  if (sum(is.na(idx)) > 0) {
    idx <- 1:nrow(X)
  }
  fit <-  glmboost(x= as.matrix(X[idx, ]), y=Y[idx],family = Gaussian())
  fit
}

predictor_glmboost <- function(fit, X_new, funcs_params = NA) {
  as.numeric(predict(fit, newdata = as.matrix(X_new)))
}

glmboost_funs <- list(fitter = fitter_glmboost,
                      predictor = predictor_glmboost,
                      loss = squared_loss,
                      name = "glmboost")

## glmnet
library(glmnet)
misclass_loss <- function(y1, y2, funcs_params = NA) {
  (y1 - y2)^2
} 

fitter_glmnet <- function(X, Y,  idx = NA, funcs_params = NA) {
  if(sum(is.na(idx)) > 0) {idx <- 1:nrow(X)}
  fit <- glmnet(X[idx, ], Y[idx], family = "gaussian") 
  
  fit
}

predictor_glmnet <- function(fit, X_new, funcs_params = NA) {
  #preds <- predict(fit, newx = as.matrix(X_new), type = "response", s = funcs_params$best_lam)
  beta_hat <- fit$beta[, funcs_params$best_lam] 
  a0_hat <- fit$a0[funcs_params$best_lam]
  preds <- (as.matrix(X_new) %*% beta_hat + a0_hat)
  preds
} 

gaussian_lasso_funs <- list(fitter = fitter_glmnet,
                            predictor = predictor_glmnet,
                            loss = misclass_loss,
                            name = "gaussian_lasso")



DATA <- data.frame('Y' = Y, 'x1' = x1, 'x2' = x2, 'A' = A, 'x1.t' = A*x1, 'x2.t' = A*x2)
n_folds <- 10
nested_cv_reps <- 5000 

set.seed(123)
#train_idx <- sample(1:nrow(DATA), round(.7 * nrow(DATA)), replace = FALSE)
#test_idx <- setdiff(1:nrow(DATA), train_idx)

Y <- DATA$Y
DATA <- DATA[ , !(names(DATA) %in% c('Y'))]

DATA <- model.matrix(Y~.-1, data = DATA)

## if we need training and test sets
#train.set <- DATA[train_idx, ]
#Y.train <- Y[train_idx]
#test.set <- DATA[test_idx, ]
#Y.test <-  Y[test_idx]

## linear
nested_cv(data.frame(DATA), as.vector(Y), linear_regression_funs, 
          n_folds = n_folds, reps  = nested_cv_reps, verbose = T, alpha)
## glmboost
nested_cv(data.frame(DATA), as.vector(Y), glmboost_funs, 
          n_folds = n_folds, reps  = nested_cv_reps, verbose = T, alpha = 0.5)
## glmnet
fit <- cv.glmnet(as.matrix(DATA), Y, family = "gaussian")
lambdas <- fit$lambda
best_lam <- match(fit$lambda.1se, lambdas) 
lambda <- lambdas[1:best_lam]


nested_cv(data.frame(DATA), as.vector(Y), gaussian_lasso_funs, 
          n_folds = n_folds, reps  = nested_cv_reps, 
          funcs_params = list("lambdas" = lambdas, "best_lam" = best_lam), verbose = T, alpha = 0.5)
######################
## reduced model
######################

## linear
library(dplyr)

squared_loss_m <- function(y1, y2, t, tau, funcs_params = NA) {
  (y1 - (y2-(tau*t)))^2
}

fitter_lm_m <- function(X, Y, Treat, tau, idx = NA, funcs_params = NA) {
  if (sum(is.na(idx)) > 0) {
    idx <- 1:nrow(X)
  }
  data.all <- cbind.data.frame(X[idx, ], Y = Y[idx]-(tau*Treat[idx]))
  fit <- lm(Y ~., data = data.all)  
  
  fit
}

predictor_lm_m <- function(fit, X_new, funcs_params = NA) {
  predict(fit, newdata = X_new)
}

linear_regression_funs_m <- list(fitter = fitter_lm_m,
                               predictor = predictor_lm_m,
                               loss = squared_loss_m,
                               name = "linear_regression_m")

## glmboost
library(mboost)

squared_loss_m <- function(y1, y2, t, tau, funcs_params = NA) {
  (y1 - (y2-(tau*t)))^2
}

fitter_glmboost_m <- function(X, Y, Treat, tau, idx = NA, funcs_params = NA) {
  if (sum(is.na(idx)) > 0) {
    idx <- 1:nrow(X)
  }
  fit <-  glmboost(x= as.matrix(X[idx, ]), y=Y[idx]-(tau*Treat[idx]),family = Gaussian())
  fit
}

predictor_glmboost_m <- function(fit, X_new, funcs_params = NA) {
  as.numeric(predict(fit, newdata = as.matrix(X_new)))
}

glmboost_funs_m <- list(fitter = fitter_glmboost_m,
                      predictor = predictor_glmboost_m,
                      loss = squared_loss_m,
                      name = "glmboost_m")
## glmnet
library(glmnet)
misclass_loss_m <- function(y1, y2, t, tau, funcs_params = NA) {
  (y1 - (y2-(tau*t)))^2
} 

fitter_glmnet_m <- function(X, Y, Treat, tau, idx = NA, funcs_params = NA) {
  if(sum(is.na(idx)) > 0) {idx <- 1:nrow(X)}
  fit <- glmnet(X[idx, ], Y[idx]-(tau*Treat[idx]), family = "gaussian") 
  
  fit
}

predictor_glmnet_m <- function(fit, X_new, funcs_params = NA) {
  #preds <- predict(fit, newx = as.matrix(X_new), type = "response", s = funcs_params$best_lam)
  beta_hat <- fit$beta[, funcs_params$best_lam] 
  a0_hat <- fit$a0[funcs_params$best_lam]
  preds <- (as.matrix(X_new) %*% beta_hat + a0_hat)
  preds
} 

gaussian_lasso_funs_m <- list(fitter = fitter_glmnet_m,
                            predictor = predictor_glmnet_m,
                            loss = misclass_loss_m,
                            name = "gaussian_lasso_m")

DATA_m <- data.frame('Y' = Y, 'x1' = x1, 'x2' = x2, 'A' = A)

n_folds <- 10
nested_cv_reps <- 5000 #average over many random splits

set.seed(123)
train_idx_m <- sample(1:nrow(DATA_m), round(.7 * nrow(DATA_m)), replace = FALSE)
test_idx_m <- setdiff(1:nrow(DATA_m), train_idx)

Y_m <- DATA_m$Y
Treat_m <- DATA_m$A
DATA_m <- DATA_m[ , !(names(DATA_m) %in% c('Y', 'A'))]

DATA_m <- model.matrix(Y_m~.-1, data = DATA_m)

## if we need training and test sets
#train.set_m <- DATA_m[train_idx_m, ]
#Y.train_m <- Y_m[train_idx_m]
#Treat.train_m <- Treat_m[train_idx_m]
#test.set_m <- DATA_m[test_idx_m, ]
#Y.test_m <-  Y_m[test_idx_m]
#Treat.test_m <- Treat_m[test_idx_m]

tau.range = seq(1,10, by =1)
## linear
nested_cv_m(data.frame(DATA_m), as.vector(Y_m), as.vector(Treat_m), tau.range, linear_regression_funs_m, 
            n_folds = n_folds, reps  = nested_cv_reps, verbose = T, alpha)

## glmboost
nested_cv_m(data.frame(DATA_m), as.vector(Y_m), as.vector(Treat_m), tau.range , glmboost_funs_m, 
            n_folds = n_folds, reps  = nested_cv_reps, verbose = T, alpha = 0.5)

## glmnet

#Fit one model to find a good lambda. This lambda will be fixed in future simulations.
fit_m <- cv.glmnet(DATA_m, Y_m, Treat_m, family = "gaussian")
lambdas_m <- fit_m$lambda
best_lam_m <- match(fit_m$lambda.1se, lambdas_m) #selected value of lambda
lambda_m <- lambdas_m[1:best_lam_m]

nested_cv_m(data.frame(DATA_m), as.vector(Y_m), as.vector(Treat_m), tau.range, gaussian_lasso_funs_m, 
            n_folds = n_folds, reps  = nested_cv_reps, 
            funcs_params = list("lambdas" = lambdas_m, "best_lam" = best_lam_m), verbose = T, alpha = 0.5)

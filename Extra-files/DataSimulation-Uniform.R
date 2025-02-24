
# Set seed for reproducibility
set.seed(123)

# Define the number of observations
n <- 5000


# Generate random values for x1 and x2 from a normal distribution
x1 <- runif(n,0,1)
x2 <- runif(n,0,1)
x3 <- runif(n,0,1)
x4 <- runif(n,0,1)
x5 <- runif(n,0,1)


# Generate treatment variable A
A <- rbinom(n, size = 1, prob = 0.5)

# Define the coefficients
beta0 <- 2
beta1 <- 3
beta2 <- -1
beta3 <- 1.5
beta4 <- 0.5
beta5 <- -2
beta6 <- 0.3
beta7 <- -0.5
beta8 <- 4
beta9 <- 0.01
beta10 <- 1.2
beta_A <- 1.5

# Generate random error terms
epsilon <- rnorm(n, mean = 0, sd = 1)

Y <- beta0 + beta1 * x1 + beta2 * x2 + beta3 * x3 + beta4 * x4 + beta5 * x5 + beta_A * A  + beta6 * A * x1 +  beta7 * A * x2 +  
  beta8 * A * x3 +  beta9 * A * x4 +  beta10 * A*x5 + epsilon
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




DATA <- data.frame('Y' = Y, 'x1' = x1, 'x2' = x2, 'x3' = x3,'x4' = x4,'x5' = x5,'A' =A ,'x1.t' = A*x1, 'x2.t' = A*x2, 
                   'x3.t' = A*x3, 'x4.t' = A*x4, 'x5.t' = A*x5)
n_folds <- 6
nested_cv_reps <- 5000 

set.seed(123)
train_idx <- sample(1:nrow(DATA), round(.7 * nrow(DATA)), replace = FALSE)
test_idx <- setdiff(1:nrow(DATA), train_idx)

Y <- DATA$Y
DATA <- DATA[ , !(names(DATA) %in% c('Y'))]

DATA <- model.matrix(Y~.-1, data = DATA)

train.set <- DATA[train_idx, ]
Y.train <- Y[train_idx]
test.set <- DATA[test_idx, ]
Y.test <-  Y[test_idx]

## linear
nested_cv(data.frame(train.set), as.vector(Y.train), linear_regression_funs, 
          n_folds = n_folds, reps  = nested_cv_reps, verbose = T, alpha = 0.01)
## glmboost
nested_cv(data.frame(train.set), as.vector(Y.train), glmboost_funs, 
          n_folds = n_folds, reps  = nested_cv_reps, verbose = T, alpha = 0.01)
## glmnet

fit <- cv.glmnet(as.matrix(train.set), Y.train, family = "gaussian")
lambdas <- fit$lambda
best_lam <- match(fit$lambda.1se, lambdas) 
lambda <- lambdas[1:best_lam]


nested_cv(data.frame(train.set), as.vector(Y.train), gaussian_lasso_funs, 
          n_folds = n_folds, reps  = nested_cv_reps, 
          funcs_params = list("lambdas" = lambdas, "best_lam" = best_lam), verbose = T, alpha = 0.01)
######################
## reduced model
######################

## linear
library(dplyr)

squared_loss <- function(y1, y2, t, tau, funcs_params = NA) {
  (y1 - (y2-(tau*t)))^2
}

fitter_lm <- function(X, Y, Treat, tau, idx = NA, funcs_params = NA) {
  if (sum(is.na(idx)) > 0) {
    idx <- 1:nrow(X)
  }
  data.all <- cbind.data.frame(X[idx, ], Y = Y[idx]-(tau*Treat[idx]))
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

squared_loss <- function(y1, y2, t, tau, funcs_params = NA) {
  (y1 - (y2-(tau*t)))^2
}

fitter_glmboost <- function(X, Y, Treat, tau, idx = NA, funcs_params = NA) {
  if (sum(is.na(idx)) > 0) {
    idx <- 1:nrow(X)
  }
  fit <-  glmboost(x= as.matrix(X[idx, ]), y=Y[idx]-(tau*Treat[idx]),family = Gaussian())
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
misclass_loss <- function(y1, y2, t, tau, funcs_params = NA) {
  (y1 - (y2-(tau*t)))^2
} 

fitter_glmnet <- function(X, Y, Treat, tau, idx = NA, funcs_params = NA) {
  if(sum(is.na(idx)) > 0) {idx <- 1:nrow(X)}
  fit <- glmnet(X[idx, ], Y[idx]-(tau*Treat[idx]), family = "gaussian") 
  
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



DATA <- data.frame('Y' = Y, 'x1' = x1, 'x2' = x2, 'x3' = x3,'x4' = x4,'x5' = x5,'A' =A)

n_folds <- 6
nested_cv_reps <- 5000 #average over many random splits

set.seed(123)
train_idx <- sample(1:nrow(DATA), round(.7 * nrow(DATA)), replace = FALSE)
test_idx <- setdiff(1:nrow(DATA), train_idx)

Y <- DATA$Y
treat <- DATA$A
DATA <- DATA[ , !(names(DATA) %in% c('Y', 'A'))]

DATA <- model.matrix(Y~.-1, data = DATA)

train.set <- DATA[train_idx, ]
Y.train <- Y[train_idx]
Treat.train <- treat[train_idx]
test.set <- DATA[test_idx, ]
Y.test <-  Y[test_idx]
Treat.test <- treat[test_idx]

tau.range = seq(1,10, by =1)
## linear
nested_cv_m(data.frame(train.set), as.vector(Y.train), as.vector(Treat.train), tau.range, linear_regression_funs, 
            n_folds = n_folds, reps  = nested_cv_reps, verbose = T, alpha = 0.01)

## glmboost
nested_cv_m(data.frame(train.set), as.vector(Y.train), as.vector(Treat.train), tau.range , glmboost_funs, 
            n_folds = n_folds, reps  = nested_cv_reps, verbose = T, alpha = 0.01)

## glmnet

#Fit one model to find a good lambda. This lambda will be fixed in future simulations.
fit <- cv.glmnet(train.set, Y.train, Treat.train, family = "gaussian")
lambdas <- fit$lambda
best_lam <- match(fit$lambda.1se, lambdas) #selected value of lambda
lambda <- lambdas[1:best_lam]

nested_cv_m(data.frame(train.set), as.vector(Y.train), as.vector(Treat.train), tau.range, gaussian_lasso_funs, 
            n_folds = n_folds, reps  = nested_cv_reps, 
            funcs_params = list("lambdas" = lambdas, "best_lam" = best_lam), verbose = T, alpha = 0.01)


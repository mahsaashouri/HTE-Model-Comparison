
misclass_loss <- function(y1, y2, t, tau, funcs_params = NA) {
  (y1 - (y2-(tau*t)))^2
  } 

fitter_glmnet <- function(X, Y, Treat, tau, idx = NA, funcs_params = NA) {
  if(sum(is.na(idx)) > 0) {idx <- 1:nrow(X)}
  fit <- glmnet(X[idx, ], Y[idx]-(tau*Treat[idx]), family = "gaussian", lambda = funcs_params$lambdas) #assumes lambda is in global env
  
  fit
}

predictor_glmnet <- function(fit, X_new, funcs_params = NA) {
  beta_hat <- fit$beta[, funcs_params$best_lam] #assumes best_lam is in global env
  a0_hat <- fit$a0[funcs_params$best_lam]
  preds <- (as.matrix(X_new) %*% beta_hat + a0_hat > 0)
  
  preds
} 

gaussian_lasso_funs <- list(fitter = fitter_glmnet,
                            predictor = predictor_glmnet,
                            loss = misclass_loss,
                            name = "gaussian_lasso")




n_folds <- 6
nested_cv_reps <- 300 #average over many random splits


DATA <- read.csv('IHDP_clean.csv', header = TRUE)[,-1]

set.seed(123)
train_idx <- sample(1:nrow(DATA), round(.7 * nrow(DATA)), replace = FALSE)
test_idx <- setdiff(1:nrow(DATA), train_idx)

# Create training and test sets using the selected index numbers
Y <- DATA$iqsb.36
treat <- DATA$treat
DATA <- DATA[ , !(names(DATA) %in% c('iqsb.36', 'treat'))]
DATA <- model.matrix(Y~.-1, data = DATA)

train.set <- DATA[train_idx, ]
Y.train <- Y[train_idx]
Treat.train <- treat[train_idx]
test.set <- DATA[test_idx, ]
Y.test <-  Y[test_idx]
Treat.test <- treat[test_idx]
n <- nrow(train.set) #number of observations
p <- ncol(train.set) #number of features
k <- 4 #number of nonzero coefficients
alpha <- .1 #nominal error rate, total across both tails.

library(glmnet)
#Fit one model to find a good lambda. This lambda will be fixed in future simulations.
fit <- cv.glmnet(train.set, Y.train, Treat.train, family = "gaussian")
lambdas <- fit$lambda
best_lam <- match(fit$lambda.1se, lambdas) #selected value of lambda
lambda <- lambdas[1:best_lam]

tau.range = seq(1,10, by =1)
nested_cv_m(data.frame(train.set), as.vector(Y.train), as.vector(Treat.train), tau.range, gaussian_lasso_funs, 
          n_folds = n_folds, reps  = nested_cv_reps, 
          funcs_params = list("lambdas" = lambdas, "best_lam" = best_lam), verbose = T)



#######################
### Comparison 
#######################
misclass_loss <- function(y1, y2, t,  funcs_params = NA) {
  (y1 - y2)^2
} 

fitter_glmnet <- function(X, Y,  idx = NA, funcs_params = NA) {
  if(sum(is.na(idx)) > 0) {idx <- 1:nrow(X)}
  fit <- glmnet(X[idx, ], Y[idx], family = "gaussian", lambda = funcs_params$lambdas) #assumes lambda is in global env
  
  fit
}

predictor_glmnet <- function(fit, X_new, funcs_params = NA) {
  beta_hat <- fit$beta[, funcs_params$best_lam] #assumes best_lam is in global env
  a0_hat <- fit$a0[funcs_params$best_lam]
  preds <- (as.matrix(X_new) %*% beta_hat + a0_hat > 0)
  
  preds
} 

gaussian_lasso_funs <- list(fitter = fitter_glmnet,
                            predictor = predictor_glmnet,
                            loss = misclass_loss,
                            name = "gaussian_lasso")




n_folds <- 6
nested_cv_reps <- 300 #average over many random splits


DATA <- read.csv('IHDP_clean.csv', header = TRUE)[,-1]

set.seed(123)
train_idx <- sample(1:nrow(DATA), round(.7 * nrow(DATA)), replace = FALSE)
test_idx <- setdiff(1:nrow(DATA), train_idx)


Y <- DATA$iqsb.36
treat <- DATA$treat
DATA <- DATA[ , !(names(DATA) %in% c('iqsb.36', 'treat'))]
Xtrt <- matrix(NA, ncol = ncol(DATA), nrow = nrow(DATA))
for(k in 1:ncol(DATA)){
  Xtrt[,k] <- treat*DATA[,k]
}
Xtrt <- model.matrix(Y~.-1, data = as.data.frame(Xtrt))
Xfull <- cbind.data.frame(DATA, treat, Xtrt)

train.set <- Xfull[train_idx, ]
Y.train <- Y[train_idx]
Treat.train <- treat[train_idx]
test.set <- Xfull[test_idx, ]
Y.test <-  Y[test_idx]
Treat.test <- treat[test_idx]
n <- nrow(train.set) #number of observations
p <- ncol(train.set) #number of features
k <- 4 #number of nonzero coefficients
alpha <- .1 #nominal error rate, total across both tails.

library(glmnet)
#Fit one model to find a good lambda. This lambda will be fixed in future simulations.
fit <- cv.glmnet(as.matrix(train.set), Y.train, family = "gaussian")
lambdas <- fit$lambda
best_lam <- match(fit$lambda.1se, lambdas) #selected value of lambda
lambda <- lambdas[1:best_lam]


nested_cv(data.frame(train.set), as.vector(Y.train), gaussian_lasso_funs, 
          n_folds = n_folds, reps  = nested_cv_reps, 
          funcs_params = list("lambdas" = lambdas, "best_lam" = best_lam), verbose = T)
#nested_cv_helper(data.frame(train.set), as.vector(Y.train), gaussian_lasso_funs, 
#                 n_folds = 10, 
#                 funcs_params = list("lambdas" = lambdas, "best_lam" = best_lam))

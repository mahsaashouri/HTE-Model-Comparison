library(glmnet)
library(mboost)
library(dbarts)

##################################
## Simulation Study 2

# Set seed for reproducibility
set.seed(12356)

# Define mu and theta 
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
  else if (choice == 4) 
    return(f4(x))
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
  } else if (choice == 5) {
  return(f2(x))
  } 
  else {
    stop("Invalid choice for tau")
  }
}


## Source fitting and nested cv functions
setwd("/Users/mahsa/Projects/HTE-Model-Comparison")  ## Change for your computer
source("CoreFunctions/CVDiffFunctions.R")
source("CoreFunctions/FitterFunctions.R")

## Define linear regression, glmnet, and boosting fitting and prediction functions.
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

bart_funs <- list(fitter = fitter_bart,
                  predictor = predictor_bart,
                  mse = mse,
                  loss = squared_loss,
                  name = "bart")

# Set the number of observations n, number of folds, and
# number of nested cv replications:
n <- 500
p <- 9   # Number of features
n_folds <- 5
nested_cv_reps <- 50 ## Use 50 or 100 for paper

## Set the number of simulation replications
nreps <- 500  ## Use nreps = 500 for paper
cover_lm <- cover_glmnet <- cover_glmboost <- cover_bart <- rep(NA, nreps)
hvalue_lm <- hvalue_glmnet <- hvalue_glmboost <- hvalue_bart <- rep(NA, nreps)
CI_lm <- CI_glmnet <- CI_glmboost <- CI_bart <- matrix(NA, nrow=nreps, ncol=2)
true_thetas <- matrix(NA, nreps, 4)
for(h in 1:nreps) {
  # Generate random values for x1 and x2 from a normal distribution
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
  A <- rbinom(n, size = 1, prob = 0.5)
  
  # Generate outcome variable Y
  mu_new <- mu(1, x) 
  theta_new <- theta(5, x) 
  
  # Generate residuals and outcomes
  epsilon <- rnorm(n, mean = 0, sd = 1)
  
  Y <- numeric(n)
  for (i in 1:n) {
    Y[i] <- mu_new[i] + A[i] * theta_new[i] + epsilon[i]
  }
  
  x.t <- matrix(NA, nrow = nrow(x), ncol = ncol(x))
  for(k in 1:ncol(x)){
    x.t[,k] <- x[,k]*A
  }
  colnames(x.t) <- paste0(colnames(x), ".t")
  ######################
  ## Compute true value of \theta_{XY} for each method
  ######################
  
  ## First fit each model with the dataset of size n
  DAT <- data.frame('Y' = Y, 'A' = A, x, x.t)
  DAT_reduced <- data.frame('Y' = Y, 'A' = A, x)
  
  Xmat_tmp <- model.matrix(Y~.-1, data = DAT)
  X0mat_tmp <- model.matrix(Y~.-1, data = DAT_reduced)
  ## Don't use treatment variable in GLMboost
  X0mat_notrt_tmp <- model.matrix(Y ~ . - A - 1, data = DAT_reduced)

  tmp_lm <- lm(Y ~., data=DAT)
  tmp_reduced_lm <- lm(Y ~., data=DAT_reduced)
  tmp_glmnet <- cv.glmnet(Xmat_tmp, DAT$Y, family = "gaussian", nfolds = 5)
  tmp_reduced_glmnet <- cv.glmnet(X0mat_tmp, DAT$Y, family = "gaussian", nfolds = 5)
  
  tmp_glmboost <- glmboost(x=Xmat_tmp, y=DAT$Y, family = Gaussian())
  
  f0fn <- function(tau) {
    Wtau <- DAT_reduced$Y - tau*DAT_reduced$A
    fit_reduced <- glmboost(x=X0mat_notrt_tmp, y=Wtau,family = Gaussian())
    mse_local <- mean((Wtau - predict(fit_reduced, newdata=X0mat_notrt_tmp))^2)
    return(mse_local)
  }
  tau.star.gboost <- optimize(f0fn, interval=c(-5, 5))$minimum
  Wtau.star <- DAT_reduced$Y - tau.star.gboost*DAT_reduced$A
  tmp_reduced_glmboost <- glmboost(x=X0mat_notrt_tmp, y = Wtau.star,family = Gaussian())
  
  tmp_bart <- bart(x.train=Xmat_tmp, y.train=DAT$Y,
                   ntree=50L, numcut=10L, nskip=100L, ndpost=500L, keeptrees=TRUE, verbose=FALSE)
  f00fn <- function(tau) {
    Wtau <- DAT_reduced$Y - tau*DAT_reduced$A
    fit_reduced <- bart(x.train=X0mat_notrt_tmp, y.train=Wtau, ntree=50L, numcut=10L,
                        nskip=100L, ndpost=500L, keeptrees=TRUE, verbose=FALSE)
    mse_local <- mean((Wtau - predict(fit_reduced, newdata=X0mat_notrt_tmp))^2)
    return(mse_local)
  }
  tau.star.bart <- optimize(f00fn, interval=c(-5, 5))$minimum
  Wtau.star <- DAT_reduced$Y - tau.star.bart*DAT_reduced$A
  tmp_reduced_bart <- bart(x.train=X0mat_notrt_tmp, y.train = Wtau.star, ntree=50L, numcut=10L, nskip=100L,
                           ndpost=500L, keeptrees=TRUE, verbose=FALSE)
  
  ## Now, evaluate MSE difference on a much larger "future" dataset
  nr <- 100000
  
  ## Generating future observations:
  epsilon <- rnorm(nr, mean = 0, sd = 1)
  xk <- matrix(0, nrow = nr, ncol = p)
  for (j in 1:p) {
    if (j %% 2 == 0) {
      xk[, j] <- rnorm(nr, 0, 1)  # Normal(0, 1) if j is even
    } else {
      xk[, j] <- rbinom(nr, 1, 0.5)  # Bernoulli(0.5) if j is odd
    }
  }
  colnames(xk) <-paste0("X", 1:p) 
  Ak <- rbinom(nr, size = 1, prob = 0.5)
  
  Yk <- numeric(nr)
  for (i in 1:n) {
    Yk[i] <- mu_new[i] + Ak[i] * theta_new[i] + epsilon[i] 
  }
  
  x.tk <- matrix(NA, nrow = nrow(xk), ncol = ncol(xk))
  for(k in 1:ncol(xk)){
    x.tk[,k] <- xk[,k]*Ak
  }
  colnames(x.tk) <- paste0(colnames(xk), ".t")
 
  
  DATk <- data.frame('Y' = Yk, 'A' = Ak, xk, x.tk)
  DATk_reduced <- data.frame('Y' = Yk, 'A' = A, xk)
  Xmat_tmpk <- model.matrix(Y~.-1, data = DATk)
  X0mat_tmpk <- model.matrix(Y~.-1, data = DATk_reduced)
  X0mat_notrt_tmpk <- model.matrix(Y ~ . - A - 1, data = DATk_reduced)
  
  ## Getting predictions for these future observations
  pp_reduced_lm <- predict(tmp_reduced_lm, newdata=DATk_reduced)
  pp_full_lm <- predict(tmp_lm, newdata=DATk)
  pp_reduced_glmnet <- as.numeric(predict(tmp_reduced_glmnet, newx=X0mat_tmpk))
  pp_full_glmnet <- as.numeric(predict(tmp_glmnet, newx=Xmat_tmpk))
  pp_reduced_glmboost <- as.numeric(predict(tmp_reduced_glmboost, newdata=X0mat_notrt_tmpk)) + tau.star.gboost*DATk_reduced$A
  pp_full_glmboost <- as.numeric(predict(tmp_glmboost, newdata = Xmat_tmpk))
  pp_reduced_bart <- colMeans(predict(tmp_reduced_bart, newdata=X0mat_notrt_tmpk, ndpost=500L, verbose=FALSE)) + tau.star.bart*DATk_reduced$A
  pp_full_bart <- colMeans(predict(tmp_bart, newdata=Xmat_tmpk, ndpost=500L, verbose=FALSE))
  
  ## Compute \theta_{XY} by looking at differences in MSE
  theta_lm <- mean((Yk - pp_reduced_lm)^2) - mean((Yk - pp_full_lm)^2)
  theta_glmnet <- mean((Yk - pp_reduced_glmnet)^2) - mean((Yk - pp_full_glmnet)^2)
  theta_glmboost <- mean((Yk - pp_reduced_glmboost)^2) - mean((Yk - pp_full_glmboost)^2)
  theta_bart <- mean((Yk - pp_reduced_bart)^2) - mean((Yk - pp_full_bart)^2)
  
  #######################################################
  ###. Getting confidence intervals and h-values
  ######################################################
  naive.trt.effect <- mean(Y[A==1]) - mean(Y[A==0])
  tau.range <- sort(c(-5*naive.trt.effect, 5*naive.trt.effect))
  XX <- data.frame(model.matrix(Y ~. - 1, data = DAT))
  XX0 <- data.frame(model.matrix(Y ~. - A - 1, data = DAT_reduced))
  ncv_lm <- nested_cv(X=XX, X0=XX0, Y=as.vector(Y), Trt=A, tau.range=tau.range, funcs=linear_regression_funs,
                      n_folds = n_folds, reps  = nested_cv_reps)
  
  ncv_net <- nested_cv(X=XX, X0=XX0, Y=as.vector(Y), Trt=A, tau.range=tau.range, funcs=glmnet_funs,
                       n_folds = n_folds, reps  = nested_cv_reps)
  
  ncv_boost <- nested_cv(X=XX, X0=XX0, Y=as.vector(Y), Trt=A, tau.range=tau.range, funcs=glmboost_funs,
                         n_folds = n_folds, reps  = nested_cv_reps)
  
  #ncv_bart <- nested_cv(X=XX, X0=XX0, Y=as.vector(Y), Trt=A, tau.range=tau.range, funcs=bart_funs,
  #                       n_folds = n_folds, reps  = nested_cv_reps)
  ############################################
  ### Record Results
  ############################################
  cover_lm[h] <- theta_lm > ncv_lm$ci_lo & theta_lm < ncv_lm$ci_hi
  CI_lm[h,1] <- ncv_lm$ci_lo
  CI_lm[h,2] <- ncv_lm$ci_hi
  true_thetas[h,] <- c(theta_lm, theta_glmnet, theta_glmboost, theta_bart)
  hvalue_lm[h] <- ncv_lm$hvalue
  
  cover_glmnet[h] <- theta_glmnet > ncv_net$ci_lo & theta_glmnet < ncv_net$ci_hi
  CI_glmnet[h,1] <- ncv_net$ci_lo
  CI_glmnet[h,2] <- ncv_net$ci_hi
  hvalue_glmnet[h] <- ncv_net$hvalue
  
  cover_glmboost[h] <- theta_glmboost > ncv_boost$ci_lo & theta_glmboost < ncv_boost$ci_hi
  CI_glmboost[h,1] <- ncv_boost$ci_lo
  CI_glmboost[h,2] <- ncv_boost$ci_hi
  hvalue_glmboost[h] <- ncv_boost$hvalue
  
  #cover_bart[h] <- theta_bart > ncv_bart$ci_lo & theta_bart < ncv_bart$ci_hi
  #CI_bart[h,1] <- ncv_bart$ci_lo
  #CI_bart[h,2] <- ncv_bart$ci_hi
  #hvalue_bart[h] <- ncv_bart$hvalue
  cat("Simulation Replication: ", h, "\n")
}
## Save: cover_lm, cover_glmnet, cover_glmboost, cover_bart
##       CI_lm, CI_glmnet, CI_glmboost, CI_bart
##      hvalue_lm, hvalue_glmnet, hvalue_glmboost, hvalue_bart


library(glmnet)
library(mboost)
library(dbarts)

##################################
## Simulation Study 1
##      Data generated from a linear model with two covariates,
##      in addition to treatment variable

# Set seed for reproducibility
set.seed(12356)

# Define the coefficients
beta0 <- 2
beta1 <- 3
beta2 <- -1
beta3 <- 1.5
beta4 <- 0.5
beta5 <- -2

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
  x1 <- rnorm(n, mean = 0, sd = 1)
  x2 <- rnorm(n, mean = 0, sd = 1)

  # Generate treatment variable A
  A <- rbinom(n, size = 1, prob = 0.5)

  # Generate residuals and outcomes
  epsilon <- rnorm(n, mean = 0, sd = 1)

  Y <- beta0 + beta1*x1 + beta2*x2 + beta3*A + beta4*A*x1 + beta5*A* x2 + epsilon

  ######################
  ## Compute true value of \theta_{XY} for each method
  ######################

  ## First fit each model with the dataset of size n
  DAT <- data.frame('Y' = Y, 'x1' = x1, 'x2' = x2, 'A' = A, 'x1.t' = A*x1, 'x2.t' = A*x2)
  DAT_reduced <- data.frame('Y' = Y, 'x1' = x1, 'x2' = x2, 'A' = A)

  Xmat_tmp <- model.matrix(Y ~ x1 + x2 + A + A:x1 + A:x2 - 1, data=DAT)
  X0mat_tmp <- model.matrix(Y ~ x1 + x2 + A - 1, data=DAT_reduced)
  ## Don't use treatment variable in GLMboost
  X0mat_notrt_tmp <- model.matrix(Y ~ x1 + x2 - 1, data=DAT_reduced)

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
  x1k <- rnorm(nr, mean = 0, sd = 1)
  x2k <- rnorm(nr, mean = 0, sd = 1)
  Ak <- rbinom(nr, size = 1, prob = 0.5)
  Yk <- beta0 + beta1*x1k + beta2*x2k + beta3*Ak + beta4*Ak*x1k + beta5*Ak*x2k + epsilon

  DATk <- data.frame('Y' = Yk, 'x1' = x1k, 'x2' = x2k, 'A' = Ak, 'x1.t' = Ak*x1k, 'x2.t' = Ak*x2k)
  DATk_reduced <- data.frame('Y' = Yk, 'x1' = x1k, 'x2' = x2k, 'A' = Ak)
  Xmat_tmpk <- model.matrix(Y ~ x1 + x2 + A + A:x1 + A:x2 - 1, data=DATk)
  X0mat_tmpk <- model.matrix(Y ~ x1 + x2 + A - 1, data=DATk_reduced)
  X0mat_notrt_tmpk <- model.matrix(Y ~ x1 + x2 - 1, data=DATk_reduced)

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
  XX <- data.frame(model.matrix(Y ~ x1 + x2 + A + A:x1 + A:x2 - 1))
  XX0 <- data.frame(model.matrix(Y ~ x1 + x2 - 1))
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


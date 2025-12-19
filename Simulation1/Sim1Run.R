library(glmnet)
library(mboost)
library(dbarts)
library(randomForest)
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
#beta4 <- 0
#beta5 <- 0

## Source fitting and nested cv functions
setwd("~/Projects/HTE-Model-Comparison")  ## Change for your computer
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

rf_funs <- list(fitter = fitter_rf,
                predictor = predictor_rf,
                mse = mse,
                loss = squared_loss,
                name = "rf")

#bart_funs <- list(fitter = fitter_bart,
#                      predictor = predictor_bart,
#                      mse = mse,
#                      loss = squared_loss,
#                      name = "bart")

# Set the number of observations n, number of folds, and
# number of nested cv replications:
n <- 1000
n_folds <- 5
nested_cv_reps <- 50 ## Use 50 or 100 for paper

## Set the number of simulation replications
nreps <- 500  ## Use nreps = 500 for paper
cover_lm <- cover_glmnet <- cover_glmboost <- cover_rf <- rep(NA, nreps)
hvalue_lm <- hvalue_glmnet <- hvalue_glmboost <- hvalue_rf <- rep(NA, nreps)
CI_lm <- CI_glmnet <- CI_glmboost <- CI_rf <- matrix(NA, nrow=nreps, ncol=2)
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
  drop_cols <- c(which(colnames(DAT_reduced)=="Y"), which(colnames(DAT_reduced)=="A"))
  DAT_red <- data.frame(Wtau = DAT_reduced$Y, x1 = x1, x2 = x2)

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
    DAT_red$Wtau <- DAT_reduced$Y - tau*DAT_reduced$A
    
    fit_reduced <- randomForest(Wtau ~., data=DAT_red, maxnodes=16, ntree=100)
    mse_local <- mean((DAT_red$Wtau - predict(fit_reduced, newdata=DAT_red))^2)
    return(mse_local)
  }
  tau.star.rf <- optimize(f00fn, interval=c(-5, 5))$minimum
  DAT_red$Wtau <- DAT_reduced$Y - tau.star.rf*DAT_reduced$A
  ## Is the line below correct?
  tmp_reduced_rf <- randomForest(Wtau ~., data=DAT_red, maxnodes=16, ntree=100)
  
  ## BART
  #f00fn <- function(tau) {
  #  Wtau <- DAT_reduced$Y - tau*DAT_reduced$A
  #  fit_reduced <- bart(x.train=X0mat_notrt_tmp, y.train=Wtau, ntree=50L, numcut=10L,
  #                      nskip=100L, ndpost=500L, keeptrees=TRUE, verbose=FALSE)
  #  mse_local <- mean((Wtau - predict(fit_reduced, newdata=X0mat_notrt_tmp))^2)
  #  return(mse_local)
  #}
  #tau.star.bart <- optimize(f00fn, interval=c(-5, 5))$minimum
  #Wtau.star <- DAT_reduced$Y - tau.star.bart*DAT_reduced$A
  #tmp_reduced_bart <- bart(x.train=X0mat_notrt_tmp, y.train = Wtau.star, ntree=50L, numcut=10L, nskip=100L,
  #                         ndpost=500L, keeptrees=TRUE, verbose=FALSE)

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
  DATk_red <- data.frame('Y' = Yk, 'x1' = x1k, 'x2' = x2k)
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
  #pp_reduced_rf <- predict(tmp_reduced_rf, newdata=DATk_red) + tau.star.rf*DATk_reduced$A
  #pp_full_rf <- predict(tmp_rf, newdata=DATk_reduced)
  new_rf <- fitter_rf(X0mat_tmp, X0mat_notrt_tmp, Y, A, tau.range=c(-5,5), idx = NA)
  new_rf_pred <- predictor_rf(new_rf, X_new=DATk_reduced, X0_new=DATk_red, Trt_new=DATk_reduced$A)
  ## BART
  #pp_reduced_bart <- colMeans(predict(tmp_reduced_bart, newdata=X0mat_notrt_tmpk, ndpost=500L, verbose=FALSE)) + tau.star.bart*DATk_reduced$A
  #pp_full_bart <- colMeans(predict(tmp_bart, newdata=Xmat_tmpk, ndpost=500L, verbose=FALSE))
  
  #A <- cbind(hh$pred_full, pp_full_rf, hh$pred_reduced, pp_reduced_rf)
  #X0mat_tmp <- model.matrix(Y~.-1, data = DAT_reduced)
  ## Don't use treatment variable in GLMboost
  #X0mat_notrt_tmp <- model.matrix(Y ~ . - A - 1, data = DAT_reduced)
  
  theta_rf <- mean(squared_loss(new_rf_pred$pred_full, new_rf_pred$pred_reduced, Yk))

  ## Compute \theta_{XY} by looking at differences in MSE
  theta_lm <- mean((Yk - pp_reduced_lm)^2) - mean((Yk - pp_full_lm)^2)
  theta_glmnet <- mean((Yk - pp_reduced_glmnet)^2) - mean((Yk - pp_full_glmnet)^2)
  theta_glmboost <- mean((Yk - pp_reduced_glmboost)^2) - mean((Yk - pp_full_glmboost)^2)
  #theta_bart <- mean((Yk - pp_reduced_bart)^2) - mean((Yk - pp_full_bart)^2)

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
  
  ncv_rf <- nested_cv(X=XX, X0=XX0, Y=as.vector(Y), Trt=A, tau.range=tau.range, funcs=rf_funs,
                      n_folds = n_folds, reps = 10, bias_reps = 0)

  #ncv_bart <- nested_cv(X=XX, X0=XX0, Y=as.vector(Y), Trt=A, tau.range=tau.range, funcs=bart_funs,
  #                       n_folds = n_folds, reps  = nested_cv_reps)
  ############################################
  ### Record Results
  ############################################
  cover_lm[h] <- theta_lm > ncv_lm$ci_lo & theta_lm < ncv_lm$ci_hi
  CI_lm[h,1] <- ncv_lm$ci_lo
  CI_lm[h,2] <- ncv_lm$ci_hi
  true_thetas[h,] <- c(theta_lm, theta_glmnet, theta_glmboost, theta_rf)
  hvalue_lm[h] <- ncv_lm$hvalue

  cover_glmnet[h] <- theta_glmnet > ncv_net$ci_lo & theta_glmnet < ncv_net$ci_hi
  CI_glmnet[h,1] <- ncv_net$ci_lo
  CI_glmnet[h,2] <- ncv_net$ci_hi
  hvalue_glmnet[h] <- ncv_net$hvalue

  cover_glmboost[h] <- theta_glmboost > ncv_boost$ci_lo & theta_glmboost < ncv_boost$ci_hi
  CI_glmboost[h,1] <- ncv_boost$ci_lo
  CI_glmboost[h,2] <- ncv_boost$ci_hi
  hvalue_glmboost[h] <- ncv_boost$hvalue
  
  cover_rf[h] <- theta_rf > ncv_rf$ci_lo & theta_rf < ncv_rf$ci_hi
  CI_rf[h,1] <- ncv_rf$ci_lo
  CI_rf[h,2] <- ncv_rf$ci_hi
  hvalue_rf[h] <- ncv_rf$hvalue

  #cover_bart[h] <- theta_bart > ncv_bart$ci_lo & theta_bart < ncv_bart$ci_hi
  #CI_bart[h,1] <- ncv_bart$ci_lo
  #CI_bart[h,2] <- ncv_bart$ci_hi
  #hvalue_bart[h] <- ncv_bart$hvalue
  cat("Simulation Replication: ", h, "\n")
}
## Save: cover_lm, cover_glmnet, cover_glmboost, cover_bart
##       CI_lm, CI_glmnet, CI_glmboost, CI_bart
##      hvalue_lm, hvalue_glmnet, hvalue_glmboost, hvalue_bart

save(
  cover_lm, cover_glmnet, cover_glmboost, cover_rf,
  CI_lm, CI_glmnet, CI_glmboost, 
  hvalue_lm, hvalue_glmnet, hvalue_glmboost, hvalue_rf,
  true_thetas,
  file = "Sim1B_1000.RData"
)

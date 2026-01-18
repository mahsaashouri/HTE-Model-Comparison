library(glmnet)
library(mboost)
library(dbarts)
library(randomForest)
library(bcf)     # For Bayesian Causal Forest
library(grf)     # For Causal Forest

##################################
## Simulation Study 1 - Extended

# Set seed for reproducibility
set.seed(12356)

# Define the coefficients
beta0 <- 2
beta1 <- 3
beta2 <- -1
beta3 <- 1.5
#beta4 <- 0.5
#beta5 <- -2
beta4 <- 0
beta5 <- 0

## Source fitting and nested cv functions
setwd("~/Projects/HTE-Model-Comparison")  ## Change for your computer
source("CoreFunctions/CVDiffFunctions.R")
#source("CoreFunctions/FitterFunctions.R")
source("CoreFunctions/NewFitterFunctions.R")  # Contains new methods

## Define all method function lists (original + new)
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

# New method function lists
#ridge_funs <- list(fitter = fitter_ridge,
#                   predictor = predictor_ridge,
#                   mse = mse,
#                   loss = squared_loss,
#                   name = "ridge")

#bcf_funs <- list(fitter = fitter_bcf,
#                 predictor = predictor_bcf,
#                 mse = mse,
#                 loss = squared_loss,
#                 name = "bcf")

#causal_forest_funs <- list(fitter = fitter_causal_forest,
#                           predictor = predictor_causal_forest,
#                           mse = mse,
#                           loss = squared_loss,
#                           name = "causal_forest")

# Set the number of observations n, number of folds, and
# number of nested cv replications:  
n <- 100
n_folds <- 5
nested_cv_reps <- 50 ## Use 50 or 100 for paper

## Set the number of simulation replications
nreps <- 500  ## Use nreps = 500 for paper

# Initialize storage for all methods (original + new)
cover_lm <- cover_glmnet <- cover_glmboost <- cover_rf <- rep(NA, nreps)
#cover_ridge <- cover_bcf <- cover_causal_forest <- rep(NA, nreps)

hvalue_lm <- hvalue_glmnet <- hvalue_glmboost <- hvalue_rf <- rep(NA, nreps)
#hvalue_ridge <- hvalue_bcf <- hvalue_causal_forest <- rep(NA, nreps)

# ADD:  Standard error storage (proper extraction from nested_cv)
se_lm <- se_glmnet <- se_glmboost <- se_rf <- rep(NA, nreps)
#se_ridge <- se_bcf <- se_causal_forest <- rep(NA, nreps)

# ADD: Additional nested CV outputs for more detailed analysis
sd_infl_lm <- sd_infl_glmnet <- sd_infl_glmboost <- sd_infl_rf <- rep(NA, nreps)
#sd_infl_ridge <- sd_infl_bcf <- sd_infl_causal_forest <- rep(NA, nreps)

bias_est_lm <- bias_est_glmnet <- bias_est_glmboost <- bias_est_rf <- rep(NA, nreps)
#bias_est_ridge <- bias_est_bcf <- bias_est_causal_forest <- rep(NA, nreps)

CI_lm <- CI_glmnet <- CI_glmboost <- CI_rf <- matrix(NA, nrow=nreps, ncol=2)
#CI_ridge <- CI_bcf <- CI_causal_forest <- matrix(NA, nrow=nreps, ncol=2)

#true_thetas <- matrix(NA, nreps, 7)  # Now 7 methods
true_thetas <- matrix(NA, nreps, 4)  # Now 7 methods


for(h in 1:nreps) {
  # Generate random values for x1 and x2 from a normal distribution
  x1 <- rnorm(n, mean = 0, sd = 1)
  x2 <- rnorm(n, mean = 0, sd = 1)
  
  # Generate treatment variable A
  A <- rbinom(n, size = 1, prob = 0.5)
  
  # Generate residuals and outcomes
  epsilon <- rnorm(n, mean = 0, sd = 1)
  
  Y <- beta0 + beta1*x1 + beta2*x2 + beta3*A + beta4*A*x1 + beta5*A*x2 + epsilon
  
  ######################
  ## Compute true value of \theta_{XY} for each method
  ######################
  
  ## First fit each model with the dataset of size n
  DAT <- data.frame('Y' = Y, 'x1' = x1, 'x2' = x2, 'A' = A, 'x1.t' = A*x1, 'x2.t' = A*x2)
  DAT_reduced <- data.frame('Y' = Y, 'x1' = x1, 'x2' = x2, 'A' = A)
  drop_cols <- c(which(colnames(DAT_reduced)=="Y"), which(colnames(DAT_reduced)=="A"))
  DAT_red <- data.frame(Wtau = DAT_reduced$Y, x1 = x1, x2 = x2)
  
  Xmat_tmp <- model.matrix(Y ~ x1 + x2 + A + A: x1 + A:x2 - 1, data=DAT)
  X0mat_tmp <- model.matrix(Y ~ x1 + x2 + A - 1, data=DAT_reduced)
  ## Don't use treatment variable in GLMboost
  X0mat_notrt_tmp <- model.matrix(Y ~ x1 + x2 - 1, data=DAT_reduced)
  
  # Original method fits
  tmp_lm <- lm(Y ~., data=DAT)
  tmp_reduced_lm <- lm(Y ~., data=DAT_reduced)
  tmp_glmnet <- cv.glmnet(Xmat_tmp, DAT$Y, family = "gaussian", nfolds = 5, alpha = 1)  # LASSO
  tmp_reduced_glmnet <- cv.glmnet(X0mat_tmp, DAT$Y, family = "gaussian", nfolds = 5, alpha = 1)
  
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
  tmp_reduced_rf <- randomForest(Wtau ~., data=DAT_red, maxnodes=16, ntree=100)
  
  # New method fits
  # Ridge
  #tmp_ridge <- cv.glmnet(Xmat_tmp, DAT$Y, family = "gaussian", nfolds = 5, alpha = 0)  # Ridge
  #tmp_reduced_ridge <- cv.glmnet(X0mat_tmp, DAT$Y, family = "gaussian", nfolds = 5, alpha = 0)
  
  # BCF (with error handling)
  
  #bcf_result <- fitter_bcf(data.frame(Xmat_tmp), data.frame(X0mat_notrt_tmp), Y, A, tau.range=c(-5,5), idx = NA)
  #tmp_bcf <- bcf_result$full
  #tmp_reduced_bcf <- bcf_result$reduced
  #tau.star.bcf <- bcf_result$tau
  
  
  # Causal Forest
  #tmp_cf <- causal_forest(X = cbind(x1, x2), Y = Y, W = A, num.trees = 1000)
  #f_cf <- function(tau) {
  #  Wtau <- Y - tau * A
  #  fit_reduced <- regression_forest(X = cbind(x1, x2), Y = Wtau, num.trees = 1000)
  #  pred_reduced <- predict(fit_reduced, cbind(x1, x2))$predictions
  #  mse_local <- mean((Wtau - pred_reduced)^2)
  #  return(mse_local)
  #}
  #tau.star.cf <- optimize(f_cf, interval=c(-5, 5))$minimum
  #Wtau.star.cf <- Y - tau.star.cf * A
  #tmp_reduced_cf <- regression_forest(X = cbind(x1, x2), Y = Wtau.star.cf, num.trees = 1000)
  
  
  
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
  
  ## Getting predictions for these future observations (original methods)
  pp_reduced_lm <- predict(tmp_reduced_lm, newdata=DATk_reduced)
  pp_full_lm <- predict(tmp_lm, newdata=DATk)
  pp_reduced_glmnet <- as.numeric(predict(tmp_reduced_glmnet, newx=X0mat_tmpk))
  pp_full_glmnet <- as.numeric(predict(tmp_glmnet, newx=Xmat_tmpk))
  pp_reduced_glmboost <- as.numeric(predict(tmp_reduced_glmboost, newdata=X0mat_notrt_tmpk)) + tau.star.gboost*DATk_reduced$A
  pp_full_glmboost <- as.numeric(predict(tmp_glmboost, newdata = Xmat_tmpk))
  
  new_rf <- fitter_rf(X0mat_tmp, X0mat_notrt_tmp, Y, A, tau.range=c(-5,5), idx = NA)
  new_rf_pred <- predictor_rf(new_rf, X_new=DATk_reduced, X0_new=DATk_red, Trt_new=DATk_reduced$A)
  
  # New method predictions
  # Ridge
  #pp_reduced_ridge <- as.numeric(predict(tmp_reduced_ridge, newx=X0mat_tmpk))
  #pp_full_ridge <- as.numeric(predict(tmp_ridge, newx=Xmat_tmpk))
  
  # BCF
  #bcf_pred <- predictor_bcf(bcf_result, X_new=DATk_reduced, X0_new=DATk_red, Trt_new=DATk_reduced$A)
  #pp_full_bcf <-  as.numeric(bcf_pred$pred_full)
  #pp_reduced_bcf <- as.numeric(bcf_pred$pred_reduced)
  
  # Causal Forest
  #mu_forest <- regression_forest(cbind(x1, x2), Y)
  #mu_pred_k <- predict(mu_forest, cbind(x1k, x2k))$predictions
  #tau_pred_k <- predict(tmp_cf, cbind(x1k, x2k))$predictions
  #pp_full_cf <- mu_pred_k + tau_pred_k * Ak
  #pp_reduced_cf <- predict(tmp_reduced_cf, cbind(x1k, x2k))$predictions + tau.star.cf * Ak
  
  ## Compute \theta_{XY} by looking at differences in MSE (original methods)
  theta_lm <- mean((Yk - pp_reduced_lm)^2) - mean((Yk - pp_full_lm)^2)
  theta_glmnet <- mean((Yk - pp_reduced_glmnet)^2) - mean((Yk - pp_full_glmnet)^2)
  theta_glmboost <- mean((Yk - pp_reduced_glmboost)^2) - mean((Yk - pp_full_glmboost)^2)
  theta_rf <- mean(squared_loss(new_rf_pred$pred_full, new_rf_pred$pred_reduced, Yk))
  
  # New method thetas
  #theta_ridge <- mean((Yk - pp_reduced_ridge)^2) - mean((Yk - pp_full_ridge)^2)
  #theta_bcf <- mean((Yk - pp_reduced_bcf)^2) - mean((Yk - pp_full_bcf)^2)
  #theta_cf <- mean((Yk - pp_reduced_cf)^2) - mean((Yk - pp_full_cf)^2)
  
  #######################################################
  ###.  Getting confidence intervals and h-values
  ######################################################
  naive.trt.effect <- mean(Y[A==1]) - mean(Y[A==0])
  tau.range <- sort(c(-5*naive.trt.effect, 5*naive.trt.effect))
  XX <- data.frame(model.matrix(Y ~ x1 + x2 + A + A:x1 + A:x2 - 1))
  XX0 <- data.frame(model.matrix(Y ~ x1 + x2 - 1))
  XX_cf <- data.frame(x1 = x1, x2 = x2)  # For causal forest
  
  # Original method nested CV
  ncv_lm <- nested_cv(X=XX, X0=XX0, Y=as.vector(Y), Trt=A, tau.range=tau.range, funcs=linear_regression_funs,
                      n_folds = n_folds, reps  = nested_cv_reps)
  
  ncv_net <- nested_cv(X=XX, X0=XX0, Y=as.vector(Y), Trt=A, tau.range=tau.range, funcs=glmnet_funs,
                       n_folds = n_folds, reps  = nested_cv_reps)
  
  ncv_boost <- nested_cv(X=XX, X0=XX0, Y=as.vector(Y), Trt=A, tau.range=tau.range, funcs=glmboost_funs,
                         n_folds = n_folds, reps  = nested_cv_reps)
  
  ncv_rf <- nested_cv(X=XX, X0=XX0, Y=as.vector(Y), Trt=A, tau.range=tau.range, funcs=rf_funs,
                      n_folds = n_folds, reps = 10, bias_reps = 0)
  
  # New method nested CV
  #ncv_ridge <- nested_cv(X=XX, X0=XX0, Y=as.vector(Y), Trt=A, tau.range=tau.range, funcs=ridge_funs,
  #                       n_folds = n_folds, reps = nested_cv_reps)
  
  
  #ncv_bcf <- nested_cv(X=XX, X0=XX0, Y=as.vector(Y), Trt=A, tau.range=tau.range, funcs=bcf_funs,
  #                     n_folds = n_folds, reps = min(10, nested_cv_reps))
  
  
  #ncv_cf <- nested_cv(X=XX_cf, X0=XX_cf, Y=as.vector(Y), Trt=A, tau.range=tau.range, funcs=causal_forest_funs,
  #                    n_folds = n_folds, reps = min(10, nested_cv_reps))
  
  ############################################
  ### Record Results 
  ############################################

  
  # Original methods
  cover_lm[h] <- theta_lm > ncv_lm$ci_lo & theta_lm < ncv_lm$ci_hi
  CI_lm[h,1] <- ncv_lm$ci_lo
  CI_lm[h,2] <- ncv_lm$ci_hi
  hvalue_lm[h] <- ncv_lm$hvalue
  se_lm[h] <- ncv_lm$std_err 
  sd_infl_lm[h] <- ncv_lm$sd_infl
  bias_est_lm[h] <- ncv_lm$bias_est
  
  cover_glmnet[h] <- theta_glmnet > ncv_net$ci_lo & theta_glmnet < ncv_net$ci_hi
  CI_glmnet[h,1] <- ncv_net$ci_lo
  CI_glmnet[h,2] <- ncv_net$ci_hi
  hvalue_glmnet[h] <- ncv_net$hvalue
  se_glmnet[h] <- ncv_net$std_err
  sd_infl_glmnet[h] <- ncv_net$sd_infl
  bias_est_glmnet[h] <- ncv_net$bias_est
  
  cover_glmboost[h] <- theta_glmboost > ncv_boost$ci_lo & theta_glmboost < ncv_boost$ci_hi
  CI_glmboost[h,1] <- ncv_boost$ci_lo
  CI_glmboost[h,2] <- ncv_boost$ci_hi
  hvalue_glmboost[h] <- ncv_boost$hvalue
  se_glmboost[h] <- ncv_boost$std_err
  sd_infl_glmboost[h] <- ncv_boost$sd_infl
  bias_est_glmboost[h] <- ncv_boost$bias_est
  
  cover_rf[h] <- theta_rf > ncv_rf$ci_lo & theta_rf < ncv_rf$ci_hi
  CI_rf[h,1] <- ncv_rf$ci_lo
  CI_rf[h,2] <- ncv_rf$ci_hi
  hvalue_rf[h] <- ncv_rf$hvalue
  se_rf[h] <- ncv_rf$std_err
  sd_infl_rf[h] <- ncv_rf$sd_infl
  bias_est_rf[h] <- ncv_rf$bias_est
  
  # New methods
  #cover_ridge[h] <- theta_ridge > ncv_ridge$ci_lo & theta_ridge < ncv_ridge$ci_hi
  #CI_ridge[h,1] <- ncv_ridge$ci_lo
  #CI_ridge[h,2] <- ncv_ridge$ci_hi
  #hvalue_ridge[h] <- ncv_ridge$hvalue
  #se_ridge[h] <- ncv_ridge$std_err
  #sd_infl_ridge[h] <- ncv_ridge$sd_infl
  #bias_est_ridge[h] <- ncv_ridge$bias_est
  
  #cover_bcf[h] <- theta_bcf > ncv_bcf$ci_lo & theta_bcf < ncv_bcf$ci_hi
  #CI_bcf[h,1] <- ncv_bcf$ci_lo
  #CI_bcf[h,2] <- ncv_bcf$ci_hi
  #hvalue_bcf[h] <- ncv_bcf$hvalue
  #se_bcf[h] <- ncv_bcf$std_err
  #sd_infl_bcf[h] <- ncv_bcf$sd_infl
  #bias_est_bcf[h] <- ncv_bcf$bias_est
  
  #cover_causal_forest[h] <- theta_cf > ncv_cf$ci_lo & theta_cf < ncv_cf$ci_hi
  #CI_causal_forest[h,1] <- ncv_cf$ci_lo
  #CI_causal_forest[h,2] <- ncv_cf$ci_hi
  #hvalue_causal_forest[h] <- ncv_cf$hvalue
  #se_causal_forest[h] <- ncv_cf$std_err
  #sd_infl_causal_forest[h] <- ncv_cf$sd_infl
  #bias_est_causal_forest[h] <- ncv_cf$bias_est
  
  # Store all true thetas
  #true_thetas[h,] <- c(theta_lm, theta_glmnet, theta_glmboost, theta_rf, theta_ridge, theta_bcf, theta_cf)
  true_thetas[h,] <- c(theta_lm, theta_glmnet, theta_glmboost, theta_rf)
  
  cat("Simulation Replication:  ", h, "\n")
}

# Save all results (UPDATED to include standard errors and additional metrics)
save(
  # Original methods
  cover_lm, cover_glmnet, cover_glmboost, cover_rf,
  CI_lm, CI_glmnet, CI_glmboost, CI_rf,
  hvalue_lm, hvalue_glmnet, hvalue_glmboost, hvalue_rf,
  # Standard errors (properly calculated)
  se_lm, se_glmnet, se_glmboost, se_rf,
  # Additional nested CV diagnostics
  sd_infl_lm, sd_infl_glmnet, sd_infl_glmboost, sd_infl_rf,
  bias_est_lm, bias_est_glmnet, bias_est_glmboost, bias_est_rf,
  # New methods
  #cover_ridge, cover_bcf, cover_causal_forest,
  #CI_ridge, CI_bcf, CI_causal_forest,
  #hvalue_ridge, hvalue_bcf, hvalue_causal_forest,
  # Standard errors for new methods
  #se_ridge, se_bcf, se_causal_forest,
  # Additional nested CV diagnostics for new methods
  #sd_infl_ridge, sd_infl_bcf, sd_infl_causal_forest,
  #bias_est_ridge, bias_est_bcf, bias_est_causal_forest,
  # True thetas for all methods
  true_thetas,
  file = "Sim1B_100.RData"
)


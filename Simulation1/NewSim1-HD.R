library(glmnet)
library(mboost)
library(dbarts)
library(randomForest)
library(bcf)     # Bayesian Causal Forest
library(grf)     # Causal Forest
library(bartCause)  
library(rlearner)   

##################################
## Simulation Study 1B: High-Dimensional Setting (p=100) 
##################################
# Set seed for reproducibility
set.seed(12356)

# Define the coefficients
beta0 <- 2
p <- 100  # Total number of predictors
beta_main <- c(3, -1, 0.8, -0.6, 1.2, rep(0, p - 5))  # Main effects
beta_trt <- 1.5  # Treatment main effect
beta_interaction <- c(0.5, -2, 0.7, 3, 0.9, rep(0, p - 5))  # Treatment interactions

## Source fitting and nested cv functions
setwd("~/Projects/HTE-Model-Comparison")  ## Change for your computer
source("CoreFunctions/NewCVDiffFunctions.R")
source("CoreFunctions/NewFitterFunctions.R")

## Define all method function lists
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

ridge_funs <- list(fitter = fitter_ridge,
                   predictor = predictor_ridge,
                   mse = mse,
                   loss = squared_loss,
                   name = "ridge")

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

causal_forest_funs <- list(fitter = fitter_causal_forest,
                           predictor = predictor_causal_forest,
                           mse = mse,
                           loss = squared_loss,
                           name = "causal_forest")

bartcause_funs <- list(fitter = fitter_bartcause,
                       predictor = predictor_bartcause,
                       mse = mse,
                       loss = squared_loss,
                       name = "bartcause")

rlearner_lasso_funs <- list(fitter = function(X, X0, Y, Trt, tau.range, idx = NA) {
  fitter_rlearner(X, X0, Y, Trt, tau.range, idx, method = "rlasso")
},
predictor = predictor_rlearner,
mse = mse,
loss = squared_loss,
name = "rlearner_lasso")

# Set the number of observations n, number of folds, and
# number of nested cv replications:  
n <- 500
n_folds <- 5
nested_cv_reps <- 50 ## Use 50 or 100 for paper

## Set the number of simulation replications
nreps <- 1  ## Use nreps = 500 for paper

# Initialize storage for all methods
cover_lm <- cover_glmnet <- cover_ridge <- cover_glmboost <- cover_rf <- cover_bcf <- cover_causal_forest <- cover_bartcause <- cover_rlearner_lasso <- rep(NA, nreps)
hvalue_lm <- hvalue_glmnet <- hvalue_ridge <- hvalue_glmboost <- hvalue_rf <- hvalue_bcf <- hvalue_causal_forest <- hvalue_bartcause <- hvalue_rlearner_lasso <- rep(NA, nreps)
CI_lm <- CI_glmnet <- CI_ridge <- CI_glmboost <- CI_rf <- CI_bcf <- CI_causal_forest <- CI_bartcause <- CI_rlearner_lasso <- matrix(NA, nrow=nreps, ncol=2)
true_thetas <- matrix(NA, nreps, 8)  

for(h in 1:nreps) {
  
  # Generate random values for x1 to x100 from a normal distribution
  X_matrix <- matrix(rnorm(n * p, mean = 0, sd = 1), nrow = n, ncol = p)
  colnames(X_matrix) <- paste0("x", 1:p)
  
  # Generate treatment variable A
  A <- rbinom(n, size = 1, prob = 0.5)
  
  # Generate residuals and outcomes
  epsilon <- rnorm(n, mean = 0, sd = 1)
  
  # Calculate outcome 
  Y <- beta0 + 
    X_matrix %*% beta_main + 
    beta_trt * A + 
    A * (X_matrix %*% beta_interaction) + 
    epsilon
  Y <- as.vector(Y)  # Convert to vector
  
  ######################
  ## Compute true value of \theta_{XY} for each method
  ######################
  
  ## Create data frames
  # Full model with all interactions
  DAT <- data.frame('Y' = Y, X_matrix, 'A' = A)
  # Add interaction terms for all predictors
  for(i in 1:p) {
    DAT[paste0("x", i, ".t")] <- A * X_matrix[, i]
  }
  
  # Reduced model (main effects + treatment)
  DAT_reduced <- data.frame('Y' = Y, X_matrix, 'A' = A)
  
  # For methods that don't use treatment variable directly
  DAT_red <- data.frame(Wtau = Y, X_matrix)
  
  # Create design matrices
  # Full model matrix (main effects + treatment + all interactions)
  formula_full <- as.formula(paste("Y ~", paste(c(paste0("x", 1:p), "A", paste0("A: x", 1:p)), collapse = " + ")))
  Xmat_tmp <- model.matrix(formula_full, data = DAT)[, -1]  # Remove intercept
  
  # Reduced model matrix (main effects + treatment, no interactions)
  formula_reduced <- as.formula(paste("Y ~", paste(c(paste0("x", 1:p), "A"), collapse = " + ")))
  X0mat_tmp <- model.matrix(formula_reduced, data = DAT_reduced)[, -1]  # Remove intercept
  
  # Matrix without treatment (for glmboost and rf)
  formula_no_trt <- as.formula(paste("Y ~", paste(paste0("x", 1:p), collapse = " + ")))
  X0mat_notrt_tmp <- model.matrix(formula_no_trt, data = data.frame(X_matrix))[, -1]  # Remove intercept
  
  ## Fit models for theta computation
  tmp_lm <- lm(formula_full, data = DAT)
  tmp_reduced_lm <- lm(formula_reduced, data = DAT_reduced)
  
  tmp_glmnet <- cv.glmnet(Xmat_tmp, DAT$Y, family = "gaussian", nfolds = 5, alpha = 1)  # LASSO
  tmp_reduced_glmnet <- cv.glmnet(X0mat_tmp, DAT$Y, family = "gaussian", nfolds = 5, alpha = 1)
  
  tmp_ridge <- cv.glmnet(Xmat_tmp, DAT$Y, family = "gaussian", nfolds = 5, alpha = 0)  # Ridge
  tmp_reduced_ridge <- cv.glmnet(X0mat_tmp, DAT$Y, family = "gaussian", nfolds = 5, alpha = 0)
  
  tmp_glmboost <- glmboost(x=Xmat_tmp, y=DAT$Y, family = Gaussian())
  
  f0fn <- function(tau) {
    Wtau <- DAT_reduced$Y - tau*DAT_reduced$A
    fit_reduced <- glmboost(x=X0mat_notrt_tmp, y=Wtau, family = Gaussian())
    mse_local <- mean((Wtau - predict(fit_reduced, newdata=X0mat_notrt_tmp))^2)
    return(mse_local)
  }
  tau.star.gboost <- optimize(f0fn, interval=c(-5, 5))$minimum
  Wtau.star <- DAT_reduced$Y - tau.star.gboost*DAT_reduced$A
  tmp_reduced_glmboost <- glmboost(x=X0mat_notrt_tmp, y = Wtau.star, family = Gaussian())
  
  f00fn <- function(tau) {
    DAT_red$Wtau <- DAT_reduced$Y - tau*DAT_reduced$A
    fit_reduced <- randomForest(Wtau ~., data=DAT_red, maxnodes=16, ntree=100)
    mse_local <- mean((DAT_red$Wtau - predict(fit_reduced, newdata=DAT_red))^2)
    return(mse_local)
  }
  tau.star.rf <- optimize(f00fn, interval=c(-5, 5))$minimum
  DAT_red$Wtau <- DAT_reduced$Y - tau.star.rf*DAT_reduced$A
  tmp_reduced_rf <- randomForest(Wtau ~., data=DAT_red, maxnodes=16, ntree=100)

    # Fit Causal Forest model
  tmp_cf <- causal_forest(X = X_matrix,
                          Y = Y,
                          W = A,
                          num.trees = 100, 
                          honesty = TRUE,
                          honesty.fraction = 0.5,
                          ci.group.size = 2)
  
  # Reduced model for Causal Forest
  f_cf <- function(tau) {
    Wtau <- Y - tau * A
    fit_reduced <- regression_forest(X = X_matrix, Y = Wtau, num.trees = 100)
    pred_reduced <- predict(fit_reduced, X_matrix)$predictions
    mse_local <- mean((Wtau - pred_reduced)^2)
    return(mse_local)
  }
  tau.star.cf <- optimize(f_cf, interval=c(-5, 5))$minimum
  Wtau.star.cf <- Y - tau.star.cf * A
  tmp_reduced_cf <- regression_forest(X = X_matrix, Y = Wtau.star.cf, num.trees = 100)
  
  # bartCause model
  tmp_bartcause <- bartc(response = Y,
                         treatment = A,
                         confounders = X_matrix,
                         estimand = "ate",
                         method.rsp = "bart",
                         method.trt = "bart",
                         verbose = FALSE)
  
  # For reduced model, create a proper fitting function
  f_bartcause <- function(tau) {
    Wtau <- Y - tau * A
    # Use regression forest or lm for reduced model
    fit_reduced <- lm(Wtau ~ ., data = data.frame(X_matrix))
    mse_local <- mean((Wtau - fitted(fit_reduced))^2)
    return(mse_local)
  }
  tau.star.bartcause <- optimize(f_bartcause, interval=c(-5, 5))$minimum
  Wtau.star.bartcause <- Y - tau.star.bartcause * A
  tmp_reduced_bartcause <- lm(Wtau.star.bartcause ~ ., data = data.frame(X_matrix))

  
  # rlearner model
  tmp_rlearner <- rlasso(X_matrix, A, Y)
  
  # Calculate baseline for rlearner
  baseline_rl <- mean(Y[A == 0])
  
  f_rlearner <- function(tau) {
    Wtau <- Y - tau * A
    # Use proper data frame with column names that match future data
    fit_reduced <- lm(Wtau ~ ., data = data.frame(X_matrix))
    mse_local <- mean((Wtau - fitted(fit_reduced))^2)
    return(mse_local)
  }
  tau.star.rlearner <- optimize(f_rlearner, interval=c(-5, 5))$minimum
  Wtau.star.rlearner <- Y - tau.star.rlearner * A
  # Create proper data frame with matching column names
  tmp_reduced_rlearner <- lm(Wtau.star.rlearner ~ ., data = data.frame(X_matrix))
  
  ## Now, evaluate MSE difference on a much larger "future" dataset
  nr <- 100000
  
  ## Generating future observations:  
  epsilon_k <- rnorm(nr, mean = 0, sd = 1)
  X_matrix_k <- matrix(rnorm(nr * p, mean = 0, sd = 1), nrow = nr, ncol = p)
  colnames(X_matrix_k) <- paste0("x", 1:p)
  Ak <- rbinom(nr, size = 1, prob = 0.5)
  
  Yk <- beta0 + 
    X_matrix_k %*% beta_main + 
    beta_trt * Ak + 
    Ak * (X_matrix_k %*% beta_interaction) + 
    epsilon_k
  Yk <- as.vector(Yk)
  
  # Create future data frames
  DATk <- data.frame('Y' = Yk, X_matrix_k, 'A' = Ak)
  for(i in 1:p) {
    DATk[paste0("x", i, ".t")] <- Ak * X_matrix_k[, i]
  }
  
  DATk_reduced <- data.frame('Y' = Yk, X_matrix_k, 'A' = Ak)
  DATk_red <- data.frame('Y' = Yk, X_matrix_k)
  
  # Future design matrices - ensure consistent column naming
  Xmat_tmpk <- model.matrix(formula_full, data = DATk)[, -1]
  X0mat_tmpk <- model.matrix(formula_reduced, data = DATk_reduced)[, -1]
  X0mat_notrt_tmpk <- X_matrix_k  # Direct assignment - no model.matrix needed
  
  ## Getting predictions for these future observations
  pp_reduced_lm <- predict(tmp_reduced_lm, newdata=DATk_reduced)
  pp_full_lm <- predict(tmp_lm, newdata=DATk)
  
  pp_reduced_glmnet <- as.numeric(predict(tmp_reduced_glmnet, newx=X0mat_tmpk))
  pp_full_glmnet <- as.numeric(predict(tmp_glmnet, newx=Xmat_tmpk))
  
  pp_reduced_ridge <- as.numeric(predict(tmp_reduced_ridge, newx=X0mat_tmpk))
  pp_full_ridge <- as.numeric(predict(tmp_ridge, newx=Xmat_tmpk))
  
  pp_reduced_glmboost <- as.numeric(predict(tmp_reduced_glmboost, newdata=X0mat_notrt_tmpk)) + tau.star.gboost*DATk_reduced$A
  pp_full_glmboost <- as.numeric(predict(tmp_glmboost, newdata = Xmat_tmpk))
  
  # Use fitter functions for consistent theta computation
  new_rf <- fitter_rf(data.frame(Xmat_tmp), data.frame(X0mat_notrt_tmp), Y, A, tau.range=c(-5,5), idx = NA)
  new_rf_pred <- predictor_rf(new_rf, X_new=data.frame(Xmat_tmpk), X0_new=data.frame(X0mat_notrt_tmpk), Trt_new=DATk_reduced$A)
  
  # Causal Forest predictions
  mu_forest <- regression_forest(X_matrix, Y)
  mu_pred_k <- predict(mu_forest, X_matrix_k)$predictions
  tau_pred_k <- predict(tmp_cf, X_matrix_k)$predictions
  pp_full_cf <- mu_pred_k + tau_pred_k * Ak
  pp_reduced_cf <- predict(tmp_reduced_cf, X_matrix_k)$predictions + tau.star.cf * Ak
  
  # bartCause predictions
  mu_0_k <- extract(tmp_bartcause, "mu.0")  # Control outcomes
  mu_1_k <- extract(tmp_bartcause, "mu.1")  # Treatment outcomes
  
  # Get posterior means
  mu_0_mean_k <- colMeans(mu_0_k)
  mu_1_mean_k <- colMeans(mu_1_k)
  
  # For future data, we need to predict.  Since bartCause doesn't have a direct predict method,
  # we'll use the average treatment effect approach
  icate_k <- extract(tmp_bartcause, "icate")
  ate_estimate <- mean(colMeans(icate_k))
  baseline_bc <- mean(mu_0_mean_k)
  pp_full_bartcause <- baseline_bc + ate_estimate * Ak
  
  # Fixed prediction for reduced model - use future data
  pp_reduced_bartcause <- predict(tmp_reduced_bartcause, newdata = data.frame(X_matrix_k)) + tau.star.bartcause * Ak
  
  # rlearner predictions
  tau_pred_rl <- predict(tmp_rlearner, X_matrix_k)
  pp_full_rlearner <- baseline_rl + tau_pred_rl * Ak
  pp_reduced_rlearner <- predict(tmp_reduced_rlearner, newdata = data.frame(X_matrix_k)) + tau.star.rlearner * Ak
  
  
  ## Compute \theta_{XY} by looking at differences in MSE
  theta_lm <- mean((Yk - pp_reduced_lm)^2) - mean((Yk - pp_full_lm)^2)
  theta_glmnet <- mean((Yk - pp_reduced_glmnet)^2) - mean((Yk - pp_full_glmnet)^2)
  theta_ridge <- mean((Yk - pp_reduced_ridge)^2) - mean((Yk - pp_full_ridge)^2)
  theta_glmboost <- mean((Yk - pp_reduced_glmboost)^2) - mean((Yk - pp_full_glmboost)^2)
  theta_rf <- mean(squared_loss(new_rf_pred$pred_full, new_rf_pred$pred_reduced, Yk))
  theta_cf <- mean((Yk - pp_reduced_cf)^2) - mean((Yk - pp_full_cf)^2)
  theta_bartcause <- mean((Yk - pp_reduced_bartcause)^2) - mean((Yk - pp_full_bartcause)^2)
  theta_rlearner <- mean((Yk - pp_reduced_rlearner)^2) - mean((Yk - pp_full_rlearner)^2)
  
  #######################################################
  ###.   Getting confidence intervals and h-values
  ######################################################
  naive.trt.effect <- mean(Y[A==1]) - mean(Y[A==0])
  tau.range <- sort(c(-5*naive.trt.effect, 5*naive.trt.effect))
  
  # Create design matrices for nested CV
  XX <- data.frame(Xmat_tmp)
  XX0 <- data.frame(X0mat_notrt_tmp)
  
  # Nested CV for all methods
  ncv_lm <- nested_cv(X=XX, X0=XX0, Y=as.vector(Y), Trt=A, tau.range=tau.range, funcs=linear_regression_funs,
                      n_folds = n_folds, reps = nested_cv_reps, n_cores = 4)
  
  ncv_net <- nested_cv(X=XX, X0=XX0, Y=as.vector(Y), Trt=A, tau.range=tau.range, funcs=glmnet_funs,
                       n_folds = n_folds, reps = nested_cv_reps, n_cores = 4)
  
  ncv_ridge <- nested_cv(X=XX, X0=XX0, Y=as.vector(Y), Trt=A, tau.range=tau.range, funcs=ridge_funs,
                         n_folds = n_folds, reps = nested_cv_reps, n_cores = 4)
  
  ncv_boost <- nested_cv(X=XX, X0=XX0, Y=as.vector(Y), Trt=A, tau.range=tau.range, funcs=glmboost_funs,
                         n_folds = n_folds, reps = nested_cv_reps, n_cores = 4)
  
  ncv_rf <- nested_cv(X=XX, X0=XX0, Y=as.vector(Y), Trt=A, tau.range=tau.range, funcs=rf_funs,
                      n_folds = n_folds, reps = 10, bias_reps = 0, n_cores = 4)
  XX_cf <- data.frame(X_matrix)  
  XX0_cf <- data.frame(X_matrix) 
  ncv_cf <- nested_cv(X=XX, X0=XX0, Y=as.vector(Y), Trt=A, tau.range=tau.range, funcs=causal_forest_funs,
                      n_folds = n_folds, reps = nested_cv_reps, n_cores = 4) 
  
  ncv_bartcause <- nested_cv(X=XX, X0=XX0, Y=as.vector(Y), Trt=A, tau.range=tau.range, 
                             funcs=bartcause_funs, n_folds = n_folds, reps = nested_cv_reps, n_cores = 4)
  
  ncv_rlearner <- nested_cv(X=XX, X0=XX0, Y=as.vector(Y), Trt=A, tau.range=tau.range, 
                            funcs=rlearner_lasso_funs, n_folds = n_folds, reps = nested_cv_reps, n_cores = 4)
  
  ############################################
  ### Record Results
  ############################################
  cover_lm[h] <- theta_lm > ncv_lm$ci_lo & theta_lm < ncv_lm$ci_hi
  CI_lm[h,1] <- ncv_lm$ci_lo
  CI_lm[h,2] <- ncv_lm$ci_hi
  hvalue_lm[h] <- ncv_lm$hvalue
  
  cover_glmnet[h] <- theta_glmnet > ncv_net$ci_lo & theta_glmnet < ncv_net$ci_hi
  CI_glmnet[h,1] <- ncv_net$ci_lo
  CI_glmnet[h,2] <- ncv_net$ci_hi
  hvalue_glmnet[h] <- ncv_net$hvalue
  
  cover_ridge[h] <- theta_ridge > ncv_ridge$ci_lo & theta_ridge < ncv_ridge$ci_hi
  CI_ridge[h,1] <- ncv_ridge$ci_lo
  CI_ridge[h,2] <- ncv_ridge$ci_hi
  hvalue_ridge[h] <- ncv_ridge$hvalue
  
  cover_glmboost[h] <- theta_glmboost > ncv_boost$ci_lo & theta_glmboost < ncv_boost$ci_hi
  CI_glmboost[h,1] <- ncv_boost$ci_lo
  CI_glmboost[h,2] <- ncv_boost$ci_hi
  hvalue_glmboost[h] <- ncv_boost$hvalue
  
  cover_rf[h] <- theta_rf > ncv_rf$ci_lo & theta_rf < ncv_rf$ci_hi
  CI_rf[h,1] <- ncv_rf$ci_lo
  CI_rf[h,2] <- ncv_rf$ci_hi
  hvalue_rf[h] <- ncv_rf$hvalue

  
  cover_causal_forest[h] <- theta_cf > ncv_cf$ci_lo & theta_cf < ncv_cf$ci_hi
  CI_causal_forest[h,1] <- ncv_cf$ci_lo
  CI_causal_forest[h,2] <- ncv_cf$ci_hi
  hvalue_causal_forest[h] <- ncv_cf$hvalue
  
  cover_bartcause[h] <- theta_bartcause > ncv_bartcause$ci_lo & theta_bartcause < ncv_bartcause$ci_hi
  CI_bartcause[h,1] <- ncv_bartcause$ci_lo
  CI_bartcause[h,2] <- ncv_bartcause$ci_hi
  hvalue_bartcause[h] <- ncv_bartcause$hvalue
  
  cover_rlearner_lasso[h] <- theta_rlearner > ncv_rlearner$ci_lo & theta_rlearner < ncv_rlearner$ci_hi
  CI_rlearner_lasso[h,1] <- ncv_rlearner$ci_lo
  CI_rlearner_lasso[h,2] <- ncv_rlearner$ci_hi
  hvalue_rlearner_lasso[h] <- ncv_rlearner$hvalue
  
  # Store true thetas for all methods
  true_thetas[h,] <- c(theta_lm, theta_glmnet, theta_ridge, theta_glmboost, theta_rf, theta_cf, theta_bartcause, theta_rlearner)
  
  cat("Simulation Replication:   ", h, "\n")
}


save(
  cover_lm, cover_glmnet, cover_ridge, cover_glmboost, cover_rf, cover_causal_forest, 
  cover_bartcause, cover_rlearner_lasso,
  CI_lm, CI_glmnet, CI_ridge, CI_glmboost, CI_rf, CI_causal_forest,
  CI_bartcause, CI_rlearner_lasso,
  hvalue_lm, hvalue_glmnet, hvalue_ridge, hvalue_glmboost, hvalue_rf, hvalue_causal_forest,
  hvalue_bartcause, hvalue_rlearner_lasso,
  true_thetas,
  beta_main, beta_interaction, beta_trt, beta0,
  file = "Sim1B_HighDim_Extended_500. RData"
)


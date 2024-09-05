
##################################
## Simulation Study 2

# Define mu and theta
f1 <- function(x) {
  return((x[, 2]*x[, 4]*x[, 6]) + (2*x[, 2]*x[, 4]*(1-x[, 6])) + (3*x[, 2]*(1-x[, 4])*x[, 6]) + (4*x[, 2]*(1-x[, 4])*(1-x[, 4])) + (5*(1-x[, 2])*x[, 4]*x[, 6]) +
           (6*(1-x[, 2])*x[, 4]*(1-x[, 6])) + (7*(1-x[, 2])*(1-x[, 4])*x[, 6]) + (8*(1-x[, 2])*(1-x[, 4])*(1-x[, 6])))
}

f4 <- function(x) {
  return((4*ifelse(x[, 1] > 1 & x[, 3] > 0, 1, 0)) + (4*ifelse(x[, 5] > 1 & x[, 7] > 1, 1, 0)) + (x[, 8]*x[, 9]))
}

f3 <- function(x) {
  return(0.5*(x[, 1]^2 + x[, 2] + x[, 3]^2 + x[, 4] + x[, 5]^2 + x[, 6] + x[, 7]^2 + x[, 8] + x[, 9]^2 -11))
}

f2 <- function(x) {
  #(1/sqrt(2))*(f1(x) + (x[, 1] + x[, 3] + x[, 5] + x[, 7] + x[, 8] + x[, 9] - 2))
  aa <- (1/sqrt(2))*(x[, 1] + x[, 3] + x[, 5] + x[, 7] + x[, 8] + x[, 9] - 2)
  return(aa)
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
    return(rep(0, nrow(x)))
  } else if (choice == 2) {
    return(rep(1, nrow(x)))
  } else if (choice == 3) {
    #return(2 + 0.1/(1 + exp(-x[, 2])))
    return(2 + 2/(1 + exp(-x[,2])))
  } else if (choice == 4) {
    return(f1(x))
  } else {
    stop("Invalid choice for tau")
  }
}

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

OneSimTwoRun <- function(n, mu_choice, theta_choice, p, n_folds, nested_cv_reps) {

nreps <- 1
cover_lm <- cover_glmnet <- cover_glmboost <- cover_rf <- rep(NA, nreps)
hvalue_lm <- hvalue_glmnet <- hvalue_glmboost <- hvalue_rf <- rep(NA, nreps)
hvalue1_lm <- hvalue1_glmnet <- hvalue1_glmboost <- hvalue1_rf <- rep(NA, nreps)
CI_lm <- CI_glmnet <- CI_glmboost <- CI_rf <- matrix(NA, nrow=nreps, ncol=2)
true_thetas <- matrix(NA, nreps, 4)
colnames(true_thetas) <- c("LM", "Lasso", "Boost", "RF")
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
  mu_new <- mu(mu_choice, x)
  theta_new <- theta(theta_choice, x)

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
  drop_cols <- c(which(colnames(DAT_reduced)=="Y"), which(colnames(DAT_reduced)=="A"))
  DAT_red <- data.frame(Wtau=DAT_reduced$Y, x)

  Xmat_tmp <- model.matrix(Y~.-1, data = DAT)
  X0mat_tmp <- model.matrix(Y~.-1, data = DAT_reduced)
  ## Don't use treatment variable in GLMboost
  X0mat_notrt_tmp <- model.matrix(Y ~ . - A - 1, data = DAT_reduced)

  tmp_lm <- lm(Y ~., data=DAT)
  tmp_reduced_lm <- lm(Y ~., data=DAT_reduced)
  tmp_glmnet <- cv.glmnet(Xmat_tmp, DAT$Y, family = "gaussian", nfolds = 5)

  f01fn <- function(tau) {
    Wtau <- DAT_reduced$Y - tau*DAT_reduced$A
    fit_reduced <- cv.glmnet(X0mat_notrt_tmp, Wtau, family = "gaussian", nfolds = 5)
    mse_local <- mean((Wtau - predict(fit_reduced, newx=X0mat_notrt_tmp))^2)
    return(mse_local)
  }
  tau.star.glmnet <- optimize(f01fn, interval=c(-20, 20))$minimum

  Wtau.star <- DAT_reduced$Y - tau.star.glmnet*DAT_reduced$A
  tmp_reduced_glmnet <- cv.glmnet(X0mat_notrt_tmp, Wtau.star, family = "gaussian", nfolds = 5)

  tmp_glmboost <- glmboost(x=Xmat_tmp, y=DAT$Y, family = Gaussian())
  f0fn <- function(tau) {
    Wtau <- DAT_reduced$Y - tau*DAT_reduced$A
    fit_reduced <- glmboost(x=X0mat_notrt_tmp, y=Wtau, family = Gaussian())
    mse_local <- mean((Wtau - predict(fit_reduced, newdata=X0mat_notrt_tmp))^2)
    return(mse_local)
  }
  tau.star.gboost <- optimize(f0fn, interval=c(-20, 20))$minimum
  Wtau.star <- DAT_reduced$Y - tau.star.gboost*DAT_reduced$A
  tmp_reduced_glmboost <- glmboost(x=X0mat_notrt_tmp, y = Wtau.star, family = Gaussian())

  f00fn <- function(tau) {
    DAT_red$Wtau <- DAT_reduced$Y - tau*DAT_reduced$A

    fit_reduced <- randomForest(Wtau ~., data=DAT_red, maxnodes=16, ntree=100)
    mse_local <- mean((DAT_red$Wtau - predict(fit_reduced, newdata=DAT_red))^2)
    return(mse_local)
  }
  tau.star.rf <- optimize(f00fn, interval=c(-20, 20))$minimum
  DAT_red$Wtau <- DAT_reduced$Y - tau.star.rf*DAT_reduced$A
  ## Is the line below correct?
  tmp_reduced_rf <- randomForest(Wtau ~., data=DAT_red, maxnodes=16, ntree=100)

  ## Now, evaluate MSE difference on a much larger "future" dataset
  nr <- 1000000

  ## Generating future observations:
  epsilon_future <- rnorm(nr, mean = 0, sd = 1)
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

  mu_future <- mu(mu_choice, xk)
  theta_future <- theta(theta_choice, xk)
  Yk <- mu_future + Ak*theta_future + epsilon_future

  x.tk <- matrix(NA, nrow = nrow(xk), ncol = ncol(xk))
  for(k in 1:ncol(xk)){
    x.tk[,k] <- xk[,k]*Ak
  }
  colnames(x.tk) <- paste0(colnames(xk), ".t")

  DATk <- data.frame('Y' = Yk, 'A' = Ak, xk, x.tk)
  DATk_reduced <- data.frame('Y' = Yk, 'A' = Ak, xk)
  DATk_red <- data.frame('Y' = Yk, xk)
  Xmat_tmpk <- model.matrix(Y~.-1, data = DATk)
  X0mat_tmpk <- model.matrix(Y~.-1, data = DATk_reduced)
  X0mat_notrt_tmpk <- model.matrix(Y ~ . - A - 1, data = DATk_reduced)

  ## Getting predictions for these future observations
  pp_reduced_lm <- predict(tmp_reduced_lm, newdata=DATk_reduced)
  pp_full_lm <- predict(tmp_lm, newdata=DATk)
  pp_reduced_glmnet <- as.numeric(predict(tmp_reduced_glmnet, newx=X0mat_notrt_tmpk)) + tau.star.glmnet*DATk_reduced$A
  pp_full_glmnet <- as.numeric(predict(tmp_glmnet, newx=Xmat_tmpk))
  pp_reduced_glmboost <- as.numeric(predict(tmp_reduced_glmboost, newdata=X0mat_notrt_tmpk)) + tau.star.gboost*DATk_reduced$A
  pp_full_glmboost <- as.numeric(predict(tmp_glmboost, newdata = Xmat_tmpk))
  #pp_reduced_rf <- predict(tmp_reduced_rf, newdata=DATk_red) + tau.star.rf*DATk_reduced$A
  #pp_full_rf <- predict(tmp_rf, newdata=DATk_reduced)
  new_rf <- fitter_rf(X0mat_tmp, X0mat_notrt_tmp, Y, A, tau.range=c(-20,20), idx = NA)
  new_rf_pred <- predictor_rf(new_rf, X_new=DATk_reduced, X0_new=DATk_red, Trt_new=DATk_reduced$A)
  
  ## Compute \theta_{XY} by looking at differences in MSE
  theta_lm <- mean((Yk - pp_reduced_lm)^2) - mean((Yk - pp_full_lm)^2)
  theta_glmnet <- mean((Yk - pp_reduced_glmnet)^2) - mean((Yk - pp_full_glmnet)^2)
  theta_glmboost <- mean((Yk - pp_reduced_glmboost)^2) - mean((Yk - pp_full_glmboost)^2)
  theta_rf <- mean(squared_loss(new_rf_pred$pred_full, new_rf_pred$pred_reduced, Yk))
  
  #######################################################
  ###. Getting confidence intervals and h-values
  ######################################################
  naive.trt.effect <- max(abs(mean(Y[A==1]) - mean(Y[A==0])), 2)
  tau.range <- c(-5*naive.trt.effect, 5*naive.trt.effect)
  XX <- data.frame(model.matrix(Y ~. - 1, data = DAT))
  XX0 <- data.frame(model.matrix(Y ~. - A - 1, data = DAT_reduced))
  ncv_lm <- nested_cv(X=XX, X0=XX0, Y=as.vector(Y), Trt=A, tau.range=tau.range, funcs=linear_regression_funs,
                      n_folds = n_folds, reps  = nested_cv_reps)

  ncv_net <- nested_cv(X=XX, X0=XX0, Y=as.vector(Y), Trt=A, tau.range=tau.range, funcs=glmnet_funs,
                       n_folds = n_folds, reps  = nested_cv_reps)

  ncv_boost <- nested_cv(X=XX, X0=XX0, Y=as.vector(Y), Trt=A, tau.range=c(-20,20), funcs=glmboost_funs,
                         n_folds = n_folds, reps  = nested_cv_reps)

  cirf_lo <- cirf_hi <- cirf_pe <- rep(NA, 10)
  for(u in 1:10) {
      tmp_ncv_rf <- nested_cv(X=XX, X0=XX0, Y=as.vector(Y), Trt=A, tau.range=c(-20, 20), funcs=rf_funs,
                              n_folds = n_folds, reps = 10, bias_reps=10)
      cirf_lo[u] <- tmp_ncv_rf$ci_lo
      cirf_hi[u] <- tmp_ncv_rf$ci_hi
      cirf_pe[u] <- (tmp_ncv_rf$ci_lo + tmp_ncv_rf$ci_hi)/2
  }
  ncv_rf <- list()
  ncv_rf$ci_lo <- min(cirf_lo)
  ncv_rf$ci_hi <- max(cirf_hi)
  rfse <- (max(cirf_hi) - min(cirf_lo))/4
  rfpe <- mean(cirf_pe)
  
  ncv_rf$hvalue <- 2*pnorm(-abs(rfpe)/rfse)
  ncv_rf$hvalue1sided <- pnorm(rfpe/rfse, lower.tail=FALSE)
  
  #ncv_rf <- nested_cv(X=XX, X0=XX0, Y=as.vector(Y), Trt=A, tau.range=c(-20,20), funcs=rf_funs,
  #                    n_folds = n_folds, reps = 10, bias_reps=10)
  ############################################
  ### Record Results
  ############################################
  cover_lm[h] <- theta_lm > ncv_lm$ci_lo & theta_lm < ncv_lm$ci_hi
  CI_lm[h,1] <- ncv_lm$ci_lo
  CI_lm[h,2] <- ncv_lm$ci_hi
  true_thetas[h,] <- c(theta_lm, theta_glmnet, theta_glmboost, theta_rf)
  hvalue_lm[h] <- ncv_lm$hvalue
  hvalue1_lm[h] <- ncv_lm$hvalue1sided
  
  cover_glmnet[h] <- theta_glmnet > ncv_net$ci_lo & theta_glmnet < ncv_net$ci_hi
  CI_glmnet[h,1] <- ncv_net$ci_lo
  CI_glmnet[h,2] <- ncv_net$ci_hi
  hvalue_glmnet[h] <- ncv_net$hvalue
  hvalue1_glmnet[h] <- ncv_net$hvalue1sided
  
  cover_glmboost[h] <- theta_glmboost > ncv_boost$ci_lo & theta_glmboost < ncv_boost$ci_hi
  CI_glmboost[h,1] <- ncv_boost$ci_lo
  CI_glmboost[h,2] <- ncv_boost$ci_hi
  hvalue_glmboost[h] <- ncv_boost$hvalue
  hvalue1_glmboost[h] <- ncv_boost$hvalue1sided
  
  cover_rf[h] <- theta_rf > ncv_rf$ci_lo & theta_rf < ncv_rf$ci_hi
  CI_rf[h,1] <- ncv_rf$ci_lo
  CI_rf[h,2] <- ncv_rf$ci_hi
  hvalue_rf[h] <- ncv_rf$hvalue
  hvalue1_rf[h] <- ncv_rf$hvalue1sided
  
  #cat("Simulation Replication: ", h, "\n")
}
ans <- list()
ans$true_thetas <- true_thetas
ans$cover_lm <- cover_lm
ans$CI_lm <- CI_lm
ans$hvalue_lm <- hvalue_lm
ans$hvalue1side_lm <- hvalue1_lm
ans$cover_glmnet <- cover_glmnet
ans$CI_glmnet <- CI_glmnet
ans$hvalue_glmnet <- hvalue_glmnet
ans$hvalue1side_glmnet <- hvalue1_glmnet
ans$cover_glmboost <- cover_glmboost
ans$CI_glmboost <- CI_glmboost
ans$hvalue_glmboost <- hvalue_glmboost
ans$hvalue1side_glmboost <- hvalue1_glmboost
ans$cover_rf <- cover_rf
ans$CI_rf <- CI_rf
ans$hvalue_rf <- hvalue_rf
ans$hvalue1side_rf <- hvalue1_rf
return(ans)
}
## Save: cover_lm, cover_glmnet, cover_glmboost, cover_rf
##       CI_lm, CI_glmnet, CI_glmboost, CI_rf
##      hvalue_lm, hvalue_glmnet, hvalue_glmboost, hvalue_bart


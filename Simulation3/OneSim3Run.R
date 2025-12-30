
##################################
## Simulation Study 3

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
                               loss = wsquared_loss,
                               name = "linear_regression")

glmnet_funs <- list(fitter = fitter_glmnet,
                    predictor = predictor_glmnet,
                    loss = wsquared_loss,
                    name = "glmnet")

glmboost_funs <- list(fitter = fitter_glmboost,
                      predictor = predictor_glmboost,
                      loss = wsquared_loss,
                      name = "glmboost")

ridge_funs <- list(fitter = fitter_ridge,
                predictor = predictor_ridge,
                loss = wsquared_loss,
                name = "ridge")

OneSimThreeRun <- function(n, mu_choice, theta_choice, p, n_folds, nested_cv_reps) {
  
  nreps <- 1
  cover_lm <- cover_glmnet <- cover_glmboost <- cover_ridge <- rep(NA, nreps)
  hvalue_lm <- hvalue_glmnet <- hvalue_glmboost <- hvalue_ridge <- rep(NA, nreps)
  hvalue1_lm <- hvalue1_glmnet <- hvalue1_glmboost <- hvalue1_ridge <- rep(NA, nreps)
  CI_lm <- CI_glmnet <- CI_glmboost <- CI_ridge <- matrix(NA, nrow=nreps, ncol=2)
  true_thetas <- matrix(NA, nreps, 4)
  colnames(true_thetas) <- c("LM", "Lasso", "Boost", "Ridge")
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
    W <- 2*A*Y - 2*(1 - A)*Y
    
    ######################
    ## Compute true value of \tilde{\theta}_{XY} for each method
    ######################
    
    ## First fit each model with the dataset of size n
    ## Change to W here.
    DAT <- data.frame('W' = W, x)
    
    Xmat_tmp <- model.matrix(W~.-1, data = DAT)
    
    tmp_lm <- lm(W ~., data=DAT)
    tmp_glmnet <- cv.glmnet(Xmat_tmp, DAT$W, family = "gaussian", nfolds = 5)
    tmp_glmboost <- glmboost(x=Xmat_tmp, y=DAT$W, family = Gaussian())
    tmp_ridge <- cv.glmnet(Xmat_tmp, DAT$W, family = "gaussian", nfolds = 5, alpha=1)
    tau.star.observed <- mean(DAT$W)
    
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
    Wk <- 2*Ak*Yk - 2*(1 - Ak)*Yk
    
    ## Change to W here.
    DATk_red <- data.frame('W' = Wk, xk)
    Xmat_tmpk <- model.matrix(W ~ . - 1, data = DATk_red)
    
    ## Getting predictions for these future observations
    pp_lm <- predict(tmp_lm, newdata=DATk_red)
    pp_glmnet <- as.numeric(predict(tmp_glmnet, newx=Xmat_tmpk))
    pp_glmboost <- as.numeric(predict(tmp_glmboost, newdata = Xmat_tmpk))
    pp_ridge <- as.numeric(predict(tmp_ridge, newx=Xmat_tmpk))
   
    ## Compute \tilde{\theta}_{XY} by looking at differences in MSE
    theta_lm <- mean((theta_future - pp_lm)^2) - mean((theta_future - tau.star.observed)^2)
    theta_glmnet <- mean((theta_future - pp_glmnet)^2) - mean((theta_future - tau.star.observed)^2)
    theta_glmboost <- mean((theta_future - pp_glmboost)^2) - mean((theta_future - tau.star.observed)^2)
    theta_ridge <- mean((theta_future - pp_ridge)^2) - mean((theta_future - tau.star.observed)^2)
    
    #######################################################
    ###. Getting confidence intervals and h-values
    ######################################################
    XX <- data.frame(model.matrix(W ~. - 1, data = DAT))
    ncv_lm <- nested_cv(X=XX, W=DAT$W, funcs=linear_regression_funs,
                        n_folds = n_folds, reps  = nested_cv_reps)
    
    ncv_net <- nested_cv(X=XX, W=DAT$W, funcs=glmnet_funs,
                         n_folds = n_folds, reps  = nested_cv_reps)
    
    ncv_boost <- nested_cv(X=XX, W=DAT$W, funcs=glmboost_funs,
                           n_folds = n_folds, reps  = nested_cv_reps)
    
    ncv_ridge <- nested_cv(X=XX, W=DAT$W, funcs=ridge_funs,
                           n_folds = n_folds, reps  = nested_cv_reps)
    
    ############################################
    ### Record Results
    ############################################
    cover_lm[h] <- theta_lm > ncv_lm$ci_lo & theta_lm < ncv_lm$ci_hi
    CI_lm[h,1] <- ncv_lm$ci_lo
    CI_lm[h,2] <- ncv_lm$ci_hi
    true_thetas[h,] <- c(theta_lm, theta_glmnet, theta_glmboost, theta_ridge)
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
    
    cover_ridge[h] <- theta_ridge > ncv_ridge$ci_lo & theta_ridge < ncv_ridge$ci_hi
    CI_ridge[h,1] <- ncv_ridge$ci_lo
    CI_ridge[h,2] <- ncv_ridge$ci_hi
    hvalue_ridge[h] <- ncv_ridge$hvalue
    hvalue1_ridge[h] <- ncv_ridge$hvalue1sided
    
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
  ans$cover_ridge <- cover_ridge
  ans$CI_ridge <- CI_ridge
  ans$hvalue_ridge <- hvalue_ridge
  ans$hvalue1side_ridge <- hvalue1_ridge
  ans$cover_glmboost <- cover_glmboost
  ans$CI_glmboost <- CI_glmboost
  ans$hvalue_glmboost <- hvalue_glmboost
  ans$hvalue1side_glmboost <- hvalue1_glmboost
  return(ans)
}
## Save: cover_lm, cover_glmnet, cover_glmboost, cover_rf
##       CI_lm, CI_glmnet, CI_glmboost, CI_rf
##      hvalue_lm, hvalue_glmnet, hvalue_glmboost, hvalue_bart


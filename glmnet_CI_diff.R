
library(glmnet)

squared_loss <- function(y1, y2, y3, tau, Trt) {
  ## y1 - full model
  ## y2 - reduced model 
  ## y3 - outcome
  (y1 - y3)^2 - (y2 - (y3 - tau*Trt))^2
}

fitter_glmnet <- function(X, X0, Y, Trt, tau.seq, idx = NA, funcs_params = NA) {
  ## X0 does not have treatment column
  if(sum(is.na(idx)) > 0) {
    idx <- 1:nrow(X)
  }
  fit <- glmnet(X[idx, ], Y[idx], family = "gaussian", lambda = funcs_params$lambdas) 
  
  mse <- rep(NA, length(tau.seq))
  for(k in 1:length(tau.seq)) {
    glmnet_tmp <- glmnet(X0[idx, ], Y[idx]-(tau.seq[k]*Treat[idx]), family = "gaussian", funcs_params$lambdas) 
    mse[k] <- mean(((Y[idx]-(tau.seq[k]*Treat[idx])) - glmnet_tmp$fitted)^2)
  }
  tau.star <- tau.seq[which.min(mse)]
  

  
  fit_reduced <- glmnet(X0[idx, ], Y[idx]-(tau.star*Treat[idx]), family = "gaussian", lambda = funcs_params$lambdas) 
  
  return(list(full=fit, reduced=fit_reduced, tau = tau.star))
}

predictor_glmnet <- function(fit, X_new, funcs_params = NA) {
  beta_hat <- fit$beta[, funcs_params$best_lam] 
  a0_hat <- fit$a0[funcs_params$best_lam]
  preds <- (as.matrix(X_new) %*% beta_hat + a0_hat)
  preds
}

gaussian_lasso_funs <- list(fitter = fitter_glmnet,
                            predictor = predictor_glmnet,
                            loss = squared_loss,
                            name = "gaussian_lasso")

n_folds <- 6
nested_cv_reps <- 5000 #average over many random splits


DATA <- read.csv('IHDP_clean.csv', header = TRUE)[,-1]
cor_matrix <- cor(DATA)
cor_with_output <- cor_matrix[, 1]
cor_with_output[1] <- 0  # Set the correlation with output to 0
sorted_correlations <- sort(abs(cor_with_output), decreasing = TRUE)
top_10_correlations <- sorted_correlations[1:10]
top_10_predictors <- names(top_10_correlations)
# Subset your data to keep only the top 10 predictor columns

selected_data <- DATA[, top_10_predictors]
output <- DATA[,c('iqsb.36', 'treat')]

## X0 - reduced model
DATA.cor.reduced <- bind_cols(output, selected_data)
# Create training and test sets using the selected index numbers
Y <- DATA.cor.reduced$iqsb.36
Treat <- DATA.cor.reduced$treat
tau <- mean(Treat == 1) - mean(Treat == 0)
tau.range <- seq(-4*tau, 4*tau, length.out = 9)
DATA.cor.reduced <- DATA.cor.reduced[ , !(names(DATA.cor.reduced) %in% c('iqsb.36', 'treat'))]
DATA.cor.reduced <- model.matrix(Y~.-1, data = DATA.cor.reduced)

## X - full model
selected_data.T <- matrix(NA, nrow = nrow(selected_data), ncol = ncol(selected_data))
for(k in 1:10){
  selected_data.T[,k] <- selected_data[,k]*output$treat
}
colnames(selected_data.T) <- paste0(colnames(selected_data), ".t")
DATA.cor <- bind_cols(output, selected_data, selected_data.T)
DATA.cor <- bind_cols('iqsb.36' = output, selected_data)

# Create training and test sets using the selected index numbers
Y <- DATA.cor$iqsb.36
DATA.cor <- DATA.cor[ , !(names(DATA.cor) %in% c('iqsb.36'))]
DATA.cor <- model.matrix(Y~.-1, data = DATA.cor)

fit <- cv.glmnet(as.matrix(DATA.cor.reduced), Y, Treat, family = "gaussian")
lambdas <- fit$lambda
best_lam <- match(fit$lambda.1se, lambdas) #selected value of lambda
lambda <- lambdas[1:best_lam]

fit <- cv.glmnet(as.matrix(DATA.cor), Y, family = "gaussian")
lambdas <- fit$lambda
best_lam <- match(fit$lambda.1se, lambdas) #selected value of lambda
lambda <- lambdas[1:best_lam]

nested_cv(data.frame(DATA.cor), data.frame(DATA.cor.reduced), as.vector(Y), as.vector(Treat), tau.seq = tau.range, gaussian_lasso_funs, 
            n_folds = n_folds, reps  = nested_cv_reps, verbose = T)

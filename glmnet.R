
library(dplyr)

######################
## Reduced model
######################

misclass_loss <- function(y1, y2, t, tau, funcs_params = NA) {
  (y1 - (y2-(tau*t)))^2
  } 

fitter_glmnet <- function(X, Y, Treat, tau, idx = NA, funcs_params = NA) {
  if(sum(is.na(idx)) > 0) {idx <- 1:nrow(X)}
  fit <- glmnet(X[idx, ], Y[idx]-(tau*Treat[idx]), family = "gaussian", lambda = funcs_params$lambdas) 
  
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
DATA.cor <- bind_cols(output, selected_data)

set.seed(123)
train_idx <- sample(1:nrow(DATA.cor), round(.7 * nrow(DATA.cor)), replace = FALSE)
test_idx <- setdiff(1:nrow(DATA.cor), train_idx)

# Create training and test sets using the selected index numbers
Y <- DATA.cor$iqsb.36
treat <- DATA.cor$treat
DATA.cor <- DATA.cor[ , !(names(DATA.cor) %in% c('iqsb.36', 'treat'))]
DATA.cor <- model.matrix(Y~.-1, data = DATA.cor)

train.set <- DATA.cor[train_idx, ]
Y.train <- Y[train_idx]
Treat.train <- treat[train_idx]
test.set <- DATA.cor[test_idx, ]
Y.test <-  Y[test_idx]
Treat.test <- treat[test_idx]


library(glmnet)
#Fit one model to find a good lambda. This lambda will be fixed in future simulations.
fit <- cv.glmnet(train.set, Y.train, Treat.train, family = "gaussian")
lambdas <- fit$lambda
best_lam <- match(fit$lambda.1se, lambdas) #selected value of lambda
lambda <- lambdas[1:best_lam]

tau.range = seq(1,10, by =1)
nested_cv_m(data.frame(train.set), as.vector(Y.train), as.vector(Treat.train), tau.range, gaussian_lasso_funs, 
          n_folds = n_folds, reps  = nested_cv_reps, 
          funcs_params = list("lambdas" = lambdas, "best_lam" = best_lam), verbose = T, alpha = 0.01)



#######################
### Full model
#######################
misclass_loss <- function(y1, y2, funcs_params = NA) {
  (y1 - y2)^2
} 

fitter_glmnet <- function(X, Y,  idx = NA, funcs_params = NA) {
  if(sum(is.na(idx)) > 0) {idx <- 1:nrow(X)}
  fit <- glmnet(X[idx, ], Y[idx], family = "gaussian", lambda = funcs_params$lambdas) 
  
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

output <- DATA[, c('iqsb.36', 'treat')]
selected_data.T <- matrix(NA, nrow = nrow(selected_data), ncol = ncol(selected_data))
for(k in 1:10){
  selected_data.T[,k] <- selected_data[,k]*output$treat
}
colnames(selected_data.T) <- paste0(colnames(selected_data), ".t")
DATA.cor <- bind_cols('iqsb.36' = output, selected_data)

set.seed(123)
train_idx <- sample(1:nrow(DATA.cor), round(.7 * nrow(DATA.cor)), replace = FALSE)
test_idx <- setdiff(1:nrow(DATA.cor), train_idx)

# Create training and test sets using the selected index numbers
Y <- DATA.cor$iqsb.36
DATA.cor <- DATA.cor[ , !(names(DATA.cor) %in% c('iqsb.36'))]
DATA.cor <- model.matrix(Y~.-1, data = DATA.cor)

train.set <- DATA.cor[train_idx, ]
Y.train <- Y[train_idx]
test.set <- DATA.cor[test_idx, ]
Y.test <-  Y[test_idx]


library(glmnet)
#Fit one model to find a good lambda. This lambda will be fixed in future simulations.
fit <- cv.glmnet(as.matrix(train.set), Y.train, family = "gaussian")
lambdas <- fit$lambda
best_lam <- match(fit$lambda.1se, lambdas) #selected value of lambda
lambda <- lambdas[1:best_lam]


nested_cv(data.frame(train.set), as.vector(Y.train), gaussian_lasso_funs, 
          n_folds = n_folds, reps  = nested_cv_reps, 
          funcs_params = list("lambdas" = lambdas, "best_lam" = best_lam), verbose = T, alpha = 0.01)


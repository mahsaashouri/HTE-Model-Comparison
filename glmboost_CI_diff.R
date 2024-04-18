
library(mboost)
library(dplyr)

squared_loss <- function(y1, y2, y3, Trt, tau) {
  ## y1 - full model
  ## y2 - reduced model 
  ## y3 - outcome
  (y1 - y3)^2 - (y2 - (y3 - tau*Trt))^2
}

fitter_glmboost <- function(X, X0, Y, Trt, tau.seq, idx = NA) {
  ## X0 does not have treatment column
  if(sum(is.na(idx)) > 0) {
    idx <- 1:nrow(X)
  }
  fit <- glmboost(x= as.matrix(X[idx, ]), y=Y[idx],family = Gaussian())
  
  mse <- rep(NA, length(tau.seq))
  for(k in 1:length(tau.seq)) {
    glmboost_tmp <- glmboost(x= as.matrix(X0[idx, ]), y=Y[idx]-(tau.seq[k]*Trt[idx]),family = Gaussian())
    mse[k] <- mean((Y[idx]-(tau.seq[k]*Trt[idx]) - glmboost_tmp$fitted())^2)
  }
  tau.star <- tau.seq[which.min(mse)]
  
  fit_reduced <- glmboost(x= as.matrix(X0[idx, ]), y=Y[idx]-(tau.star*Trt[idx]),family = Gaussian())
  
  return(list(full=fit, reduced=fit_reduced, tau = tau.star))
}

predictor_glmboost <- function(fit, X_new) {
  as.numeric(predict(fit, newdata = as.matrix(X_new)))
}

glmboost_funs <- list(fitter = fitter_glmboost,
                      predictor = predictor_glmboost,
                      loss = squared_loss,
                      name = "glmboost")
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


nested_cv(data.frame(DATA.cor), data.frame(DATA.cor.reduced), as.vector(Y), as.vector(Treat),tau.seq = tau.range, glmboost_funs, 
          n_folds = n_folds, reps  = nested_cv_reps, verbose = T)

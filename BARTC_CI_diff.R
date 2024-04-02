
library(dplyr)
library(bartCause)

# Function for squared loss
squared_loss <- function(y1, y2, y3, tau, Trt, funcs_params = NA) {
  ## y1 - full model
  ## y2 - reduced model 
  ## y3 - outcome
  (y1 - y3)^2 - (y2 - (y3 - tau*Trt))^2
}

# Function for fitting BART Cause model
fitter_bartc <- function(X, X0, Y, Trt, tau.seq, idx = NA, funcs_params = NA) {
  if (sum(is.na(idx)) > 0) {
    idx <- 1:nrow(X)
  }
  fit <- bartc(Y[idx], Trt[idx], X[idx, ], keepTrees = TRUE, n.samples = 100L, n.burn = 15L, n.chains = 2L)  
  
  mse <- rep(NA, length(tau.seq))
  for(k in 1:length(tau.seq)) {
    bartC_tmp <- bartc(Y[idx] - tau.seq[k]*Trt[idx], Trt[idx], X0[idx, ], keepTrees = TRUE, n.samples = 100L, n.burn = 15L, n.chains = 2L)  
    
    fitted <- colMeans(apply(bartC_tmp$mu.hat.obs, c(3), function(x) rowMeans(x)))
    mse[k] <- mean(((Y[idx] - tau.seq[k]*Trt[idx]) - fitted)^2)
  }
  tau.star <- tau.seq[which.min(mse)]

  fit_reduced <- bartc(Y[idx] - tau.star*Trt[idx], Trt[idx], X0[idx, ], keepTrees = TRUE, n.samples = 100L, n.burn = 15L, n.chains = 2L)  
  
  return(list(full=fit, reduced=fit_reduced, tau = tau.star))
  
  fit
}

# Function for making predictions using BART Cause model
predictor_bartc <- function(fit, X_new, Treat_new, funcs_params = NA) {
  colMeans(predict(fit, newdata = cbind.data.frame(X_new, z = Treat_new)))
}

# Linear regression functions
bartC_funs <- list(fitter = fitter_bartc,
                   predictor = predictor_bartc,
                   loss = squared_loss,
                   name = "bartC")




n_folds <- 6
nested_cv_reps <- 10 #average over many random splits


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

nested_cv(data.frame(DATA.cor), data.frame(DATA.cor.reduced), as.vector(Y), as.vector(Treat),tau.seq = tau.range, bartC_funs, 
          n_folds = n_folds, reps  = nested_cv_reps, verbose = T)






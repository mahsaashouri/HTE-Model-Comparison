
library(dplyr)
library(bartCause)
######
## full model
#####
# Function for squared loss
squared_loss <- function(y1, y2, funcs_params = NA) {
  (y1 - y2)^2
}

# Function for fitting BART Cause model
fitter_bartc <- function(X, Y, Treat, idx = NA, funcs_params = NA) {
  if (sum(is.na(idx)) > 0) {
    idx <- 1:nrow(X)
  }
  fit <- bartc(Y[idx], Treat[idx], X[idx, ], keepTrees = TRUE, n.samples = 100L, n.burn = 15L, n.chains = 2L)  
  
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

#set.seed(123)
#train_idx <- sample(1:nrow(DATA.cor), round(.7 * nrow(DATA.cor)), replace = FALSE)
#test_idx <- setdiff(1:nrow(DATA.cor), train_idx)

# Create training and test sets using the selected index numbers
Y <- DATA.cor$iqsb.36
Treat <- DATA.cor$treat
DATA.cor <- DATA.cor[ , !(names(DATA.cor) %in% c('iqsb.36', 'treat'))]
DATA.cor <- model.matrix(Y~.-1, data = DATA.cor)

#train.set <- DATA.cor[train_idx, ]
#Y.train <- Y[train_idx]
#test.set <- DATA.cor[test_idx, ]
#Y.test <-  Y[test_idx]

n_folds <- 6
nested_cv_reps <- 5000 #average over many random splits

nested_cv_BART(data.frame(DATA.cor), as.vector(Y), as.vector(Treat), bartC_funs, 
               n_folds = n_folds, reps  = nested_cv_reps, verbose = T, alpha = 0.01)
#nested_cv_helper_BART(data.frame(DATA), as.vector(Y), as.vector(Treat), bartC_funs, 
#                 n_folds = n_folds)
#naive_cv_BART(data.frame(DATA), as.vector(Y), as.vector(Treat),bartC_funs, n_folds = n_folds, alpha = 0.5,
#                          trans = list(identity), funcs_params = NULL, fold_id = NULL)


######
## reduced model
#####

squared_loss <- function(y1, y2, t, tau, funcs_params = NA) {
  (y1 - (y2-(tau*t)))^2
}

fitter_bartc <- function(X, Y, Treat, tau, idx = NA, funcs_params = NA) {
  if (sum(is.na(idx)) > 0) {
    idx <- 1:nrow(X)
  }
  fit <- bartc(Y[idx]-(tau*Treat[idx]), Treat[idx], X[idx, ], keepTrees = TRUE, n.samples = 100L, n.burn = 15L, n.chains = 2L)  
  
  fit
}

predictor_bartc <- function(fit, X_new, Treat_new, funcs_params = NA) {
  colMeans(predict(fit, newdata = cbind.data.frame(X_new, z = Treat_new)))
}

bartC_funs <- list(fitter = fitter_bartc,
                   predictor = predictor_bartc,
                   loss = squared_loss,
                   name = "bartC")

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
#train_idx <- sample(1:nrow(DATA.cor), round(.7 * nrow(DATA.cor)), replace = FALSE)
#test_idx <- setdiff(1:nrow(DATA.cor), train_idx)

# Create training and test sets using the selected index numbers
Y <- DATA.cor$iqsb.36
Treat <- DATA.cor$treat
DATA.cor <- DATA.cor[ , !(names(DATA.cor) %in% c('iqsb.36', 'treat'))]

DATA.cor <- model.matrix(Y~.-1, data = DATA.cor)

#train.set <- DATA.cor[train_idx, ]
#Y.train <- Y[train_idx]
#Treat.train <- treat[train_idx]
#test.set <- DATA.cor[test_idx, ]
#Y.test <-  Y[test_idx]
#Treat.test <- treat[test_idx]

tau.range = seq(1,10, by =1)
nested_cv_m(data.frame(DATA.cor), as.vector(Y), as.vector(Treat), tau.range, bartC_funs, 
            n_folds = n_folds, reps  = nested_cv_reps, verbose = T, alpha = 0.01)


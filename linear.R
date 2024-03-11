
library(dplyr)


######################
## Reduced model
######################

squared_loss <- function(y1, y2, t, tau, funcs_params = NA) {
  (y1 - (y2-(tau*t)))^2
}

fitter_lm <- function(X, Y, Treat, tau, idx = NA, funcs_params = NA) {
  if (sum(is.na(idx)) > 0) {
    idx <- 1:nrow(X)
  }
  data.all <- cbind.data.frame(X[idx, ], Y = Y[idx]-(tau*Treat[idx]))
  fit <- lm(Y ~., data = data.all) 
  
  fit
}

predictor_lm <- function(fit, X_new, funcs_params = NA) {
  predict(fit, newdata = X_new)
}

linear_regression_funs <- list(fitter = fitter_lm,
                               predictor = predictor_lm,
                               loss = squared_loss,
                               name = "linear_regression")


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

tau.range = seq(1,10, by =1)
nested_cv_m(data.frame(train.set), as.vector(Y.train), as.vector(Treat.train), tau.range, linear_regression_funs, 
            n_folds = n_folds, reps  = nested_cv_reps, verbose = T, alpha = 0.5)



######################
## Full model
######################

squared_loss <- function(y1, y2, funcs_params = NA) {
  (y1 - y2)^2
}

fitter_lm <- function(X, Y, idx = NA, funcs_params = NA) {
  if (sum(is.na(idx)) > 0) {
    idx <- 1:nrow(X)
  }
  data.all <- cbind.data.frame(X[idx, ], Y = Y[idx])
  fit <- lm(Y ~., data = data.all)  
  
  fit
}

predictor_lm <- function(fit, X_new, funcs_params = NA) {
  predict(fit, newdata = X_new)
}

linear_regression_funs <- list(fitter = fitter_lm,
                               predictor = predictor_lm,
                               loss = squared_loss,
                               name = "linear_regression")


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

DATA.cor <- bind_cols(output, selected_data, selected_data.T)
DATA.cor <- bind_cols('iqsb.36' = output, selected_data)

set.seed(123)
#train_idx <- sample(1:nrow(DATA.cor), round(.7 * nrow(DATA.cor)), replace = FALSE)
#test_idx <- setdiff(1:nrow(DATA.cor), train_idx)

# Create training and test sets using the selected index numbers
Y <- DATA.cor$iqsb.36
DATA.cor <- DATA.cor[ , !(names(DATA.cor) %in% c('iqsb.36'))]
DATA.cor <- model.matrix(Y~.-1, data = DATA.cor)

#train.set <- DATA.cor[train_idx, ]
#Y.train <- Y[train_idx]
#test.set <- DATA.cor[test_idx, ]
#Y.test <-  Y[test_idx]


nested_cv(data.frame(DATA.cor), as.vector(Y), linear_regression_funs, 
                    n_folds = n_folds, reps  = nested_cv_reps, verbose = T, alpha = 0.5)


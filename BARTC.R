
# Function for squared loss
squared_loss <- function(y1, y2, funcs_params = NA) {
  (y1 - y2)^2
}

# Function for fitting linear regression model
fitter_bartc <- function(X, Y, TR, idx = NA, funcs_params = NA) {
  if (sum(is.na(idx)) > 0) {
    idx <- 1:nrow(X)
  }
  fit <- bartc(Y[idx], TR[idx], X[idx, ], keepTrees = TRUE)  # Assumes X and Y are matrices/data.frames
  
  fit
}

# Function for making predictions using linear regression model
predictor_bartc <- function(fit, X_new, TR_new, funcs_params = NA) {
  colMeans(predict(fit, newdata = cbind.data.frame(X_new, TR_new), type = "mu"))
}

# Linear regression functions
bartC_funs <- list(fitter = fitter_bartc,
                               predictor = predictor_bartc,
                               loss = squared_loss,
                               name = "bartC")


n_folds <- 6
nested_cv_reps <- 300 #average over many random splits


DATA <- read.csv('IHDP_clean.csv', header = TRUE)[,-1]

train_idx <- sample(1:nrow(DATA), round(.7 * nrow(DATA)), replace = FALSE)
test_idx <- setdiff(1:nrow(DATA), train_idx)

# Create training and test sets using the selected index numbers
Y <- DATA$iqsb.36
TR <- DATA$treat
DATA <- DATA[ , !(names(DATA) %in% c('iqsb.36', 'treat'))]
DATA <- model.matrix(Y~.-1, data = DATA)

train.set <- DATA[train_idx, ]
Y.train <- Y[train_idx]
TR.train <- TR[train_idx]
test.set <- DATA[test_idx, ]
Y.test <-  Y[test_idx]
TR.test <- TR[test_idx]

library(bartCause)
fit <- bartc(as.vector(Y.train), as.vector(TR.train), data.frame(train.set), keepTrees = TRUE)
predict(fit, newdata = cbind.data.frame(test.set, z = TR.test), type = "mu")


nested_cv(data.frame(train.set), as.vector(Y.train), as.vector(TR.train), bartC_funs, 
          n_folds = n_folds, reps  = nested_cv_reps, verbose = T)
#nested_cv_helper(data.frame(train.set), as.vector(Y.train), bartC_funs, 
 #                n_folds = n_folds)


library(mboost)

# Function for squared loss
squared_loss <- function(y1, y2, funcs_params = NA) {
  (y1 - y2)^2
}

# Function for fitting linear regression model
fitter_glmboost <- function(X, Y, idx = NA, funcs_params = NA) {
  if (sum(is.na(idx)) > 0) {
    idx <- 1:nrow(X)
  }
  fit <-  glmboost(x= as.matrix(X[idx, ]), y=Y[idx],family = Gaussian())
  fit
}

# Function for making predictions using linear regression model
predictor_glmboost <- function(fit, X_new, funcs_params = NA) {
 as.numeric(predict(fit, newdata = as.matrix(X_new)))
}

# Linear regression functions
glmboost_funs <- list(fitter = fitter_glmboost,
                   predictor = predictor_glmboost,
                   loss = squared_loss,
                   name = "glmboost")


n_folds <- 6
nested_cv_reps <- 300 #average over many random splits


DATA <- read.csv('IHDP_clean.csv', header = TRUE)[,-1]

train_idx <- sample(1:nrow(DATA), round(.7 * nrow(DATA)), replace = FALSE)
test_idx <- setdiff(1:nrow(DATA), train_idx)

# Create training and test sets using the selected index numbers
Y <- DATA$iqsb.36
DATA <- DATA[ , !(names(DATA) %in% c('iqsb.36'))]
DATA <- model.matrix(Y~.-1, data = DATA)

train.set <- DATA[train_idx, ]
Y.train <- Y[train_idx]
test.set <- DATA[test_idx, ]
Y.test <-  Y[test_idx]


nested_cv(data.frame(train.set), as.vector(Y.train), glmboost_funs, 
          n_folds = n_folds, reps  = nested_cv_reps, verbose = T)
nested_cv_helper(data.frame(train.set), as.vector(Y.train), glmboost_funs, 
                n_folds = n_folds)

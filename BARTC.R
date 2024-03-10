
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

#train_idx <- sample(1:nrow(DATA), round(.7 * nrow(DATA)), replace = FALSE)
#test_idx <- setdiff(1:nrow(DATA), train_idx)

# Create training and test sets using the selected index numbers
Y <- DATA$iqsb.36
Treat <- DATA$treat
DATA <- DATA[ , !(names(DATA) %in% c('iqsb.36', 'treat'))]
DATA <- model.matrix(Y~.-1, data = DATA)

library(bartCause)

## if we want to test how to run the bartCause
#train.set <- DATA[train_idx, ]
#Y.train <- Y[train_idx]
#Treat.train <- Treat[train_idx]
#test.set <- DATA[test_idx, ]
#Y.test <-  Y[test_idx]
#Treat.test <- Treat[test_idx]
#fit <- bartc(as.vector(Y.train), as.vector(Treat.train), data.frame(train.set), keepTrees = TRUE)
#pred <- predict(fit, newdata = cbind.data.frame(test.set, z = Treat.test))


alpha <- 0.01
n_folds <- 6
nested_cv_reps <- 5000 #average over many random splits

nested_cv_BART(data.frame(DATA), as.vector(Y), as.vector(Treat), bartC_funs, 
          n_folds = n_folds, reps  = nested_cv_reps, verbose = T)
#nested_cv_helper_BART(data.frame(DATA), as.vector(Y), as.vector(Treat), bartC_funs, 
#                 n_folds = n_folds)
#naive_cv_BART(data.frame(DATA), as.vector(Y), as.vector(Treat),bartC_funs, n_folds = n_folds, alpha = 0.5,
#                          trans = list(identity), funcs_params = NULL, fold_id = NULL)

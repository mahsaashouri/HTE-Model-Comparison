

# Load necessary libraries
library(glmnet)
library(boot)
ihdp_clean <- read.csv('IHDP_clean.csv')[,-1]

# Create vector of random indices
indices <- sample(1:nrow(ihdp_clean))
# Calculate number of observations in training set
n_train <- floor(0.8 * nrow(ihdp_clean))
# Subset data into training and testing sets
train <- ihdp_clean[indices[1:n_train], ]
test <- ihdp_clean[indices[(n_train+1):nrow(ihdp_clean)], ]

# Fit linear model and perform cross-validation
linear_model <- glm(iqsb.36 ~ ., data  = train)
## cv.glm: This function calculates the estimated K-fold cross-validation prediction error for generalized linear models.
linear_cv <- cv.glm(train, linear_model)
#print(paste("Linear model RMSE:", round(sqrt(mean((linear_cv$seed)^2)), 2)))

# Fit glmnet model and perform cross-validation
glmnet_model <- glmnet(as.matrix(train[, -which(names(train) == "iqsb.36")]), 
                       train$iqsb.36, alpha = 1)
glmnet_cv <- cv.glmnet(as.matrix(train[, -which(names(train) == "iqsb.36")]), 
                       train$iqsb.36, alpha = 1)
#print(paste("GLMNET model RMSE:", round(sqrt(mean(glmnet_cv$cvm)^2), 2)))
## bartCause package
library(bartCause)
model <- bartc(treatment = train$treat, response  = train$iqsb.36, confounders = train[, -c(1, 2)], num_trees = 2000, n_burn = 1000, n_post = 1000,
               prior_t = "auto", prior_int = "auto", method_ps = "glm")

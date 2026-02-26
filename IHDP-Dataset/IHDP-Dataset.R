library(glmnet)
library(mboost)
library(dbarts)
library(tidyverse)
library(randomForest)
# Set seed for reproducibility
set.seed(12356)
##################################
## IHDP dataset
DATA <- read.csv('IHDP_clean.csv', header = TRUE)[,-1]
cor_matrix <- cor(DATA)
cor_with_output <- cor_matrix[, 1]
cor_with_output[1] <- 0  # Set the correlation with output to 0
sorted_correlations <- sort(abs(cor_with_output), decreasing = TRUE)
top_10_correlations <- sorted_correlations[1:10]
top_10_predictors <- c(names(top_10_correlations), "momage")
# Subset your data to keep only the top 10 predictor columns
selected_data <- DATA[, top_10_predictors]
## reduced model
output <- DATA[,c('iqsb.36', 'treat')]
DAT_reduced <- bind_cols(output, selected_data)
## full model
selected_data.T <- matrix(NA, nrow = nrow(selected_data), ncol = ncol(selected_data))
for(k in 1:11){
  selected_data.T[,k] <- selected_data[,k]*output$treat
}
colnames(selected_data.T) <- paste0(colnames(selected_data), ".t")
DAT <- bind_cols(output, selected_data, selected_data.T)

drop_cols <- c(which(colnames(DAT_reduced)=="iqsb.36"), which(colnames(DAT_reduced)=="treat"))
DAT_red <- data.frame(Wtau=DAT_reduced$iqsb.36, selected_data)

## Source fitting and nested cv functions
setwd("~/Projects/HTE-Model-Comparison")  ## Change for your computer
source("CoreFunctions/CVDiffFunctions.R")
source("CoreFunctions/FitterFunctions.R")

## Define linear regression, glmnet, and boosting fitting and prediction functions.
linear_regression_funs <- list(fitter = fitter_lm,
                               predictor = predictor_lm,
                               mse = mse,
                               loss = squared_loss,
                               name = "linear_regression")

glmnet_funs <- list(fitter = fitter_glmnet,
                    predictor = predictor_glmnet,
                    mse = mse,
                    loss = squared_loss,
                    name = "glmnet")

glmboost_funs <- list(fitter = fitter_glmboost,
                      predictor = predictor_glmboost,
                      mse = mse,
                      loss = squared_loss,
                      name = "glmboost")

rf_funs <- list(fitter = fitter_rf,
                predictor = predictor_rf,
                mse = mse,
                loss = squared_loss,
                name = "rf")

# Set the number of folds, and
# number of nested cv replications:
n_folds <- 5
nested_cv_reps <- 10 ## Use 50 or 100 for paper

## Set the number of simulation replications
nreps <- 500  ## Use nreps = 500 for paper
cover_lm <- cover_glmnet <- cover_glmboost <- cover_rf <- rep(NA, nreps)
hvalue_lm <- hvalue_glmnet <- hvalue_glmboost <- hvalue_rf <- rep(NA, nreps)
CI_lm <- CI_glmnet <- CI_glmboost <- CI_rf <- matrix(NA, nrow=nreps, ncol=2)
true_thetas <- matrix(NA, nreps, 4)
for(h in 1:nreps) {
  
  ######################
  ## Compute true value of \theta_{XY} for each method
  ######################

  Xmat_tmp <- model.matrix(iqsb.36 ~ . - 1, data=DAT)
  X0mat_tmp <- model.matrix(iqsb.36 ~ . - 1, data=DAT_reduced)
  ## Don't use treatment variable in GLMboost
  X0mat_notrt_tmp <- model.matrix(iqsb.36 ~ .- treat - 1, data=DAT_reduced)
  
  tmp_lm <- lm(iqsb.36 ~., data=DAT)
  tmp_reduced_lm <- lm(iqsb.36 ~., data=DAT_reduced)
  tmp_glmnet <- cv.glmnet(Xmat_tmp, DAT$iqsb.36, family = "gaussian", nfolds = 5)
  tmp_reduced_glmnet <- cv.glmnet(X0mat_tmp, DAT$iqsb.36, family = "gaussian", nfolds = 5)
  
  tmp_glmboost <- glmboost(x=Xmat_tmp, y=DAT$iqsb.36, family = Gaussian())
  
  f0fn <- function(tau) {
    Wtau <- DAT_reduced$iqsb.36 - tau*DAT_reduced$treat
    fit_reduced <- glmboost(x=X0mat_notrt_tmp, y=Wtau,family = Gaussian())
    mse_local <- mean((Wtau - predict(fit_reduced, newdata=X0mat_notrt_tmp))^2)
    return(mse_local)
  }
  tau.star.gboost <- optimize(f0fn, interval=c(-5, 5))$minimum
  Wtau.star <- DAT_reduced$iqsb.36 - tau.star.gboost*DAT_reduced$treat
  tmp_reduced_glmboost <- glmboost(x=X0mat_notrt_tmp, y = Wtau.star,family = Gaussian())
  
  
  f00fn <- function(tau) {
    DAT_red$Wtau <- DAT_reduced$iqsb.36 - tau*DAT_reduced$treat
    
    fit_reduced <- randomForest(Wtau ~., data=DAT_red, maxnodes=16, ntree=100)
    mse_local <- mean((DAT_red$Wtau - predict(fit_reduced, newdata=DAT_red))^2)
    return(mse_local)
  }
  tau.star.rf <- optimize(f00fn, interval=c(-5, 5))$minimum
  DAT_red$Wtau <- DAT_reduced$iqsb.36 - tau.star.rf*DAT_reduced$treat
  ## Is the line below correct?
  tmp_reduced_rf <- randomForest(Wtau ~., data=DAT_red, maxnodes=16, ntree=100)

  
  #######################################################
  ###. Getting confidence intervals and h-values
  ######################################################
  naive.trt.effect <- mean(DAT$iqsb.36[DAT$treat==1]) - mean(DAT$iqsb.36[DAT$treat==0])
  tau.range <- sort(c(-5*naive.trt.effect, 5*naive.trt.effect))
  XX <- data.frame(model.matrix(iqsb.36 ~. - 1, data = DAT))
  XX0 <- data.frame(model.matrix(iqsb.36 ~. - treat - 1, data = DAT_reduced))
  
  ncv_lm <- nested_cv(X=XX, X0=XX0, Y=as.vector(DAT$iqsb.36), Trt=DAT$treat, tau.range=tau.range, funcs=linear_regression_funs,
                      n_folds = n_folds, reps  = nested_cv_reps)
  
  ncv_net <- nested_cv(X=XX, X0=XX0, Y=as.vector(DAT$iqsb.36), Trt=DAT$treat, tau.range=tau.range, funcs=glmnet_funs,
                       n_folds = n_folds, reps  = nested_cv_reps)
  
  ncv_boost <- nested_cv(X=XX, X0=XX0, Y=as.vector(DAT$iqsb.36), Trt=DAT$treat, tau.range=tau.range, funcs=glmboost_funs,
                         n_folds = n_folds, reps  = nested_cv_reps)
  
  ncv_rf <- nested_cv(X=XX, X0=XX0, Y=as.vector(DAT$iqsb.36), Trt=DAT$treat, tau.range=tau.range, funcs=rf_funs,
                      n_folds = n_folds, reps = 10, bias_reps = 0)
  #ncv_bart <- nested_cv(X=XX, X0=XX0, Y=as.vector(DAT$iqsb.36), Trt=DAT$treat, tau.range=tau.range, funcs=bart_funs,
  #                       n_folds = n_folds, reps  = nested_cv_reps)
  ############################################
  ### Record Results
  ############################################
  #cover_lm[h] <- theta_lm > ncv_lm$ci_lo & theta_lm < ncv_lm$ci_hi
  CI_lm[h,1] <- ncv_lm$ci_lo
  CI_lm[h,2] <- ncv_lm$ci_hi
  #true_thetas[h,] <- c(theta_lm, theta_glmnet, theta_glmboost, theta_bart)
  hvalue_lm[h] <- ncv_lm$hvalue
  
  #cover_glmnet[h] <- theta_glmnet > ncv_net$ci_lo & theta_glmnet < ncv_net$ci_hi
  CI_glmnet[h,1] <- ncv_net$ci_lo
  CI_glmnet[h,2] <- ncv_net$ci_hi
  hvalue_glmnet[h] <- ncv_net$hvalue
  
  #cover_glmboost[h] <- theta_glmboost > ncv_boost$ci_lo & theta_glmboost < ncv_boost$ci_hi
  CI_glmboost[h,1] <- ncv_boost$ci_lo
  CI_glmboost[h,2] <- ncv_boost$ci_hi
  hvalue_glmboost[h] <- ncv_boost$hvalue
  
  #cover_rf[h] <- theta_rf > ncv_rf$ci_lo & theta_rf < ncv_rf$ci_hi
  CI_rf[h,1] <- ncv_rf$ci_lo
  CI_rf[h,2] <- ncv_rf$ci_hi
  hvalue_rf[h] <- ncv_rf$hvalue
  cat("Simulation Replication: ", h, "\n")
}
## Save: cover_lm, cover_glmnet, cover_glmboost, cover_bart
##       CI_lm, CI_glmnet, CI_glmboost, CI_bart
##      hvalue_lm, hvalue_glmnet, hvalue_glmboost, hvalue_bart


# histogram - h-values
#data_hvalues <- data.frame(
#  value = c(hvalue_lm, hvalue_glmnet, hvalue_glmboost, hvalue_rf),
#  group = factor(rep(c("linear", "Lasso", "boosting", "random forest"), each = 50), levels = c("linear", "Lasso", "boosting", "random forest"))
#)

data_hvalues <- data.frame(
  value = c(hvalue_lm, hvalue_glmnet, hvalue_glmboost),
  group = factor(rep(c("linear", "Lasso", "boosting"), each = 500), levels = c("linear", "Lasso", "boosting"))
)

# Plot the overlapping histograms
ggplot(data_hvalues, aes(x = value, fill = group)) +
  geom_histogram(position = "identity", alpha = 0.7, color = 'darkgray', binwidth = 0.05) +
  labs(x = "h-value", y = "Frequency") +
  #scale_fill_manual(values = c("lightblue", "lightpink3", "darkolivegreen4", "darkgoldenrod"), name = "Methods") +
  scale_fill_manual(values = c("lightblue", "lightpink3", "darkolivegreen4"), name = "Methods") +
  coord_cartesian(xlim = c(0,1)) +
  theme_minimal()+
  theme(
    text = element_text(size = 20),  
    axis.title = element_text(size = 22), 
    axis.text = element_text(size = 20),
    plot.title = element_text(size = 16, face = "bold"), 
    legend.title = element_text(size = 24), 
    legend.text = element_text(size = 20),
    legend.position = "bottom" 
  )


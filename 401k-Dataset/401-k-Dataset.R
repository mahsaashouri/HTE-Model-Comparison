
library(glmnet)
library(mboost)
library(dbarts)
library(tidyverse)
library(randomForest)

# Set seed for reproducibility
set.seed(12356)

##################################
## SOLVDT dataset
#install.packages("DoubleML")
library(DoubleML)
DATA = fetch_401k(return_type = "data.table")
head(DATA)

##################################

## reduced model
Cov_DAT <- DATA[, -c('net_tfa', 'e401')]
output <- DATA[, c('net_tfa', 'e401')]
DAT_reduced <- bind_cols(output, Cov_DAT)
#DAT_reduced <- bind_cols(output, selected_data)

## full model - create interaction terms

Cov_DAT.T <- matrix(NA, nrow = nrow(Cov_DAT), ncol = ncol(Cov_DAT))
for(k in 1:ncol(Cov_DAT)){
  Cov_DAT.T[,k] <- as.matrix(Cov_DAT)[,k] * output$e401
}
colnames(Cov_DAT.T) <- paste0(colnames(Cov_DAT), ".t")
DAT <- bind_cols(output, Cov_DAT, Cov_DAT.T)

# Create data for random forest
drop_cols <- c(which(colnames(DAT_reduced)=="net_tfa"), which(colnames(DAT_reduced)=="e401"))
DAT_red <- data.frame(Wtau=DAT_reduced$net_tfa, Cov_DAT)

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
nested_cv_reps <- 500 ## Use 50 or 100 for paper

## Set the number of simulation replications
nreps <- 5  ## Use nreps = 50 for paper
cover_lm <- cover_glmnet <- cover_glmboost <- cover_rf <- rep(NA, nreps)
hvalue_lm <- hvalue_glmnet <- hvalue_glmboost <- hvalue_rf <- rep(NA, nreps)
CI_lm <- CI_glmnet <- CI_glmboost <- CI_rf <- matrix(NA, nrow=nreps, ncol=2)
true_thetas <- matrix(NA, nreps, 4)

for(h in 1:nreps) {
  
  ######################
  ## Compute true value of \theta_{XY} for each method
  ######################
  
  Xmat_tmp <- model.matrix(net_tfa ~ . - 1, data=DAT)
  X0mat_tmp <- model.matrix(net_tfa ~ .  - 1, data=DAT_reduced)
  ## Don't use treatment variable in GLMboost
  X0mat_notrt_tmp <- model.matrix(net_tfa ~ .- e401 - 1, data=DAT_reduced)
  
  tmp_lm <- lm(net_tfa ~., data=DAT)
  tmp_reduced_lm <- lm(net_tfa ~., data=DAT_reduced)
  tmp_glmnet <- cv.glmnet(Xmat_tmp, DAT$net_tfa, family = "gaussian", nfolds = 5)
  tmp_reduced_glmnet <- cv.glmnet(X0mat_tmp, DAT$net_tfa, family = "gaussian", nfolds = 5)
  
  tmp_glmboost <- glmboost(x=Xmat_tmp, y=DAT$net_tfa, family = Gaussian())
  
  f0fn <- function(tau) {
    Wtau <- DAT_reduced$net_tfa - tau*DAT_reduced$e401
    fit_reduced <- glmboost(x=X0mat_notrt_tmp, y=Wtau, family = Gaussian())
    mse_local <- mean((Wtau - predict(fit_reduced, newdata=X0mat_notrt_tmp))^2)
    return(mse_local)
  }
  tau.star.gboost <- optimize(f0fn, interval=c(-5, 5))$minimum
  Wtau.star <- DAT_reduced$net_tfa - tau.star.gboost*DAT_reduced$e401
  tmp_reduced_glmboost <- glmboost(x=X0mat_notrt_tmp, y = Wtau.star, family = Gaussian())
  
  
  f00fn <- function(tau) {
    DAT_red$Wtau <- DAT_reduced$net_tfa - tau*DAT_reduced$e401
    
    fit_reduced <- randomForest(Wtau ~., data=DAT_red, maxnodes=16, ntree=100)
    mse_local <- mean((DAT_red$Wtau - predict(fit_reduced, newdata=DAT_red))^2)
    return(mse_local)
  }
  tau.star.rf <- optimize(f00fn, interval=c(-5, 5))$minimum
  DAT_red$Wtau <- DAT_reduced$net_tfa - tau.star.rf*DAT_reduced$e401
  ## Is the line below correct?
  tmp_reduced_rf <- randomForest(Wtau ~., data=DAT_red, maxnodes=16, ntree=100)
  
  
  #######################################################
  ###. Getting confidence intervals and h-values
  ######################################################
  naive.trt.effect <- mean(DAT$net_tfa[DAT$e401==1]) - mean(DAT$net_tfa[DAT$e401==0])
  tau.range <- sort(c(-5*naive.trt.effect, 5*naive.trt.effect))
  XX <- data.frame(model.matrix(net_tfa ~. - 1, data = DAT))
  XX0 <- data.frame(model.matrix(net_tfa ~. - e401 - 1, data = DAT_reduced))
  
  ncv_lm <- nested_cv(X=XX, X0=XX0, Y=as.vector(DAT$net_tfa), Trt=DAT$e401, tau.range=tau.range, funcs=linear_regression_funs,
                      n_folds = n_folds, reps  = nested_cv_reps)
  
  ncv_net <- nested_cv(X=XX, X0=XX0, Y=as.vector(DAT$net_tfa), Trt=DAT$e401, tau.range=tau.range, funcs=glmnet_funs,
                       n_folds = n_folds, reps  = nested_cv_reps)
  
  ncv_boost <- nested_cv(X=XX, X0=XX0, Y=as.vector(DAT$net_tfa), Trt=DAT$e401, tau.range=tau.range, funcs=glmboost_funs,
                         n_folds = n_folds, reps  = nested_cv_reps)
  
  ncv_rf <- nested_cv(X=XX, X0=XX0, Y=as.vector(DAT$net_tfa), Trt=DAT$e401, tau.range=tau.range, funcs=rf_funs,
                      n_folds = n_folds, reps = 10, bias_reps = 0)
  
  ############################################
  ### Record Results
  ############################################
  CI_lm[h,1] <- ncv_lm$ci_lo
  CI_lm[h,2] <- ncv_lm$ci_hi
  hvalue_lm[h] <- ncv_lm$hvalue
  
  CI_glmnet[h,1] <- ncv_net$ci_lo
  CI_glmnet[h,2] <- ncv_net$ci_hi
  hvalue_glmnet[h] <- ncv_net$hvalue
  
  CI_glmboost[h,1] <- ncv_boost$ci_lo
  CI_glmboost[h,2] <- ncv_boost$ci_hi
  hvalue_glmboost[h] <- ncv_boost$hvalue
  
  CI_rf[h,1] <- ncv_rf$ci_lo
  CI_rf[h,2] <- ncv_rf$ci_hi
  hvalue_rf[h] <- ncv_rf$hvalue
  
  cat("Simulation Replication: ", h, "\n")
}

# Create histogram data for h-values
data_hvalues <- data.frame(
  value = c(hvalue_lm, hvalue_glmnet, hvalue_glmboost),
  group = factor(rep(c("linear", "Lasso", "boosting"), each = nreps), levels = c("linear", "Lasso", "boosting"))
)

# Plot the overlapping histograms
ggplot(data_hvalues, aes(x = value, fill = group)) +
  geom_histogram(position = "identity", alpha = 0.7, color = 'darkgray', binwidth = 0.03) +
  labs(x = "h-value", y = "Frequency") +
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



save(
  cover_lm, cover_glmnet, cover_glmboost, cover_rf,
  CI_lm, CI_glmnet, CI_glmboost, CI_rf,
  hvalue_lm, hvalue_glmnet, hvalue_glmboost, hvalue_rf,
  true_thetas,
  file = "401k_nrep5.RData"
)

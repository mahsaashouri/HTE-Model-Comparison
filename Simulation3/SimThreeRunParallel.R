library(parallel)
library(glmnet)
library(mboost)
library(randomForest)

# Set seed for reproducibility
set.seed(123567)

## Source fitting and nested cv functions
#setwd("/Users/mahsa/Projects/HTE-Model-Comparison")  ## Change for your computer

setwd("~/Documents/HTEevaluation/HTE-Model-Comparison")
source("CoreFunctions/CVDiffFunctionsTransformed.R")
source("CoreFunctions/FitterFunctionsTransformed.R")
source("Simulation3/OneSim3Run.R")


# Set a seed for reproducibility
set.seed(123)

# Define the number of cores to use
num_cores <- detectCores() - 1  # Use all cores except one

## change muvec and thetavec if you don't want to
## run as many simulation settings
#muvec <- rep(1:3, each=4)
#thetavec <- rep(1:4, 3)
muvec <- rep(c(1,2,3), each=4)
thetavec <- rep(1:4, 3)

# Set up the cluster
for(j in 1:length(muvec)) {
  cl <- makeCluster(num_cores)

  # Export the necessary objects to the cluster
  clusterExport(cl, c("OneSimThreeRun", "mu", "theta", "f1", "f2", "f3", "f4", "cv.glmnet", "glmboost",
                      "predict.glmnet", "predict.glmboost", "Gaussian",
                      "fitter_lm", "fitter_glmnet", "fitter_glmboost", "fitter_ridge", "wsquared_loss",
                      "predictor_lm", "predictor_glmnet", "predictor_glmboost", "predictor_ridge",
                      "nested_cv", "nested_cv_helper", "naive_cv", "linear_regression_funs",
                      "glmnet_funs", "glmboost_funs", "ridge_funs", "muvec", "thetavec", "j"))

  # Run the simulation
  ## first argument in rep(500, 500) is sample size
  ## second argument in rep(500, 500) is number of simulation replications
  results <- parLapply(cl, rep(5000, 500), function(n) { OneSimThreeRun(n, mu_choice=muvec[j], theta_choice=thetavec[j], p=9, n_folds=5, nested_cv_reps=10) } )

  cover_lm <- cover_glmnet <- cover_glmboost <- cover_ridge <- rep(NA, length(results))
  hvalue_lm <- hvalue_glmnet <- hvalue_glmboost <- hvalue_ridge <- rep(NA, length(results))
  hvalue1_lm <- hvalue1_glmnet <- hvalue1_glmboost <- hvalue1_ridge <- rep(NA, length(results))
  CI_lm <- CI_glmnet <- CI_glmboost <- CI_ridge <- matrix(NA, nrow=length(results), ncol=2)
  TrueThetas <- matrix(NA, nrow=length(results), ncol=4)
  for(k in 1:length(results)) {
    cover_lm[k] <- results[[k]]$cover_lm
    cover_glmnet[k] <- results[[k]]$cover_glmnet
    cover_glmboost[k] <- results[[k]]$cover_glmboost
    cover_ridge[k] <- results[[k]]$cover_ridge

    CI_lm[k,] <- results[[k]]$CI_lm
    CI_glmnet[k,] <- results[[k]]$CI_glmnet
    CI_glmboost[k,] <- results[[k]]$CI_glmboost
    CI_ridge[k,] <- results[[k]]$CI_ridge

    hvalue_lm[k] <- results[[k]]$hvalue_lm
    hvalue1_lm[k] <- results[[k]]$hvalue1side_lm
    hvalue_glmnet[k] <- results[[k]]$hvalue_glmnet
    hvalue1_glmnet[k] <- results[[k]]$hvalue1side_glmnet
    hvalue_glmboost[k] <- results[[k]]$hvalue_glmboost
    hvalue1_glmboost[k] <- results[[k]]$hvalue1side_glmboost
    hvalue_ridge[k] <- results[[k]]$hvalue_ridge
    hvalue1_ridge[k] <- results[[k]]$hvalue1side_ridge

    TrueThetas[k,] <- results[[k]]$true_thetas
  }
  ## Change file name here:
  fname <- paste("~/Documents/HTEevaluation/HTE-Model-Comparison/SimulationResults/SimulationThreeResults", muvec[j],"Theta",thetavec[j],"n5000.RData", sep="")
  save(cover_lm, cover_glmnet, cover_glmboost, cover_ridge, CI_lm, CI_glmnet,
       CI_glmboost, CI_ridge, hvalue_lm, hvalue_glmnet, hvalue_glmboost, hvalue_ridge,
       hvalue1_lm, hvalue1_glmnet, hvalue1_glmboost, hvalue1_ridge, TrueThetas,
       file=fname)
  # Stop the cluster
  stopCluster(cl)
}


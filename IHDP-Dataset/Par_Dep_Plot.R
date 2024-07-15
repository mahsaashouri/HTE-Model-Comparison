##########################################
## Partial Dependence plots
###########################################
library(mboost)
## IHDP dataset
DATA <- read.csv('IHDP_clean.csv', header = TRUE)[,-1]
cor_matrix <- cor(DATA)
cor_with_output <- cor_matrix[, 1]


Col_ParDep <- c('ppvt.imp', 'momwhite', 'momed4F', 'mom.scoll', 'momblack', 'mom.lths', 'b.marry', 'birth.o', 'parity', 'site4')
Y <-  DATA[,c('iqsb.36')]
# Subset your data to keep only the top 10 predictor columns
selected_data <- DATA[, Col_ParDep]
## reduced model
output <- DATA[,c('iqsb.36', 'treat')]
DAT_reduced <- bind_cols(output, selected_data)
## full model
selected_data.T <- matrix(NA, nrow = nrow(selected_data), ncol = ncol(selected_data))
for(k in 1:10){
  selected_data.T[,k] <- selected_data[,k]*output$treat
}
colnames(selected_data.T) <- paste0(colnames(selected_data), ".t")
DAT <- bind_cols(output, selected_data, selected_data.T)

fit_boost_func <- function(data_full, data_reduced){
  Xmat_tmp <- model.matrix(iqsb.36 ~ . - 1, data=data_full)
  X0mat_notrt_tmp <- model.matrix(iqsb.36 ~ . - treat - 1, data=data_reduced)
  tmp_glmboost <- glmboost(x=Xmat_tmp, y=data_full$iqsb.36, family = Gaussian())
  f0fn <- function(tau) {
    Wtau <- data_reduced$iqsb.36 - tau*data_reduced$treat
    fit_reduced <- glmboost(x=X0mat_notrt_tmp, y=Wtau,family = Gaussian())
    mse_local <- mean((Wtau - predict(fit_reduced, newdata=X0mat_notrt_tmp))^2)
    return(mse_local)
  }
  tau.star.gboost <- optimize(f0fn, interval=c(-5, 5))$minimum
  Wtau.star <- data_reduced$iqsb.36 - tau.star.gboost*data_reduced$treat
  tmp_reduced_glmboost <- glmboost(x=X0mat_notrt_tmp, y = Wtau.star,family = Gaussian())
  return(list(tmp_glmboost, tmp_reduced_glmboost, tau.star.gboost))
}

fit_func <- fit_boost_func(DAT, DAT_reduced)

## Setup 

ngrid <- 40
GridEnd <-  rbind(c(38, 160), c(0, 1), c(1, 4), c(0, 1), c(0, 1), c(0, 1), c(1, 4), c(1, 8), c(1, 8), c(0, 1))
ff_full <- matrix(NA, nrow = ngrid, ncol = 3)
ff_reduced <- matrix(NA, nrow = ngrid, ncol = 3)
partial_results_full <- partial_results_reduced <- list()
for(i in 1:length(Col_ParDep)){
  pp <- seq(GridEnd[i,1], GridEnd[i,2], length.out=ngrid)
  for(k in 1:ngrid){
    
    #Xmat_full <- bind_cols(selected_data, selected_data.T)
    Xmat_full <- model.matrix(iqsb.36 ~ .-1, data = DAT)
    Xmat_full[,Col_ParDep[i]] <- rep(pp[k], nrow(DAT))
    
    #Xmat_reduced <- selected_data
    Xmat_reduced <- model.matrix(iqsb.36 ~ .-treat-1, data = DAT_reduced)
    Xmat_reduced[,Col_ParDep[i]] <- rep(pp[k], nrow(DAT_reduced))
    
   full_glmboost <- as.numeric(predict(fit_func[[1]], newdata = Xmat_full))
   reduced_glmboost <- as.numeric(predict(fit_func[[2]], newdata=Xmat_reduced)) + fit_func[[3]]*DAT_reduced$treat
    
    
    ff_full[k,1] <- mean(full_glmboost)
    ff_full[k,2] <- pp[k]
    ff_full[k,3] <- Col_ParDep[i]
    ff_reduced[k,1] <- mean(reduced_glmboost)
    ff_reduced[k,2] <- pp[k]
    ff_reduced[k,3] <- Col_ParDep[i]
    print(c(i, k))
  }
  partial_results_full[[length(partial_results_full)+1]] <- ff_full
  partial_results_reduced[[length(partial_results_reduced)+1]] <- ff_reduced
}


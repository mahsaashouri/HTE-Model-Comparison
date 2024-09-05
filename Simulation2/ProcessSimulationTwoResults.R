
muvec <- rep(1:3, each=4)
thetavec <- rep(1:4, 3)

LMTab <- LassoTab <- BoostTab <- RFTab <- matrix(NA, nrow=length(muvec), ncol=7)
LMTab[,1:2] <- cbind(muvec, thetavec)
LassoTab[,1:2] <- cbind(muvec, thetavec)
BoostTab[,1:2] <- cbind(muvec, thetavec)
RFTab[,1:2] <- cbind(muvec, thetavec)
for(j in 1:length(muvec)) {
  fname <- paste("~/Documents/HTEevaluation/HTE-Model-Comparison/SimulationResults", muvec[j],"Theta",thetavec[j],".RData", sep="")
  load(fname)

  LMTab[j,3] <- mean(TrueThetas[,1])
  LMTab[j,4] <- mean(cover_lm)
  LMTab[j,5] <- mean(CI_lm[,2] - CI_lm[,1])
  LMTab[j,6] <- median(hvalue_lm)
  LMTab[j,7] <- median(hvalue1_lm)
  
  LassoTab[j,3] <- mean(TrueThetas[,2])
  LassoTab[j,4] <- mean(cover_glmnet)
  LassoTab[j,5] <- mean(CI_glmnet[,2] - CI_glmnet[,1])
  LassoTab[j,6] <- median(hvalue_glmnet)
  LassoTab[j,7] <- median(hvalue1_glmnet)
  
  BoostTab[j,3] <- mean(TrueThetas[,3])
  BoostTab[j,4] <- mean(cover_glmboost)
  BoostTab[j,5] <- mean(CI_glmboost[,2] - CI_glmboost[,1])
  BoostTab[j,6] <- median(hvalue_glmboost)
  BoostTab[j,7] <- median(hvalue1_glmboost)
  ## mean(TrueThetas[,3] < 0 & hvalue_glmboost < 0.05)/mean(TrueThetas[,3] < 0)


  RFTab[j,3] <- mean(TrueThetas[,4])
  RFTab[j,4] <- mean(cover_rf)
  RFTab[j,5] <- mean(CI_rf[,2] - CI_rf[,1])
  RFTab[j,6] <- median(hvalue_rf)
  RFTab[j,7] <- median(hvalue1_rf)

}


LMTab
LassoTab
BoostTab
RFTab









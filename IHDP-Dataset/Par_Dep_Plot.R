library(glmnet)
library(gbm)
library(tidyverse)

# Set seed for reproducibility
set.seed(12356)
##################################
## IHDP dataset
DATA <- read.csv('IHDP_clean.csv', header = TRUE)[,-1]

cor_matrix <- cor(DATA)
cor_with_output <- cor_matrix[, 1]

VarKeep <- c('momage', 'ppvt.imp', 'momwhite', 'momed4F', 'mom.scoll', 'momblack', 'mom.lths', 'b.marry', 'birth.o', 'parity', 'site4')
Col_ParDep <- c('momage', 'ppvt.imp', 'birth.o', 'parity')
Y <-  DATA[,c('iqsb.36')]

# Subset your data to keep only the top 10 predictor columns
selected_data <- DATA[, VarKeep]

## reduced model
output <- DATA[,c('iqsb.36', 'treat')]
DAT_reduced <- bind_cols(output, selected_data)

DAT <- bind_cols(output, selected_data)

tst <- gbm(iqsb.36 ~ ., data=DAT)

fit_boost_func <- function(data_full, data_reduced){
  #Xmat_tmp <- model.matrix(iqsb.36 ~ . - 1, data=data_full)
  #X0mat_notrt_tmp <- model.matrix(iqsb.36 ~ . - treat - 1, data=data_reduced)
  
  tmp_glmboost <- gbm(iqsb.36 ~ ., data=data_full)
  f0fn <- function(tau) {
    data_reduced$Wtau <- data_reduced$iqsb.36 - tau*data_reduced$treat
    fit_reduced <- gbm(Wtau ~ . - treat - iqsb.36, data=data_reduced)
    mse_local <- mean((data_reduced$Wtau - predict(fit_reduced, newdata=data_reduced))^2)
    return(mse_local)
  }
  tau.star.gboost <- optimize(f0fn, interval=c(-20, 20))$minimum
  data_reduced$Wtau <- data_reduced$iqsb.36 - tau.star.gboost*data_reduced$treat
  tmp_reduced_glmboost <- gbm(Wtau ~ . - treat - iqsb.36, data=data_reduced)
  return(list(tmp_glmboost, tmp_reduced_glmboost, tau.star.gboost))
}

fit_func <- fit_boost_func(DAT, DAT_reduced)

## Setup

ngrid <- 200
GridEnd <-  rbind(c(13,43), c(38, 160), c(1, 8), c(1, 8))
ff_full <- matrix(NA, nrow = ngrid, ncol = 4)
ff_reduced <- matrix(NA, nrow = ngrid, ncol = 4)
partial_results_full <- partial_results_reduced <- list()
DAT_tmp <- DAT
DAT_reduced_tmp <- DAT_reduced
for(i in 1:nrow(GridEnd)){
  pp <- seq(GridEnd[i,1], GridEnd[i,2], length.out=ngrid)
  for(k in 1:ngrid){
    
    DAT_tmp[,colnames(DAT_tmp)==Col_ParDep[i]] <- rep(pp[k], nrow(DAT))
    DAT_reduced_tmp[,colnames(DAT_reduced_tmp)==Col_ParDep[i]] <- rep(pp[k], nrow(DAT_reduced))
    
    full_glmboost <- as.numeric(predict(fit_func[[1]], newdata = DAT_tmp))
    full_glmboost_0 <- full_glmboost[DAT_tmp$treat == 0]
    full_glmboost_1 <- full_glmboost[DAT_tmp$treat == 1]
    
    reduced_glmboost <- as.numeric(predict(fit_func[[2]], newdata=DAT_reduced_tmp)) + fit_func[[3]]*DAT_reduced_tmp$treat
    reduced_glmboost_0 <- reduced_glmboost[Xmat_full[,1] == 0]
    reduced_glmboost_1 <- reduced_glmboost[Xmat_full[,1] == 1]
    
    ff_full[k,1] <- mean(full_glmboost_0)
    ff_full[k,2] <- mean(full_glmboost_1)
    ff_full[k,3] <- pp[k]
    #ff_full[k,4] <- Col_ParDep[i]
    ff_reduced[k,1] <- mean(reduced_glmboost_0)
    ff_reduced[k,2] <- mean(reduced_glmboost_1)
    ff_reduced[k,3] <- pp[k]
    #ff_reduced[k,4] <- Col_ParDep[i]
    print(c(i, k))
  }
  partial_results_full[[length(partial_results_full)+1]] <- ff_full
  partial_results_reduced[[length(partial_results_reduced)+1]] <- ff_reduced
}

df_trt1 <- subset(DAT, treat==1)
df_trt0 <- subset(DAT, treat==0)

### Age
fit_trt_1 <- as.data.frame(partial_results_reduced[[1]])
fit_trt_full_1 <- as.data.frame(partial_results_full[[1]])
p1 <- ggplot() +
  geom_point(data = df_trt1, aes(x = momage, y = iqsb.36), color = "blue", alpha = 0.3, size = 2) +
  geom_point(data = df_trt0, aes(x = momage, y = iqsb.36), color = "red", alpha = 0.3, size = 2) +
  geom_line(data = fit_trt_1, aes(x = as.numeric(fit_trt_1[,3]), y = as.numeric(fit_trt_1[,2]), color = "Treatment 1", linetype = "Unrestricted"), size = 1) +
  geom_line(data = fit_trt_1, aes(x = as.numeric(fit_trt_1[,3]), y = as.numeric(fit_trt_1[,1]), color = "Treatment 0", linetype = "Unrestricted"), size = 1) +
  geom_line(data = fit_trt_full_1, aes(x = as.numeric(fit_trt_full_1[,3]), y = as.numeric(fit_trt_full_1[,2]), color = "Treatment 1", linetype = "Restricted"), size = 1) +
  geom_line(data = fit_trt_full_1, aes(x = as.numeric(fit_trt_full_1[,3]), y = as.numeric(fit_trt_full_1[,1]), color = "Treatment 0", linetype = "Restricted"), size = 1) +
  labs(x = "Age", y = "Predicted IQ at 36 months", color = "Treatment", linetype = "Method") +
  scale_color_manual(values = c("Treatment 1" = "blue", "Treatment 0" = "red")) +
  scale_linetype_manual(values = c("Restricted" = "dotted", "Unrestricted" = "solid")) +
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

### ppvt.imp
fit_trt_2 <- as.data.frame(partial_results_reduced[[2]])
fit_trt_full_2 <- as.data.frame(partial_results_full[[2]])
p2 <- ggplot() +
  geom_point(data = df_trt1, aes(x = ppvt.imp, y = iqsb.36), color = "blue", alpha = 0.3, size = 2) +
  geom_point(data = df_trt0, aes(x = ppvt.imp, y = iqsb.36), color = "red", alpha = 0.3, size = 2) +
  geom_line(data = fit_trt_2, aes(x = as.numeric(fit_trt_2[,3]), y = as.numeric(fit_trt_2[,2]), color = "Treatment 1", linetype = "Unrestricted"), size = 1) +
  geom_line(data = fit_trt_2, aes(x = as.numeric(fit_trt_2[,3]), y = as.numeric(fit_trt_2[,1]), color = "Treatment 0", linetype = "Unrestricted"), size = 1) +
  geom_line(data = fit_trt_full_2, aes(x = as.numeric(fit_trt_full_2[,3]), y = as.numeric(fit_trt_full_2[,2]), color = "Treatment 1", linetype = "Restricted"), size = 1) +
  geom_line(data = fit_trt_full_2, aes(x = as.numeric(fit_trt_full_2[,3]), y = as.numeric(fit_trt_full_2[,1]), color = "Treatment 0", linetype = "Restricted"), size = 1) +
  labs(x = "PPVT", y = "", color = "Treatment", linetype = "Method") +
  scale_color_manual(values = c("Treatment 1" = "blue", "Treatment 0" = "red")) +
  scale_linetype_manual(values = c("Restricted" = "dotted", "Unrestricted" = "solid")) +
  theme_minimal()+
  theme(
    text = element_text(size = 20),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 20),
    plot.title = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 24),
    legend.text = element_text(size = 20),
    axis.text.y = element_blank(),
    legend.position = "bottom"
  )

## birth.o
fit_trt_3 <- as.data.frame(partial_results_reduced[[3]])
fit_trt_full_3 <- as.data.frame(partial_results_full[[3]])
p3 <- ggplot() +
  geom_point(data = df_trt1, aes(x = birth.o, y = iqsb.36), color = "blue", alpha = 0.3, size = 2) +
  geom_point(data = df_trt0, aes(x = birth.o, y = iqsb.36), color = "red", alpha = 0.3, size = 2) +
  geom_line(data = fit_trt_3, aes(x = as.numeric(fit_trt_3[,3]), y = as.numeric(fit_trt_3[,2]), color = "Treatment 1", linetype = "Unrestricted"),  size = 1) +
  geom_line(data = fit_trt_3, aes(x = as.numeric(fit_trt_3[,3]), y = as.numeric(fit_trt_3[,1]), color = "Treatment 0", linetype = "Unrestricted"), size = 1) +
  geom_line(data = fit_trt_full_3, aes(x = as.numeric(fit_trt_full_3[,3]), y = as.numeric(fit_trt_full_3[,2]), color = "Treatment 1", linetype = "Restricted"),  size = 1) +
  geom_line(data = fit_trt_full_3, aes(x = as.numeric(fit_trt_full_3[,3]), y = as.numeric(fit_trt_full_3[,1]), color = "Treatment 0", linetype = "Restricted"),  size = 1) +
  labs(x = "Birth order", y = "Predicted IQ at 36 months", color = "Treatment", linetype = "Method") +
  scale_color_manual(values = c("Treatment 1" = "blue", "Treatment 0" = "red")) +
  scale_linetype_manual(values = c("Restricted" = "dotted", "Unrestricted" = "solid")) +
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

## parity

fit_trt_4 <- as.data.frame(partial_results_reduced[[4]])
fit_trt_full_4 <- as.data.frame(partial_results_full[[4]])
p4 <- ggplot() +
  geom_point(data = df_trt1, aes(x = parity, y = iqsb.36), color = "blue", alpha = 0.3, size = 2) +
  geom_point(data = df_trt0, aes(x = parity, y = iqsb.36), color = "red", alpha = 0.3, size = 2) +
  geom_line(data = fit_trt_4, aes(x = as.numeric(fit_trt_4[,3]), y = as.numeric(fit_trt_4[,2]), color = "Treatment 1", linetype = "Unrestricted"),  size = 1) +
  geom_line(data = fit_trt_4, aes(x = as.numeric(fit_trt_4[,3]), y = as.numeric(fit_trt_4[,1]), color = "Treatment 0", linetype = "Unrestricted"),  size = 1) +
  geom_line(data = fit_trt_full_4, aes(x = as.numeric(fit_trt_full_4[,3]), y = as.numeric(fit_trt_full_4[,2]), color = "Treatment 1", linetype = "Restricted"),  size = 1) +
  geom_line(data = fit_trt_full_4, aes(x = as.numeric(fit_trt_full_4[,3]), y = as.numeric(fit_trt_full_4[,1]), color = "Treatment 0", linetype = "Restricted"), size = 1) +
  labs(x = "Number of childeren", y = "", color = "Treatment", linetype = "Method") +
  scale_color_manual(values = c("Treatment 1" = "blue", "Treatment 0" = "red")) +
  scale_linetype_manual(values = c("Restricted" = "dotted", "Unrestricted" = "solid")) +
  theme_minimal()+
  theme(
    text = element_text(size = 22),
    axis.title = element_text(size = 24),
    axis.text = element_text(size = 22),
    plot.title = element_text(size = 18, face = "bold"),
    legend.title = element_text(size = 26),
    legend.text = element_text(size = 22),
    axis.text.y = element_blank(),
    legend.position = "bottom"
  )

# Install and load the patchwork package if you haven't already
# install.packages("patchwork")
library(patchwork)


# Remove the legends from individual plots
p1 <- p1 + theme(legend.position = "none")
p2 <- p2 + theme(legend.position = "none")
p3 <- p3 + theme(legend.position = "none")
p4 <- p4 + theme(legend.position = "none")

# Extract legend from one of the plots
get_legend <- function(plot) {
  g <- ggplotGrob(plot)
  legend <- g$grobs[which(g$layout$name == "guide-box")]
  legend
}

legend <- get_legend(p1)

# Combine plots without individual legends
combined_plot <- (p1 + p2) / (p3 + p4) +
  plot_layout(guides = 'collect') &
  theme(
    legend.position = 'bottom',
    legend.title = element_text(size = 26),
    legend.text = element_text(size = 22)
  )

# Print the combined plot with the legend
print(combined_plot)




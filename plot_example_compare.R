
################
### On the simulated dataset
################

set.seed(123)

# Define the number of observations
n <- 50

# Generate random values for x1 and x2 from a normal distribution
x1 <- rnorm(n, mean = 0, sd = 1)
x2 <- rnorm(n, mean = 0, sd = 1)

# Generate treatment variable A
A <- rbinom(n, size = 1, prob = 0.5)

# Define the coefficients
beta0 <- 2
beta1 <- 3
beta2 <- -1
beta3 <- 1.5
beta4 <- 0.5
beta5 <- -2

# Generate random error terms
epsilon <- rnorm(n, mean = 0, sd = 1)

Y <- beta0 + beta1 * x1 + beta2 * x2 + beta3 * A + beta4 * A * x1 + beta5 * A * x2 + epsilon


DATA_full <- data.frame('Y' = Y, 'x1' = x1, 'x2' = x2, 'A' = A, 'x1.t' = A*x1, 'x2.t' = A*x2)
DATA_reduced <- data.frame('Y' = Y, 'x1' = x1, 'x2' = x2, 'A' = A)


Y <- DATA_full$Y
DATA_full <- DATA_full[ , !(names(DATA_full) %in% c('Y'))]
DATA_full <- model.matrix(Y~.-1, data = DATA_full)


Treat <- DATA_reduced$A
DATA_reduced <- DATA_reduced[ , !(names(DATA_reduced) %in% c('Y', 'A'))]
DATA_reduced <- model.matrix(Y~.-1, data = DATA_reduced)
tau <- mean(Treat == 1) - mean(Treat == 0)
tau.range <- seq(-4*tau, 4*tau, length.out = 9)
#tau.range = seq(1,10, by =1)


data.all <- cbind.data.frame(DATA_full, Y = Y)
fit <- lm(Y ~., data = data.all)
pred <- predict(fit)

plot_data <- data.frame(Observed = data.all$Y, Predicted = pred, Trt = Treat)
# Subset dataframe based on Trt values
df_trt1 <- plot_data[plot_data$Trt == 1, ]
df_trt0 <- plot_data[plot_data$Trt == 0, ]


mse <- rep(NA, length(tau.range))
for(k in 1:length(tau.range)) {
  data.reduced <- cbind.data.frame(DATA_reduced, Y= Y - tau.range[k]*Treat)
  lm_tmp <- lm(Y ~., data = data.reduced)
  mse[k] <- mean((data.reduced$Y - lm_tmp$fitted)^2)
}
tau.star <- tau.range[which.min(mse)]

data.reduced <- cbind.data.frame(DATA_reduced, Y= Y - tau.star*Treat)
fit_reduced <- lm(Y ~., data = data.reduced)
pred_reduced <- predict(fit_reduced)

plot_data_reduced <- data.frame(Observed = data.reduced$Y, Predicted = pred_reduced, Trt = Treat)
# Subset dataframe based on Trt values
df_reduced_trt1 <- plot_data_reduced[plot_data_reduced$Trt == 1, ]
df_reduced_trt0 <- plot_data_reduced[plot_data_reduced$Trt == 0, ]



library(ggplot2)
ggplot() +
  geom_line(data = df_trt1, aes(x = Predicted, y = Observed), color = "blue") +
  geom_point(data = df_reduced_trt1, aes(y = Predicted, x = Observed), color = "blue", alpha = 0.5) +
  
  # Add df_trt0 as a line with points from df_reduced_trt0 around it
  geom_line(data = df_trt0, aes(x = Predicted, y = Observed), color = "red") +
  geom_point(data = df_reduced_trt0, aes(y = Predicted, x = Observed), color = "red", alpha = 0.5) +
  
  # Add labels to the lines
  geom_text(data = df_trt1[1, ], aes(x = Predicted, y = Observed, label = "Treat = 1"), hjust = -1.2, vjust = 0.2, color = "blue") +
  geom_text(data = df_trt0[1, ], aes(x = Predicted, y = Observed, label = "Treat = 0"), hjust = 3.0, vjust = -5.5, color = "red") +
  
  # Customize plot aesthetics
  labs(y = "Predicted", x = "Observed", color = "Treatment") +
  scale_color_manual(values = c("blue", "red"), labels = c("Trt = 1", "Trt = 0")) +
  theme_minimal()





################
### On the example dataset
################
DATA <- read.csv('IHDP_clean.csv', header = TRUE)[,-1]
cor_matrix <- cor(DATA)
cor_with_output <- cor_matrix[, 1]
cor_with_output[1] <- 0  # Set the correlation with output to 0
sorted_correlations <- sort(abs(cor_with_output), decreasing = TRUE)
top_10_correlations <- sorted_correlations[1:10]
top_10_predictors <- names(top_10_correlations)
# Subset your data to keep only the top 10 predictor columns

selected_data <- DATA[, top_10_predictors]
output <- DATA[,c('iqsb.36', 'treat')]

## X0 - reduced model
DATA.cor.reduced <- bind_cols(output, selected_data)
# Create training and test sets using the selected index numbers
Y <- DATA.cor.reduced$iqsb.36
Treat <- DATA.cor.reduced$treat
tau <- mean(Treat == 1) - mean(Treat == 0)
tau.range <- seq(-4*tau, 4*tau, length.out = 9)
DATA.cor.reduced <- DATA.cor.reduced[ , !(names(DATA.cor.reduced) %in% c('iqsb.36', 'treat'))]
DATA.cor.reduced <- model.matrix(Y~.-1, data = DATA.cor.reduced)

## X - full model
selected_data.T <- matrix(NA, nrow = nrow(selected_data), ncol = ncol(selected_data))
for(k in 1:10){
  selected_data.T[,k] <- selected_data[,k]*output$treat
}
colnames(selected_data.T) <- paste0(colnames(selected_data), ".t")
DATA.cor <- bind_cols(output, selected_data, selected_data.T)
DATA.cor <- bind_cols('iqsb.36' = output, selected_data)

# Create training and test sets using the selected index numbers
Y <- DATA.cor$iqsb.36
DATA.cor <- DATA.cor[ , !(names(DATA.cor) %in% c('iqsb.36'))]
DATA.cor <- model.matrix(Y~.-1, data = DATA.cor)


data.all <- cbind.data.frame(DATA.cor, Y = Y)
fit <- lm(Y ~., data = data.all)
pred <- predict(fit)

plot_data <- data.frame(Observed = data.all$Y, Predicted = as.vector(pred), Trt = data.all$treat)
# Subset dataframe based on Trt values
df_trt1 <- plot_data[plot_data$Trt == 1, ]
df_trt0 <- plot_data[plot_data$Trt == 0, ]


mse <- rep(NA, length(tau.range))
for(k in 1:length(tau.range)) {
  data.reduced <- cbind.data.frame(DATA.cor.reduced, Y= Y - tau.range[k]*Treat)
  lm_tmp <- lm(Y ~., data = data.reduced)
  mse[k] <- mean((data.reduced$Y - lm_tmp$fitted)^2)
}
tau.star <- tau.range[which.min(mse)]

data.reduced <- cbind.data.frame(DATA.cor.reduced, Y= Y - tau.star*Treat)
fit_reduced <- lm(Y ~., data = data.reduced)
pred_reduced <- predict(fit_reduced)

plot_data_reduced <- data.frame(Observed = data.reduced$Y, Predicted = pred_reduced, Trt = Treat)
# Subset dataframe based on Trt values
df_reduced_trt1 <- plot_data_reduced[plot_data_reduced$Trt == 1, ]
df_reduced_trt0 <- plot_data_reduced[plot_data_reduced$Trt == 0, ]



library(ggplot2)
ggplot() +
  geom_line(data = df_trt1, aes(x = Predicted, y = Observed), color = "blue") +
  geom_point(data = df_reduced_trt1, aes(y = Predicted, x = Observed), color = "blue", alpha = 0.5) +
  
  # Add df_trt0 as a line with points from df_reduced_trt0 around it
  geom_line(data = df_trt0, aes(x = Predicted, y = Observed), color = "red") +
  geom_point(data = df_reduced_trt0, aes(y = Predicted, x = Observed), color = "red", alpha = 0.5) +
  
  # Add labels to the lines
  geom_text(data = df_trt1[1, ], aes(x = Predicted, y = Observed, label = "Treat = 1"), hjust = -1.2, vjust = 0.2, color = "blue") +
  geom_text(data = df_trt0[1, ], aes(x = Predicted, y = Observed, label = "Treat = 0"), hjust = 3.0, vjust = -5.5, color = "red") +
  
  # Customize plot aesthetics
  labs(y = "Predicted", x = "Observed", color = "Treatment") +
  scale_color_manual(values = c("blue", "red"), labels = c("Trt = 1", "Trt = 0")) +
  theme_minimal()








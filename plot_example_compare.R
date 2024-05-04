
################
### On the simulated dataset
################

set.seed(123)

# Define the number of observations
n <- 500

# Generate random values for x1 and x2 from a normal distribution
x <- runif(n)


# Generate treatment variable A
A <- rbinom(n, size = 1, prob = 0.5)

# Define the coefficients
beta0 <- 2
beta1 <- 5
beta2 <- -1
beta3 <- 5.5


# Generate random error terms
epsilon <- rnorm(n, mean = 0, sd = 0.25)

Y <- beta0 + beta1 * (x-0.5)^3 + beta2 * A + beta3 * A * (x-0.5)^3 + epsilon


DATA_full <- data.frame('Y' = Y, 'x' = x, 'A' = A, 'x.t' = A *x)
DATA_reduced <- data.frame('Y' = Y, 'x' = x, 'A' = A)


Y <- DATA_full$Y
#DATA_full <- DATA_full[ , !(names(DATA_full) %in% c('Y'))]
#DATA_full <- model.matrix(Y~.-1, data = DATA_full)


Treat <- DATA_reduced$A
DATA_reduced <- DATA_reduced[ , !(names(DATA_reduced) %in% c('Y', 'A'))]
DATA_reduced <- model.matrix(Y~.-1, data = as.data.frame(DATA_reduced))
tau <- mean(Y[Treat == 1]) - mean(Y[Treat == 0])
tau.range <- seq(-4*tau, 4*tau, length.out = 20)



#data.all <- cbind.data.frame(DATA_full, Y = Y)
fit <- lm(Y ~., data = DATA_full)
pred <- predict(fit)

plot_data <- data.frame(coef = DATA_full$x, Observed = Y, Predicted = pred, Trt = Treat)
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
pred_reduced <- predict(fit_reduced) + tau.star*Treat

plot_data_reduced <- data.frame(coef = data.reduced$DATA_reduced, Observed = Y, Predicted = pred_reduced, Trt = Treat)
# Subset dataframe based on Trt values
df_reduced_trt1 <- plot_data_reduced[plot_data_reduced$Trt == 1, ]
df_reduced_trt0 <- plot_data_reduced[plot_data_reduced$Trt == 0, ]



library(ggplot2)
ggplot() +
  geom_line(data = df_trt1, aes(x = coef, y = pred), color = "blue") +
  #geom_point(data = df_reduced_trt1, aes(x = Predicted, y = Observed), color = "blue", alpha = 0.5) +
  
  # Add df_trt0 as a line with points from df_reduced_trt0 around it
  geom_line(data = df_trt0, aes(x = coef, y = pred), color = "red") +
  #geom_point(data = df_reduced_trt0, aes(x = Predicted, y = Observed), color = "red", alpha = 0.5) +
  
  # Add labels to the lines
  geom_text(data = df_trt1[1, ], aes(x = Predicted, y = Observed, label = "Treat = 1"), hjust = -1.2, vjust = 0.2, color = "blue") +
  geom_text(data = df_trt0[1, ], aes(x = Predicted, y = Observed, label = "Treat = 0"), hjust = 3.0, vjust = -5.5, color = "red") +
  
  # Customize plot aesthetics
  labs(y = "Predicted", x = "Observed", color = "Treatment") +
  scale_color_manual(values = c("blue", "red"), labels = c("Trt = 1", "Trt = 0")) +
  theme_minimal()






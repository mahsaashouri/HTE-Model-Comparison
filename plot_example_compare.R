set.seed(123)
library(splines)
library(ggplot2)

# Define the number of observations
n <- 500

# Generate random values for x1 and x2 from a normal distribution
x <- runif(n)


# Generate treatment variable A
A <- rbinom(n, size = 1, prob = 0.5)

# Define the coefficients
beta0 <- 1
beta1 <- 5
beta2 <- 1
beta3 <- 5.5


# Generate random error terms
epsilon <- rnorm(n, mean = 0, sd = 0.3)

Y <- beta0 + beta1 * (x-0.5)^3 + beta2 * A + beta3 * A * (x-0.5)^3 + A*(x > 0.5)*(x - 0.5)^3 + epsilon


DATA_full <- data.frame('Y' = Y, 'x' = x, 'A' = A, 'x.t' = A *x)
DATA_reduced <- data.frame('Y' = Y, 'x' = x, 'A' = A)

bs_fit <- lm(Y ~ bs(x, df=5) + A*bs(x, df=5), data=DATA_full)
test_points <- data.frame(x=rep(seq(0.01, 0.99, length.out=500), 2), A=rep(c(0, 1), each=500))
test_points$fitted_vals <- predict(bs_fit, newdata=test_points)


tau <- mean(Y[A == 1]) - mean(Y[A == 0])
tau.range <- seq(-4*tau, 4*tau, length.out = 20)

mse <- rep(NA, length(tau.range))
for(k in 1:length(tau.range)) {
  tau.tmp <- tau.range[k]
  data.reduced <- data.frame('Y' = Y - tau.tmp*A, 'x' = x)
  lm_tmp <- lm(Y ~ bs(x), data = data.reduced)
  ghat <- lm_tmp$fitted.values + tau.tmp*A
  mse[k] <- mean((Y - ghat)^2)
}
tau.star <- tau.range[which.min(mse)]

data.reduced <- data.frame('Y' = Y - tau.star*A, 'x' = x)
lm_star <- lm(Y ~ bs(x), data = data.reduced)

test_points_reduced0 <- test_points_reduced1 <- data.frame(x=rep(seq(0.01, 0.99, length.out=500), 2))
test_points_reduced0$fitted_vals <- predict(lm_star, newdata=test_points)
test_points_reduced1$fitted_vals <- test_points_reduced0$fitted_vals + tau.star

df_trt1 <- subset(DATA_full, A==1)
df_trt0 <- subset(DATA_full, A==0)
fit_trt1 <- subset(test_points, A==1)
fit_trt0 <- subset(test_points, A==0)


ggplot() + 
  geom_point(data = df_trt1, aes(x = x, y = Y), color = "blue", alpha = 0.5) +
  geom_point(data = df_trt0, aes(x = x, y = Y), color = "red", alpha = 0.5) +
  geom_line(data = fit_trt1, aes(x = x, y = fitted_vals), color = "blue", alpha = 0.5) +
  geom_line(data = fit_trt0, aes(x = x, y = fitted_vals), color = "red", alpha = 0.5) +
  geom_line(data = test_points_reduced1, aes(x = x, y = fitted_vals), color = "blue", alpha = 0.5, linetype = 2) +
  geom_line(data = test_points_reduced0, aes(x = x, y = fitted_vals), color = "red", alpha = 0.5, linetype = 2) +
  geom_text(data = data.frame(x = 0.0, Y = 4.5), aes(x = x, y = Y), label = "Treatment = 1", hjust = -0.5, vjust = 0.2, color = "blue", size = 7) +
  geom_text(data = data.frame(x = 0.0, Y = 4.2), aes(x = x, y = Y), label = "Treatment = 0", hjust = -0.5, vjust = 0.2, color = "red", size = 7) +
  annotate("segment", x = 0.25, xend = 0.3, y = 3.8, yend = 3.8, color = "gray40", size = 0.8) +
  annotate("text", x = 0.015, y = 3.8, label = "Restricted - ", hjust = -0.5, vjust = 0.2, color = "gray40", size = 7) +
  annotate("text", x = 0.15, y = 3.8, label = expression(hat(f)(bold(x), A)), hjust = -0.5, vjust = 0.2, color = "gray40", size = 7) +
  annotate("segment", x = 0.27, xend = 0.32, y = 3.5, yend = 3.5, color = "gray40", size = 0.8, linetype = "dashed") +
  annotate("text", x = 0.0051, y = 3.5, label = "Unrestricted - ", hjust = -0.5, vjust = 0.2, color = "gray40", size = 7) +
  annotate("text", x = 0.17, y = 3.5, label = expression(hat(g)(bold(x), A)), hjust = -0.5, vjust = 0.2, color = "gray40", size = 7) +
  labs(y = "Outcome", x = "x", color = "Treatment") +
  scale_color_manual(values = c("blue", "red"), labels = c("Treatment = 1", "Trt = 0")) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 25)) +
  geom_segment(aes(x = 0.88, xend = 0.88, y = 1.39, yend = 2.50), 
               arrow = arrow(type = "open", length = unit(0.2, "inches")), 
               color = "black", lineend = "round", size = 0.8) +
  geom_segment(aes(x = 0.88, xend = 0.88, y = 2.50, yend = 1.39), 
               arrow = arrow(type = "open", length = unit(0.2, "inches")), 
               color = "black", lineend = "round", size = 0.8) +
  annotate("text", x = 0.91, y = 1.6, label = expression(tau^'*'), color = "black", size = 13, hjust = 0.5, vjust = -0.5)


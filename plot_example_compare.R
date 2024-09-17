library(splines)
library(ggplot2)
library(gridExtra)

set.seed(123)
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
epsilon <- rnorm(n, mean = 0, sd = 0.5)

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

test_points_reduced0 <- test_points_reduced1 <- data.frame(x=seq(0.01, 0.99, length.out=500))
test_points_reduced0$fitted_vals <- predict(lm_star, newdata=test_points[test_points$A==0,])
test_points_reduced1$fitted_vals <- test_points_reduced0$fitted_vals + tau.star
g0 <- approxfun(test_points_reduced0$x, test_points_reduced0$fitted_vals)
g1 <- approxfun(test_points_reduced1$x, test_points_reduced1$fitted_vals)

df_trt1 <- subset(DATA_full, A==1)
df_trt0 <- subset(DATA_full, A==0)

fit_trt1 <- subset(test_points, A==1)
fit_trt0 <- subset(test_points, A==0)


txt_size <- 4
pt_size <- 1.5
xtau <- 0.94
ytau0 <- g0(xtau)
ytau1 <- g1(xtau)


p1 <- ggplot() +
  geom_point(data = df_trt1, aes(x = x, y = Y), color = "blue", alpha = 0.5, size = pt_size) +
  geom_point(data = df_trt0, aes(x = x, y = Y), color = "red", alpha = 0.5, size = pt_size) +
  geom_line(data = fit_trt1, aes(x = x, y = fitted_vals), color = "blue", alpha = 1, linewidth = 1.5) +
  geom_line(data = fit_trt0, aes(x = x, y = fitted_vals), color = "red", alpha = 1, linewidth =1.5) +
  geom_line(data = test_points_reduced1, aes(x = x, y = fitted_vals), color = "blue", alpha = 1, linetype = 2, linewidth = 1.5) +
  geom_line(data = test_points_reduced0, aes(x = x, y = fitted_vals), color = "red", alpha = 1, linetype = 2, linewidth = 1.5) +
  geom_text(data = data.frame(x = 0.0, Y = 4.5), aes(x = x, y = Y), label = "Treatment = 1", hjust = -0.5, vjust = 0.2, color = "blue", size = txt_size) +
  geom_text(data = data.frame(x = 0.0, Y = 4.2), aes(x = x, y = Y), label = "Treatment = 0", hjust = -0.5, vjust = 0.2, color = "red", size = txt_size) +
  annotate("segment", x = 0.01, xend = 0.09, y = 3.8, yend = 3.8, color = "gray40", size = 0.8, linetype = "dashed") +
  annotate("text", x = 0.29, y = 3.8, label = expression(paste("Restricted: ", hat(g)(bold(x), A))), color = "gray40", size = txt_size) +
  #annotate("text", x = 0.15, y = 3.8, label = expression(hat(g)(bold(x), A)), hjust = -0.5, vjust = 0.2, color = "gray40", size = txt_size) +
  annotate("segment", x = 0.01, xend = 0.09, y = 3.5, yend = 3.5, color = "gray40", size = 1) +
  annotate("text", x = 0.305, y = 3.52, label = expression(paste("Unrestricted: ", hat(f)(bold(x), A))), color = "gray40", size = txt_size) +
  #annotate("text", x = 0.17, y = 3.5, label = expression(hat(f)(bold(x), A)), hjust = -0.5, vjust = 0.2, color = "gray40", size = txt_size) +
  labs(y = "Outcome", x = "x", color = "Treatment") +
  scale_color_manual(values = c("blue", "red"), labels = c("Treatment = 1", "Trt = 0")) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 14)) +
  geom_segment(aes(x = xtau, xend = xtau, y = ytau0, yend = ytau1),
               arrow = arrow(type = "open", length = unit(0.2, "inches")),
               color = "black", lineend = "round", size = 0.8) +
  geom_segment(aes(x = xtau, xend = xtau, y = ytau1, yend = ytau0),
               arrow = arrow(type = "open", length = unit(0.2, "inches")),
               color = "black", lineend = "round", size = 0.8) +
  annotate("text", x = xtau + 0.07, y = (ytau0 + ytau1)/2 + 0.1, label = expression(tau^'*'), color = "black", size = 12) +
  ylim(-0.7,4.5)

p2 <- ggplot() +
  geom_point(data = df_trt1, aes(x = x, y = Y), color = "blue", alpha = 0.5, size = pt_size) +
  geom_line(data = fit_trt0, aes(x = x, y = fitted_vals), color = "darkorange2", alpha = 1, linewidth =1.5) +
  geom_point(data = data.reduced, aes(x = x, y = Y), color = "chartreuse3", alpha = 0.5, size = pt_size) +
  annotate("segment", x = 0.01, xend = 0.09, y = 3.6, yend = 3.6, color = "darkorange2", size = 1) +
  annotate("text", x = 0.50, y = 3.62, label = expression(paste("Restricted baseline risk estimate: ", hat(f)[tau](bold(x)))),
           color = "gray40", size = txt_size) +
  annotate(geom="point", x = 0.055, y = 3.92, colour = "chartreuse3", size = txt_size) +
  annotate("text", x = 0.42, y = 3.92, label=expression(paste(M[tau], ": Modified-outcome dataset")),
           color = "gray40", size = txt_size) +
  annotate(geom="point", x = 0.055, y = 4.28, colour = "blue", size = txt_size) +
  annotate("text", x = 0.44, y = 4.28, label = "Outcomes under Treatment = 1",
           color = "gray40", size = txt_size) +
  labs(y = "Outcome", x = "x") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 14)) +
  ylim(-0.7,4.5)

pcombined <- grid.arrange(p1, p2, nrow = 1)
pcombined
## Save plot in a 1.5 x 1 aspect ratio


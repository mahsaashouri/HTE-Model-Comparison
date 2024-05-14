
library(ggplot2)
library(gridExtra)

## Sample size 100
mu <- c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4)
tau <- c(1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5)
## linear
linear <- c(0.6, 0.6, 0.6, 0.8, 0.6, 0.8, 0.8, 0.8, 0.8, 1, 0.8, 0.8, 0.8, 0.8, 1, 0.8, 0.8, 0.8, 0.8, 1)
## glmnet
glmnet <- c(1, 1, 1, 1, 1, 0.8, 0.6, 0.2, 0.8, 0.8, 0.8, 1, 0.6, 0.8, 1, 1, 0.8, 0.4, 0.8, 0.6)
## glmboost
glmboost <- c(0.8, 0.8, 0.8, 0.8, 0.8, 1, 1, 1, 0.8, 1, 0.8, 0.8, 0.8, 0.8, 1, 1, 1, 1, 0.8, 1)

data_100 <- data.frame(
  Method = rep(c("linear", "glmnet", "glmboost"), each = length(linear)),
  Value = c(1-linear, 1-glmnet, 1-glmboost),
  SampleSize = rep("100", each = length(linear))
)

## Sample size 500
mu <- c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4)
tau <- c(1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5)
## linear
linear <- c(0.6, 0.6, 0.6, 0.8, 0.6, 0.8, 0.8, 0.8, 1, 0.8, 0.4, 0.4, 0.4, 1, 0.4, 0.6, 0.6, 0.6, 1, 0.8)
## glmnet
glmnet <- c(0.8, 0.8, 1, 1, 0.8, 0.8, 0.8, 0, 0.4, 0.8, 0.4, 0.8, 1, 0.4, 0.4, 0.6, 0.6, 0, 0.2, 1)
## glmboost
glmboost <- c(0.8, 0.8, 0.8, 1, 1, 0.8, 0.8, 0.8, 1, 0.8, 0.4, 0.2, 0.2, 1, 0.4, 0.8, 0.6, 0.6, 1, 0.8)

data_500 <- data.frame(
  Method = rep(c("linear", "glmnet", "glmboost"), each = length(linear)),
  Value = c(1-linear, 1-glmnet, 1-glmboost),
  SampleSize = rep("500", each = length(linear))
)

method_order <- c("linear", "glmnet", "glmboost")

# Combine the data for both sample sizes
combined_data <- rbind(data_100, data_500)

# Plot with facets
ggplot(combined_data, aes(x = Method, y = Value, fill = Method)) +
  geom_violin(fill = 'gray') +
  #geom_boxplot() +
  labs(x = "Method",
       y = "Value") +
  facet_wrap(~ SampleSize, ncol = 2, labeller = as_labeller(c("100" = "Sample Size = 100", "500" = "Sample Size = 500"))) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 22), 
        strip.text = element_text(size = 18),
        text = element_text(size = 20)) +
  scale_x_discrete(limits = method_order)


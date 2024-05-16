
library(ggplot2)
library(gridExtra)

## Sample size 100
mu <- c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4)
theta <- c(1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5)
## linear
linear <- c(0.6, 0.6, 0.6, 0.8, 0.6, 0.8, 0.8, 0.8, 0.8, 1, 0.8, 0.8, 0.8, 0.8, 1, 0.8, 0.8, 0.8, 0.8, 1)
## glmnet
glmnet <- c(1, 1, 1, 1, 1, 0.8, 0.6, 0.2, 0.8, 0.8, 0.8, 1, 0.6, 0.8, 1, 1, 0.8, 0.4, 0.8, 0.6)
## glmboost
glmboost <- c(0.8, 0.8, 0.8, 0.8, 0.8, 1, 1, 1, 0.8, 1, 0.8, 0.8, 0.8, 0.8, 1, 1, 1, 1, 0.8, 1)

data_100 <- data.frame(
  mu = mu,
  theta = theta, 
  Method = rep(c("linear", "glmnet", "glmboost"), each = length(linear)),
  Value = c(1-linear, 1-glmnet, 1-glmboost),
  SampleSize = rep("100", each = length(linear))
)

## Sample size 500
mu <- c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4)
theta <- c(1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5)
## linear
linear <- c(0.6, 0.6, 0.6, 0.8, 0.6, 0.8, 0.8, 0.8, 1, 0.8, 0.4, 0.4, 0.4, 1, 0.4, 0.6, 0.6, 0.6, 1, 0.8)
## glmnet
glmnet <- c(0.8, 0.8, 1, 1, 0.8, 0.8, 0.8, 0, 0.4, 0.8, 0.4, 0.8, 1, 0.4, 0.4, 0.6, 0.6, 0, 0.2, 1)
## glmboost
glmboost <- c(0.8, 0.8, 0.8, 1, 1, 0.8, 0.8, 0.8, 1, 0.8, 0.4, 0.2, 0.2, 1, 0.4, 0.8, 0.6, 0.6, 1, 0.8)

data_500 <- data.frame(
  mu = mu,
  theta = theta, 
  Method = rep(c("linear", "glmnet", "glmboost"), each = length(linear)),
  Value = c(1-linear, 1-glmnet, 1-glmboost),
  SampleSize = rep("500", each = length(linear))
)

method_order <- c("linear", "glmnet", "glmboost")

combined_data <- rbind(data_100, data_500)

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

#####################
## Separate theta equal to 0 and 1
####################
condition <- data_100$theta <= 2
data_100_1 <- cbind(data_100[condition, ], sampleCon = '100-1')
data_100_2 <- cbind(data_100[!condition, ], sampleCon = '100-2')

condition <- data_500$theta <= 2
data_500_1 <- cbind(data_500[condition, ], sampleCon = '500-1')
data_500_2 <- cbind(data_500[!condition, ], sampleCon = '500-2')

combined_data <- rbind(data_100_1, data_100_2, data_500_1, data_500_2)
method_order <- c("linear", "glmnet", "glmboost")


ggplot(combined_data, aes(x = Method, y = Value, fill = Method)) +
  geom_violin(fill = 'gray') +
  #geom_boxplot() +
  labs(x = "Method",
       y = "Value") +
  facet_wrap(~ sampleCon, ncol = 2, labeller = as_labeller(c("100-1" = "Sample Size = 100 - θ1, θ2", 
                                                              "100-2" = "Sample Size = 100 - θ3, θ4, θ5",
                                                              "500-1" = "Sample Size = 500 - θ1, θ2",
                                                              "500-2" = "Sample Size = 500 - θ3, θ4, θ5"))) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 22), 
        strip.text = element_text(size = 18),
        text = element_text(size = 20)) +
  scale_x_discrete(limits = method_order)


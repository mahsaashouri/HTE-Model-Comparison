
# Define mu and tau 

f1 <- function(x) {
  return((x[, 2]*x[, 4]*x[, 6]) + (2*x[, 2]*x[, 4]*(1-x[, 6])) + (3*x[, 2]*(1-x[, 4])*x[, 6]) + (4*x[, 2]*(1-x[, 4])*(1-x[, 4])) + (5*(1-x[, 2])*x[, 4]*x[, 6]) + 
    (6*(1-x[, 2])*x[, 4]*(1-x[, 6])) + (7*(1-x[, 2])*(1-x[, 4])*x[, 6]) + (8*(1-x[, 2])*(1-x[, 4])*(1-x[, 6])))
}

f2 <- function(x) {
  return((4*ifelse(x[, 1] > 1 & x[, 3] > 0, 1, 0)) + (4*ifelse(x[, 5] > 1 & x[, 7] > 1, 1, 0)) + (x[, 8]*x[, 9]))
}

f3 <- function(x) {
  return(0.5*(x[, 1]^2 + x[, 2] + x[, 3]^2 + x[, 4] + x[, 5]^2 + x[, 6] + x[, 7]^2 + x[, 8] + x[, 9]^2 -11))
}

f4 <- function(x) {
  (1/sqrt(2))*(f1(x) + (x[, 1] + x[, 3] + x[, 5] + x[, 7] + x[, 8] + x[, 9] - 2))
  return(x[, 7])
}

mu <- function(choice, x){
  if (choice == 1) 
    return(f1(x))
  else if (choice == 2) 
    return(f2(x))
  else if (choice == 3) 
    return(f3(x))
  else if (choice == 4) 
    return(f4(x))
  else
    stop("Invalid choice for mu")
}

tau <- function(choice, x) {
  if (choice == 1) {
    return(0)
  } else if (choice == 2) {
    return(1)
  } else if (choice == 3) {
    return(2 + 0.1/(1 + exp(-x[, 2])))
  } else if (choice == 4) {
    return(f1(x))
  } else if (choice == 5) {
    return(f2(x))
  } else {
    stop("Invalid choice for tau")
  }
}

# Generate the dataset
set.seed(123) 

n <- 100  # Number of observations
p <- 10   # Number of features

# Generate x
x <- matrix(0, nrow = n, ncol = p)
for (j in 1:p) {
  if (j %% 2 == 0) {
    x[, j] <- rnorm(n, 0, 1)  # Normal(0, 1) if j is even
  } else {
    x[, j] <- rbinom(n, 1, 0.5)  # Bernoulli(0.5) if j is odd
  }
}

# Generate treatment indicator A 
A <- rbinom(n, 1, 0.5)

# Generate outcome variable Y
mu_new <- mu(4, x) 
tau_new <- tau(3, x) 

Y <- numeric(n)
for (i in 1:n) {
  Y[i] <- mu_new[i] + A[i] * tau_new[i] + rnorm(1)  
}

data <- data.frame(Y = Y, A = A, x)



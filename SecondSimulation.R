
# Define mu and tau 

f4 <- function(x) {
  return((x[, 2]*x[, 4]*x[, 6]) + (2*x[, 2]*x[, 4]*(1-x[, 6])) + (3*x[, 2]*(1-x[, 4])*x[, 6]) + (4*x[, 2]*(1-x[, 4])*(1-x[, 4])) + (5*(1-x[, 2])*x[, 4]*x[, 6]) + 
    (6*(1-x[, 2])*x[, 4]*(1-x[, 6])) + (7*(1-x[, 2])*(1-x[, 4])*x[, 6]) + (8*(1-x[, 2])*(1-x[, 4])*(1-x[, 6])))
}

f6 <- function(x) {
  return((4*ifelse(x[, 1] > 1 & x[, 3] > 0, 1, 0)) + (4*ifelse(x[, 5] > 1 & x[, 7] > 1, 1, 0)) + (x[, 8]*x[, 9]))
}

f7 <- function(x) {
  return(0.5*(x[, 1]^2 + x[, 2] + x[, 3]^2 + x[, 4] + x[, 5]^2 + x[, 6] + x[, 7]^2 + x[, 8] + x[, 9]^2 -11))
}

f8 <- function(x) {
  (1/sqrt(2))*(f4(x) + (x[, 1] + x[, 3] + x[, 5] + x[, 7] + x[, 8] + x[, 9] - 2))
  return(x[, 7])
}

mu <- function(choice, x){
  if (choice == 1) 
    return(f4(x))
  else if (choice == 2) 
    return(f6(x))
  else if (choice == 3) 
    return(f7(x))
  else if (choice == 4) 
    return(f6(x))
  else if (choice == 5) 
    return(f8(x))
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
    return(f6(x))
  } else if (choice == 5) {
    return(f4(x))
  } else {
    stop("Invalid choice for tau")
  }
}


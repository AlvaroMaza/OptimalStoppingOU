library(goffda)


NegLogLikOU <- function (theta, delta, X) {
  n <- length(X)
  
  alpha <- exp(theta[1])
  mu  <- theta[2]
  sigma <- exp(theta[3])
  
  term1 <- -n/2 * log(sigma / (2 * alpha))
  term2 <- -0.5 * n * log1p(-exp(-2 * alpha * delta))
  term3 <- -(alpha / sigma) * sum(((X[2:n] - mu - (X[1:(n-1)] - mu) * exp(-alpha * delta))^2) / (1 - exp(-2 * alpha * delta)))
  
  log_likelihood_value <- term1 + term2 + term3
  return(- log_likelihood_value)
  
}


# Function to compute confidence intervals using Fisher information
compute_CI <- function(theta_hat, fisher_info, n, alpha = 0.05) {

  
  standard_errors <- sqrt(diag(solve(fisher_info))) / sqrt(n)
  
  z_critical <- qnorm(1 - alpha / 2)
  
  margin_of_error <- z_critical * standard_errors
  lower_bounds <- theta_hat - margin_of_error
  upper_bounds <- theta_hat + margin_of_error
  
  lower_bounds[1] <- exp(lower_bounds[1])
  theta_hat[1] <- exp(theta_hat[1])
  upper_bounds[1] <- exp(upper_bounds[1])
  lower_bounds[3] <- exp(lower_bounds[3])
  theta_hat[3] <- exp(theta_hat[3])
  upper_bounds[3] <- exp(upper_bounds[3])
  
  results <- data.frame(
    Estimate = theta_hat,
    Lower_CI = lower_bounds,
    Upper_CI = upper_bounds
  )
  return(results)
}


est.OU <- function (X, delta) {
  
  
  theta <- c(log(1), mean(X), log(1))
  
  n <- length(X)
  optim.results <- optim( par = theta, fn = NegLogLikOU, delta = delta, X = X, hessian=TRUE)
  theta_hat <- optim.results$par
  hessian <- optim.results$hessian
  
  alpha <- exp(theta_hat[1])
  mu  <- theta_hat[2]
  sigma <- exp(theta_hat[3])
  fisher_info<-hessian/n
  
  ci <- compute_CI(theta_hat, fisher_info, n)
  return( ci )
}


X <- r_ou(n = 1, t = seq(0, 1, len = 201), x0 = 10, mu = 15, sigma = 1, alpha = 4)$data
X <- as.vector(X)
delta <- 0.005

result <- est.OU(X, delta)






############### TEST CI ###############

# Set parameters
n_iterations <- 100
percentage_alpha <- numeric(n_iterations)
percentage_mu <- numeric(n_iterations)
percentage_sigma <- numeric(n_iterations)

T <- 10
delta <- 0.01

mu <- 15
sigma <- 1
alpha <- 4


# Loop through iterations
for (i in 1:n_iterations) {
  X <- r_ou(n = 1, t = seq(0, T, len = (T/delta + 1)), x0 = 15, mu = mu, sigma = sigma, alpha = alpha)$data
  X <- as.vector(X)

  result <- est.OU(X, delta)
  percentage_alpha[i] <- result[1, 2] < alpha & alpha < result[1, 3]
  percentage_mu[i] <- result[2, 2] < mu & mu < result[2, 3]
  percentage_sigma[i] <- result[3, 2] < sigma & sigma < result[3, 3]
}


# Calculate percentages
percentage_alpha <- mean(percentage_alpha)
percentage_mu <- mean(percentage_mu)
percentage_sigma <- mean(percentage_sigma)

print(paste("Percentage for alpha:", percentage_alpha))
print(paste("Percentage for mu:", percentage_mu))
print(paste("Percentage for sigma:", percentage_sigma))


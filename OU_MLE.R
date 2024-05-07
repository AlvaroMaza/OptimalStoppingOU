library(goffda)
library(GGally)



NegLogLikOU <- function (theta, delta, X) {
  n <- length(X)
  
  alpha <- exp(theta[1])
  mu  <- theta[2]
  sigma <- exp(theta[3])
  
  term1 <- -n/2 * log(sigma / (2 * alpha))
  term2 <- -0.5 * n * log1p(-exp(-2 * alpha * delta))
  term3 <- -(alpha / sigma) * sum(((X[-1] - mu - (X[-length(X)] - mu) *
                                      exp(-alpha * delta))^2) / (1 - exp(-2 * alpha * delta)))
  
  log_likelihood_value <- term1 + term2 + term3
  return(- log_likelihood_value)
  
}


# Function to compute confidence intervals using Fisher information
compute_CI <- function(theta_hat, fisher_info, n, alpha = 0.05) {

  
  standard_errors <- sqrt(diag(solve(fisher_info)))
  
  z_critical <- qnorm(1 - alpha / 2)
  
  margin_of_error <- z_critical * standard_errors / sqrt(n)
  lower_bounds <- theta_hat - margin_of_error
  upper_bounds <- theta_hat + margin_of_error
  
  lower_bounds[-2] <- exp(lower_bounds[-2])
  theta_hat[-2] <- exp(theta_hat[-2])
  upper_bounds[-2] <- exp(upper_bounds[-2])

  
  results <- data.frame(
    Estimate = theta_hat,
    Lower_CI = lower_bounds,
    Upper_CI = upper_bounds
  )
  return(results)
}


est_OU <- function (X, delta, theta0 = c(log(1), mean(X), log(1))) {
  
  
  n <- length(X)
  optim.results <- optim( par = theta0, fn = NegLogLikOU, delta = delta, X = X, hessian=TRUE)
  theta_hat <- optim.results$par
  hessian <- optim.results$hessian
  
  # alpha <- exp(theta_hat[1])
  alpha <- theta_hat[1]
  mu  <- theta_hat[2]
  # sigma <- exp(theta_hat[3])
  sigma <- theta_hat[3]
  
  fisher_info<-hessian/n
  
  ci <- compute_CI(theta_hat, fisher_info, n)
  return(list(parameters = c(alpha, mu, sigma), fisher_info = fisher_info), ci = ci)
}





############### PAIRWISE PLOT ###############

#Parameters
sample_size <- 1000
T <- 5
delta <- 0.02
n <-(T/delta + 1)
mu <- 1
sigma <- 1
alpha <- 1.4
theta <- c(log(alpha), mu, log(sigma))


#Estimation
data <- data.frame(alpha = numeric(length = sample_size),
                               mu = numeric(length = sample_size),
                               sigma = numeric(length = sample_size))
for (i in 1:sample_size) {
  X <- r_ou(n = 1, t = seq(0, T, len = (T/delta + 1)), x0 = 15, mu = mu, sigma = sigma, alpha = alpha)$data
  
  
  result <- est_OU(X, delta)
  R <- chol(result$fisher_info)
  
  point <- sqrt(n) * (R %*% (result$parameters - theta))
  data[i,] <- point
}

#Plot
density_with_normal <- function(data, mapping, ...) {
  p <- ggplot(data = data, mapping = mapping) +
    geom_density(color = "turquoise") + 
    stat_function(fun = dnorm, args = list(mean = 0, sd = 1)) +
    theme_minimal() +
    xlim(c(-6,6))
}

scatterplots <- function(data, mapping, ...) {
  p <- ggplot(data = data, mapping = mapping) +
    geom_point(color = "turquoise") +
    theme_minimal() +
    xlim(c(-6,6)) +
    ylim(c(-6,6))
}

ggpairs(data, 
        upper = list(continuous = wrap("cor", size = 4)),
        lower = list(continuous = scatterplots),
        diag = list(continuous = density_with_normal)
)





######### THEORETICAL FISHER INFORMATION #########
mu_variance <- function (n, alpha, sigma, delta) {
  return( 2*n*alpha/sigma^2 * (((1-exp(-alpha*delta))^2)/(1-exp(-2*alpha*delta))) )
}






############### MARGINAL CONFIDENCE INTERVALS ###############

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

  result <- est_OU(X, delta)
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

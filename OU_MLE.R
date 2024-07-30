# NegLogLikOU function
# 
# Description:
# This function computes the negative log-likelihood for an Ornstein-Uhlenbeck (OU) process given
# the parameters (theta), the time step (delta), and the observed data (X).
# 
# Parameters:
# - theta: A vector of parameters for the OU process where
#    - theta[1]: logarithm of the rate parameter alpha (log(alpha))
#    - theta[2]: mean parameter mu
#    - theta[3]: logarithm of the variance parameter sigma^2 (log(sigma^2))
# - delta: The time step between observations.
# - X: A numeric vector of observed data points.
# 
# Returns:
# - The negative log-likelihood value for the given parameters and data.
# 
# Usage:
# neg_loglik_value <- NegLogLikOU(theta = c(log_alpha, mu, log_sigma2), delta = delta, X = data)

NegLogLikOU <- function (theta, delta, X) {
  n <- length(X)
  
  alpha <- exp(theta[1])
  mu  <- theta[2]
  sigma2 <- exp(theta[3])
  
  term1 <- -n/2 * log(sigma2 / (2 * alpha))
  term2 <- -0.5 * n * log1p(-exp(-2 * alpha * delta))
  term3 <- -(alpha / sigma2) * sum(((X[-1] - mu - (X[-length(X)] - mu) *
                                       exp(-alpha * delta))^2) / (1 - exp(-2 * alpha * delta)))
  
  log_likelihood_value <- term1 + term2 + term3
  return(- log_likelihood_value)
  
}

# est_OU function
# 
# Description:
# This function estimates the parameters of an Ornstein-Uhlenbeck (OU) process by minimizing
# the negative log-likelihood using the `optim` function with the "L-BFGS-B" method. 
# It returns the estimated parameters.
# 
# Parameters:
# - X: A numeric vector of observed data points.
# - delta: The time step between observations.
# - theta0: Initial guesses for the parameters in the form of a vector where
#    - theta0[1]: initial guess for log(alpha)
#    - theta0[2]: initial guess for mu
#    - theta0[3]: initial guess for log(sigma^2)
#    Default value is c(log(1), mean(X), log(1)).
# 
# Returns:
# - A list containing the estimated parameters:
#    - alpha: estimated rate parameter alpha
#    - mu: estimated mean parameter mu
#    - sigma2: estimated variance parameter sigma^2
# 
# Usage:
# parameter_estimates <- est_OU(X = data, delta = delta, theta0 = initial_guesses)

est_OU <- function (X, delta, theta0 = c(log(1), mean(X), log(1))) {
  
  optim.results <- optim(par = theta0, fn = NegLogLikOU, delta = delta, X = X, method = "L-BFGS-B")
  theta_hat <- optim.results$par
  
  alpha <- exp(theta_hat[1])
  mu  <- theta_hat[2]
  sigma2 <- exp(theta_hat[3])
  
  return(list(alpha = alpha, mu = mu, sigma2 = sigma2))
}

####
## - The function "boundary" computes the optimal stopping boundary of the (exponentially) 
## discounted optimal stopping problem with finite horizon and gain function 
## G(x) = (strike - x)^+. The underlying process with parameters is a 
## time-dependent Ornstein-Uhlenbeck given by the SDE
##            dX_t = slope(t)(pull(t) - X_t) + vol(t)dW_t.
## - The boundary is computed by running a Picard iteration algorithm,
## that stops when the L2 distance between consecutive boundaries is less than
## tol, and solves the the free-boundary equation.
## - The boundary is comp uted at the time points provided in time_line
## - If errors == TRUE, a vector of the errors of the Picard scheme is provided
## alongside the boundary
####
boundary <- function (tol = 1e-3, strike = 0, time_line, discount = 0, 
                      boundary_kernel, slope, pull, vol, errors = FALSE) {
  
  N <- length(time_line)       # Partition length
  expiration <- time_line[N]        # Expiration date
  delta <- time_line[2:N] - time_line[1:(N-1)]  # Length step vector
  
  # Creating errors vector if required
  if (errors) er <- c()
  
  # Pre-allocating boundary
  bnd <- rep(min((slope(expiration) * pull(expiration) + discount * strike) /
                   (slope(expiration) + discount), strike), N)
  
  # Auxiliary pre-computations for marginal_mean and marginal_var
  f1 <- function(s) Vectorize(function(t) {
    exp(-integrate(slope, lower = t, upper = s, 
                   subdivisions = 10000, rel.tol = 1e-10)$value)
  })
  I1 <- f1(expiration)(time_line)
  f2 <- function(t) slope(t) * pull(t) * f1(expiration)(t)
  I2 <- sapply(time_line, function(s) {
    integrate(f2, lower = s, upper = expiration, 
              subdivisions = 10000, rel.tol = 1e-10)$value
  })
  f3 <- function(t) (f1(expiration)(t) * vol(t))^2
  I3 <- sapply(time_line, function(s) {
    integrate(f3, lower = s, upper = expiration, 
              subdivisions = 10000, rel.tol = 1e-10)$value
  })
  
  # Kernel definition
  boundary_kernel <- function(c1, c2, i1, x1, i2, x2) {
    
    # Compute the marginal mean
    marginal_mean <- (x1 * I1[i1] + (I2[i1] -  I2[i2])) / I1[i2]
    # Compute the marginal standard deviation
    marginal_sd <- sqrt((I3[i1] - I3[i2]) / I1[i2]^2)
    # Compute standardized values
    x2 <- (x2 - marginal_mean) / marginal_sd
    # Compute normal distribution and density
    normal_dist <- pnorm(x2, mean = 0, sd = 1, lower.tail = T)
    normal_dens <- dnorm(x2, mean = 0, sd = 1)
    # Evaluate Kernel
    K <- exp(-discount * (time_line[i2] - time_line[i1])) * 
      ((c1 - c2 * marginal_mean) * normal_dist + c2 * marginal_sd * normal_dens) 
    
    return(K)
    
  }
  
  # Boundary computation
  e <- 1  # error in the while loop
  # Fixed point algorithm
  while (e > tol) {
    
    bnd_old <- bnd
    
    # 
    for (i in (N - 1):1) {  
      
      # print(paste0("updating boundary at t_", i, " = ", time_line[i]))
      # Evaluate the kernel
      K1 <- boundary_kernel(c1 = strike, c2 = 1, 
                            i1 = i, x1 = bnd_old[i], 
                            i2 = N, x2 = strike)
      K2 <- boundary_kernel(c1 = discount * strike + slope(time_line[(i + 1):N]) * pull(time_line[(i + 1):N]),
                            c2 = discount + slope(time_line[(i + 1):N]),
                            i1 = i, x1 = bnd_old[i], 
                            i2 = (i + 1):N, x2 = bnd_old[(i + 1):N])
      
      # Update the boundary at t_present
      bnd[i] <- strike - K1 - sum(K2 * delta[i:(N - 1)])
      
      if(any(bnd[i] > strike)){
        print("Imposible: Boundary above the strike price")
      }
      
    }
    
    # absolute L2 error
    e <- sum((bnd - bnd_old)^2)
    #print(paste0("error: ", e))
    if (errors) er <- c(er, e)
    
  }
  
  if (errors) return(list(boundary = bnd, errors = er))
  
  return(bnd)
  
}



####
boundary_wrapper <- function (theta,  strike, expiration, partition_length, discount=0) {
  alpha <- theta['alpha']
  mu <- theta['mu']
  sigma2 <- theta['sigma2']
  

  time_line <- seq(0, expiration, l = partition_length)

  slope <- function(t) rep(alpha, length(t))
  pull <- function(t) rep(mu, length(t))
  vol <- function(t) rep(sigma2, length(t))
  
  bnd <- boundary(strike = strike, time_line = time_line,
                  discount = discount, boundary_kernel = boundary_kernel,
                  slope = slope, pull = pull, vol = vol, errors = TRUE)
  
  return(bnd$boundary)
}

infer_boundary <- function(X, delta, z_alpha = 0.1, partition_length = 100, strike = 0) {
  
  est_theta <- est_OU(X, delta)
  expiration <- partition_length * delta
  
  
  ## Uncomment to replicate GLD-SLV study
  strike <- est_theta$mu
  
  theta_list <- c(alpha = est_theta$alpha, mu = est_theta$mu, sigma2 = est_theta$sigma2)
  boundary <- boundary_wrapper(theta_list, strike = strike, partition_length = partition_length, expiration = expiration)
  
  h_alpha <- 1e-5 * abs(est_theta$alpha)
  h_mu <- 1e-5 * abs(est_theta$mu)
  h_sigma2 <- 1e-5 * abs(est_theta$sigma2)
  
  g_prime_theta1 <- (boundary_wrapper(c(alpha = est_theta$alpha + h_alpha, mu = est_theta$mu, sigma2 = est_theta$sigma2), strike = strike, partition_length = partition_length, expiration=expiration) -
                       boundary_wrapper(c(alpha = est_theta$alpha - h_alpha, mu = est_theta$mu, sigma2 = est_theta$sigma2), strike = strike, partition_length = partition_length, expiration=expiration)) / (2 * h_alpha)
  
  g_prime_theta2 <- (boundary_wrapper(c(alpha = est_theta$alpha, mu = est_theta$mu + h_mu, sigma2 = est_theta$sigma2), strike = strike, partition_length = partition_length, expiration=expiration) -
                       boundary_wrapper(c(alpha = est_theta$alpha, mu = est_theta$mu - h_mu, sigma2 = est_theta$sigma2), strike = strike, partition_length = partition_length, expiration=expiration)) / (2 * h_mu)
  
  g_prime_theta3 <- (boundary_wrapper(c(alpha = est_theta$alpha, mu = est_theta$mu, sigma2 = est_theta$sigma2 + h_sigma2), strike = strike, partition_length = partition_length, expiration=expiration) -
                       boundary_wrapper(c(alpha = est_theta$alpha, mu = est_theta$mu, sigma2 = est_theta$sigma2 - h_sigma2), strike = strike, partition_length = partition_length, expiration=expiration)) / (2 * h_sigma2)
  
  grad_g <- data.frame(
    g_prime_theta1 = g_prime_theta1,
    g_prime_theta2 = g_prime_theta2,
    g_prime_theta3 = g_prime_theta3
  )
  
  fisher_info <- fisher_info_gradient(est_theta$alpha, est_theta$mu, est_theta$sigma2,
                                      X, delta, 
                                      d_logL_dalpha, d_logL_dmu, d_logL_dsigma2)
  fisher_info_inv <- solve(fisher_info)
  
  std_dev <- numeric(length = nrow(grad_g))
  for (i in 1:nrow(grad_g)) {
    grad <- as.numeric(grad_g[i, ])
    std_dev[i] <- sqrt(t(grad) %*% fisher_info_inv %*% grad)/sqrt(length(X))
  }
  
  upper_bound <- boundary + qnorm(1 - z_alpha/2) * std_dev
  lower_bound <- boundary - qnorm(1 - z_alpha/2) * std_dev
  
  return(list(lower_bound = lower_bound, boundary_est = boundary, upper_bound = upper_bound, est_theta = est_theta))
}


infer_boundary_gold_silver <- function(X, delta, z_alpha=0.1, partition_length = 100, strike = 0){
  
  est_theta <- est_OU(X, delta)
  #print(est_theta)
  expiration <- partition_length * delta
  
  ## Uncomment to replicate GLD-SLV study
  strike <- est_theta$mu
  
  theta_list <- c(alpha = est_theta$alpha, mu = est_theta$mu, sigma2 = est_theta$sigma2)
  boundary <- boundary_wrapper(theta_list, strike = strike, partition_length = partition_length, expiration = expiration)
  
  ### Delta method ###
  
  # Numerical derivation
  h<- 0.01
  g_prime_theta1 <- (boundary_wrapper(c(alpha = est_theta$alpha + h, mu = est_theta$mu, sigma2 = est_theta$sigma2), strike = strike, partition_length = partition_length, expiration=expiration) - boundary) / h
  g_prime_theta2 <- (boundary_wrapper(c(alpha = est_theta$alpha, mu = est_theta$mu + h, sigma2 = est_theta$sigma2), strike = strike, partition_length = partition_length, expiration=expiration) - boundary) / h
  g_prime_theta3 <- (boundary_wrapper(c(alpha = est_theta$alpha, mu = est_theta$mu, sigma2 = est_theta$sigma2 + h), strike = strike, partition_length = partition_length, expiration=expiration) - boundary) / h
  
  
  grad_g <- data.frame(
    g_prime_theta1 = g_prime_theta1,
    g_prime_theta2 = g_prime_theta2,
    g_prime_theta3 = g_prime_theta3
  )
  
  #std_dev of the boundary function
  fisher_info <- fisher_info_gradient(est_theta$alpha, est_theta$mu, est_theta$sigma2,
                                      X, delta, 
                                      d_logL_dalpha, d_logL_dmu, d_logL_dsigma2)
  fisher_info_inv <- solve(fisher_info)
  
  std_dev <- numeric(length = nrow(grad_g))
  for (i in 1:nrow(grad_g)) {
    grad <- as.numeric(grad_g[i, ])
    std_dev[i] <- sqrt(t(grad) %*% fisher_info_inv %*% grad)/sqrt(length(X))
  }
  
  # Calculate the confidence intervals
  upper_bound <- boundary + qnorm(1-z_alpha/2) * std_dev
  lower_bound <- boundary - qnorm(1-z_alpha/2) * std_dev
  
  return(list(lower_bound = lower_bound, boundary_est = boundary, upper_bound = upper_bound, est_theta = est_theta))
}


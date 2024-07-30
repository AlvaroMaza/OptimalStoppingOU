#### First derivative functions ####

# Function to compute the first derivative of the log-likelihood with respect to alpha
d_logL_dalpha <- function(X, t, alpha, mu, sigma2, delta) {
  n <- length(X)
  gamma <- exp(-alpha * delta)
  

  Y_i <- X[t] - mu - (X[t - 1] - mu) * gamma
  Z_i <- X[t] - mu
  Z_i_minus_1 <- X[t - 1] - mu
  
  
  term1 <- delta * gamma^2 / (1 - gamma^2)
  term2 <- Y_i * (Y_i*(gamma^2 - 1) + 
           2 * alpha * delta * gamma * (Z_i * gamma - Z_i_minus_1))/ ((1 - gamma^2)^2)

  
  result <- 1 / (2 * alpha) - term1 + (1 / sigma2) * term2
  return(result)
}

# Function to compute the first derivative of the log-likelihood with respect to mu
d_logL_dmu <- function(X, t, alpha, mu, sigma2, delta) {
  n <- length(X)
  gamma <- exp(-alpha * delta)
  

  Y_i <- X[t] - mu - (X[t - 1] - mu) * gamma
  term <- Y_i * (1 - gamma) / (1 - gamma^2)

  
  result <- 2 * alpha * term / sigma2
  return(result)
}

# Function to compute the first derivative of the log-likelihood with respect to sigma^2
d_logL_dsigma2 <- function(X, t, alpha, mu, sigma2, delta) {
  gamma <- exp(-alpha * delta)
  Y_i <- X[t] - mu - (X[t - 1] - mu) * gamma

  result <- -1 / (2 * sigma2) + (alpha / sigma2^2) * Y_i^2/(1-gamma^2)
  return(result)
}

# Function to compute the Hessian matrix using the gradient
fisher_info_gradient <- function(alpha, mu, sigma2, X, delta, 
                                  d_logL_dalpha, d_logL_dmu, d_logL_dsigma2) {
  sum <- 0
  
  for (t in 2:length(X)){
    grad <- c(d_logL_dalpha(X, t, alpha, mu, sigma2, delta),
              d_logL_dmu(X, t, alpha, mu, sigma2, delta),
              d_logL_dsigma2(X, t, alpha, mu, sigma2, delta))
    
    sum <- sum + tcrossprod(grad)
    
  }

  return(sum/length(X))
}

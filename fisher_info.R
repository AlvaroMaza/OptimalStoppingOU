# Function to compute the Hessian matrix of the OU process
fisher_info_second_derivatives <- function(alpha, mu, sigma2, X, delta,
                               d2_logL_dalpha2, d2_logL_dmu2, d2_logL_dsigma2,
                               d2_logL_dalphadmu, d2_logL_dalphadsigma, d2_logL_dmudsigma) {
  
  # Compute the elements of the Hessian matrix
  H11 <- d2_logL_dalpha2(alpha, mu, sigma2, X, delta)
  H22 <- d2_logL_dmu2(alpha, mu, sigma2, X, delta)
  H33 <- d2_logL_dsigma2(alpha, mu, sigma2, X, delta)
  H12 <- d2_logL_dalphadmu(alpha, mu, sigma2, X, delta)
  H13 <- d2_logL_dalphadsigma(alpha, mu, sigma2, X, delta)
  H23 <- d2_logL_dmudsigma(alpha, mu, sigma2, X, delta)
  
  # Construct the Hessian matrix
  hessian_matrix <- matrix(c(H11, H12, H13,
                             H12, H22, H23,
                             H13, H23, H33), 
                           nrow = 3, ncol = 3, byrow = TRUE)
  
  return(-hessian_matrix/length(X))
}

#### Second derivative functions ####

# Function to compute the second derivative of the log-likelihood with respect to alpha
d2_logL_dalpha2 <- function(alpha, mu, sigma2, X, delta) {
  n <- length(X)
  gamma <- exp(-alpha * delta)

  
  sum_term1 <- 0
  sum_term2 <- 0
  sum_term3 <- 0
  sum_term4 <- 0
  sum_term5 <- 0
  sum_term6 <- 0
  sum_term7 <- 0
  sum_term8 <- 0
  
  for (i in 2:n) {
    Y_i <- X[i] - mu - (X[i - 1] - mu) * gamma
    Z_i_minus_1 <- X[i - 1] - mu
    
    sum_term1 <- sum_term1 + (delta^2 * gamma^2) / (1 - gamma^2)^2
    sum_term2 <- sum_term2 + (delta^2 * gamma^2 * Y_i^2) / (1 - gamma^2)^2
    sum_term3 <- sum_term3 + (delta^2 * gamma^4 * Y_i^2) / (1 - gamma^2)^3
    sum_term4 <- sum_term4 + (delta^2 * gamma * Z_i_minus_1 * Y_i) / (1 - gamma^2)
    sum_term5 <- sum_term5 + (delta^2 * gamma^3 * Z_i_minus_1 * Y_i) / (1 - gamma^2)^2
    sum_term6 <- sum_term6 + (delta^2 * gamma^2 * Z_i_minus_1^2) / (1 - gamma^2)
    sum_term7 <- sum_term7 + (delta * gamma^2 * Y_i^2) / (1 - gamma^2)^2
    sum_term8 <- sum_term8 + (delta^2 * gamma * Z_i_minus_1 * Y_i) / (1 - gamma^2)
  }
  
  result <- -n / (2 * alpha^2) +
    2 * sum_term1 -
    4 * alpha/sigma2 * sum_term2 -
    8 * alpha/sigma2 * sum_term3 +
    2 * alpha/sigma2 * sum_term4 +
    8 * alpha/sigma2 * sum_term5 -
    2 * alpha/sigma2 * sum_term6 +
    4/sigma2 * sum_term7 -
    4/sigma2 * sum_term8
  
  return(result)
}


# Function to compute the second derivative of the log-likelihood with respect to mu
d2_logL_dmu2 <- function(alpha, mu, sigma2, X, delta) {
  n <- length(X)
  
  result <- - (2*n*alpha)/sigma2 * (1-exp(-alpha*delta))^2/(1-exp(-2*alpha*delta))
  return(result)
}


# Function to compute the second derivative of the log-likelihood with respect to sigma^2
d2_logL_dsigma2 <- function(alpha, mu, sigma2, X, delta) {
  n <- length(X)

  sum_term <- 0
  for (i in 2:n) {
    sum_term <- sum_term + ((X[i] - mu - (X[i - 1] - mu) * exp(-alpha * delta))^2 /
                              (1 - exp(-2 * alpha * delta)))
  }
  
  result <- (n / (2 * sigma2^2)) - (2 * alpha / sigma2^3) * sum_term
  return(result)
}


# Function to compute the mixed second derivative of the log-likelihood with respect to alpha and mu
d2_logL_dalphadmu <- function(alpha, mu, sigma2, X, delta) {
  n <- length(X)
  gamma <- exp(-alpha*delta)
  
  sum_term_1 <- 0
  for (i in 2:n){
    Y_i <- (X[i] - mu - (X[i - 1] - mu) * gamma)
    
    sum_term_1  <- sum_term_1 + (Y_i * (1-gamma))/(1-gamma^2)
  }
  
  sum_term_2 <- 0
  for (i in 2:n){
    Y_i <- (X[i] - mu - (X[i - 1] - mu) * gamma)
    
    num <- delta*gamma*(1-gamma^2)*((X[i] - mu)*(1-gamma) + Y_i) - 2*delta*gamma^2*(1-gamma)*Y_i
    den <- (1-gamma^2)^2
    
    sum_term_2  <- sum_term_2 + num/den
  }
  
  result <- 2/sigma2 * sum_term_1 + 2*alpha/sigma2 * sum_term_2
  return(result)
}


# Function to compute the mixed second derivative of the log-likelihood with respect to alpha and sigma^2
d2_logL_dalphadsigma <- function(alpha, mu, sigma2, X, delta) {
  n <- length(X)
  gamma <- exp(-alpha*delta)
  
  
  sum_term <- 0
  for (i in 2:n){
    Y_i <- (X[i] - mu - (X[i - 1] - mu) * gamma)
            
    num <- Y_i * (Y_i*(gamma^2-1) + 2*alpha*delta*gamma*((X[i] - mu)*gamma - (X[i-1] - mu)))
    den <- (1-gamma^2)^2
    
    sum_term  <- sum_term + num/den

  }
  result <- - 1/(sigma2^2) * sum_term
  return(result)
}


# Function to compute the mixed second derivative of the log-likelihood with respect to alpha and sigma^2
d2_logL_dmudsigma <- function(alpha, mu, sigma2, X, delta) {
  n <- length(X)
  
  sum_term <- 0
  for (i in 2:n) {
    num <- (X[i] - mu - (X[i - 1] - mu) * exp(-alpha * delta)) * (1 - exp(-alpha * delta)) 
    den <- 1 - exp(-2 * alpha * delta)
    sum_term <- sum_term + num/den
  }
  result <- -2*alpha/(sigma2^2) * sum_term
  return(result)
}





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


######### THEORETICAL FISHER INFORMATION #########
mu_fisher <- function (alpha, sigma2, delta) {
  return( 2*alpha/sigma2 * (((1-exp(-alpha*delta))^2)/(1-exp(-2*alpha*delta))) )
}

sigma2_fisher <- function(alpha, sigma2, delta) {
  return( 1/(2*(sigma2^2)) )
}
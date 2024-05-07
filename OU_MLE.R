library(goffda)

negative_log_likelihood <- function(params, delta_t, X) {
  alpha <- params[1]
  sigma_sq <- params[2]
  mu <- params[3]
  
  n <- length(X)
  term1 <- -n/2 * log(sigma_sq / (2 * alpha))
  term2 <- -0.5 * n * log1p(-exp(-2 * alpha * delta_t))
  
  sum <- 0
  for (i in 2:n) {
    nominator <- (X[i] - mu - (X[i - 1] - mu) * exp(-alpha * delta_t))^2
    denominator <- 1 - exp(-2 * alpha * delta_t)
    sum <- sum + nominator/denominator
  }
  term3 <- -(alpha / sigma_sq) * sum
  
  log_likelihood_value <- term1 + term2 + term3
  print(c(term1, term2, term3))
  return(- log_likelihood_value)
}




delta_t <- 0.01 
T <- 1

initial_params <- c(1, 1, 1) 
lower_bounds <- c(0.05, 0.05, -Inf)  
upper_bounds <- c(Inf, Inf, Inf)    

result_list <- list()
num_iterations <- 50


for (i in 1:num_iterations) {
  X <- r_ou(n = 1, t = seq(0, T, len = T/delta_t), x0 = 20, mu = 20, sigma = 2, alpha = 4)$data
  
  result <- optim(method = "L-BFGS-B", par = initial_params, fn = negative_log_likelihood, delta_t = delta_t, X = X,
                  lower = lower_bounds, upper = upper_bounds)
  result_list[[i]] <- result$par
}


result_matrix <- do.call(rbind, result_list)
means <- colMeans(result_matrix)
print(means)




sd_val <- num_iterations/(2*(16)^2)
hist(result_matrix[, 2], main = paste("Histogram of sigma^2"),
     xlab = paste("Sigma^2"), breaks = 20, probability = T, ylim = c(0, 1))
curve(dnorm(x, mean = 16, sd = sd_val), add = TRUE, col = "blue")



hist(result_matrix[, 1], main = paste("Histogram of alpha"),
     xlab = paste("Alpha"), breaks = 20)



hist(result_matrix[, 3], main = paste("Histogram of mu"),
     xlab = paste("Mu"), breaks = 20)






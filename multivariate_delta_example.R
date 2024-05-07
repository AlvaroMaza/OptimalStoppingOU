library(mvtnorm)

g <- function(theta) {
  (sin(theta[1]) * theta[2]) / (theta[3] * 2 + 3)
}  
  
  
# Set parameters
n <- 10000  # Sample size
mu <- c(10, 10, 10)  # Mean 
sigma <- matrix(c(2, 0, 0, 0, 3, 0, 0, 0, 4), nrow = 3)  # Covariance matrix


# Generate random sample
set.seed(123)
samples <- rmvnorm(n, mu, sigma)
g_values <- apply(samples, 1, g)





##########GRADIENT NUMERICALLY##################
h <- 0.001

mean_theta1 <- mean(samples[,1])
mean_theta2 <- mean(samples[,2])
mean_theta3 <- mean(samples[,3])

num_g_prime_theta1 <- (g(c(mean_theta1 + h, mean_theta2, mean_theta3)) - g(c(mean_theta1, mean_theta2, mean_theta3))) / h
num_g_prime_theta2 <- (g(c(mean_theta1, mean_theta2 + h, mean_theta3)) - g(c(mean_theta1, mean_theta2, mean_theta3))) / h
num_g_prime_theta3 <- (g(c(mean_theta1, mean_theta2, mean_theta3 + h)) - g(c(mean_theta1, mean_theta2, mean_theta3))) / h


num_grad_g <- c(num_g_prime_theta1, num_g_prime_theta2, num_g_prime_theta3)


##########GRADIENT ANALITICALLY##################
g_prime_theta1 <- function(theta) {
  (cos(theta[1]) * theta[2]) / (theta[3] * 2 + 3)
}

g_prime_theta2 <- function(theta) {
  sin(theta[1]) / (theta[3] * 2 + 3)
}

g_prime_theta3 <- function(theta) {
  -2*(sin(theta[1]) * theta[2]) / (theta[3] * 2 + 3)^2
}

grad_g <- c(g_prime_theta1(mu), g_prime_theta2(mu), g_prime_theta3(mu))


##########GRADIENT WITH LIBRARY##################
grad_g <- grad(g, x = mu)






asymptotic_mean <- g(c(mu[1], mu[2], mu[3]))
asymptotic_variance <- t(grad_g) %*% sigma %*% grad_g

# Plot
hist(g_values, breaks = 30, freq = FALSE, col = "lightblue", main = "Distribution of g(theta1, theta2, theta3)",
     xlab = "g(theta1, theta2, theta3)")
curve(dnorm(x, mean = asymptotic_mean, sd =  sqrt(asymptotic_variance)), 
      col = "blue", lwd = 2, add = TRUE)



      
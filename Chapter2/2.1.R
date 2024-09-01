library(ggplot2)
library(goffda)

parameter_estimation_table <- function(mu, sigma2, alpha, sample_size = 10000, T = 1, delta = 0.01) {
  n <- T/delta
  theta <- c(alpha, mu, sigma2)
  
  # Estimation
  data <- data.frame(alpha = numeric(length = sample_size),
                     mu = numeric(length = sample_size),
                     sigma2 = numeric(length = sample_size))
  
  for (i in 1:sample_size) {
    X <- r_ou(n = 1, t = seq(0, T, len = n), x0 = 26, mu = mu, sigma = sqrt(sigma2), alpha = alpha)$data
    
    result <- est_OU(X, delta)
    
    data[i,] <- c(result$alpha, result$mu, result$sigma2)
  }
  
  return(data)
}

# Usage example:
results <- parameter_estimation_table(mu = 20, sigma2 = 2, alpha = 3)

# Calculate statistics
means <- colMeans(results)
medians <- apply(results, 2, median)
sds <- apply(results, 2, sd)

# Print results
cat("Mean of alpha: ", means['alpha'], "\n")
cat("Median of alpha: ", medians['alpha'], "\n")
cat("Standard deviation of alpha: ", sds['alpha'], "\n\n")

cat("Mean of mu: ", means['mu'], "\n")
cat("Median of mu: ", medians['mu'], "\n")
cat("Standard deviation of mu: ", sds['mu'], "\n\n")

cat("Mean of sigma2: ", means['sigma2'], "\n")
cat("Median of sigma2: ", medians['sigma2'], "\n")
cat("Standard deviation of sigma2: ", sds['sigma2'], "\n")


# Save the alpha density plot as a PDF
pdf("plots/alpha_density_plot.pdf", width = 12, height = 7)
ggplot(results, aes(x = alpha)) + 
  geom_density(fill = "lightblue", color = "black", alpha = 0.7) +
  geom_vline(xintercept = 3, color = "red", linetype = "dashed", size = 1) +
  xlab(expression(hat(alpha))) + 
  ylab("Density") +
  xlim(1, 5)
dev.off()

# Save the mu density plot as a PDF
pdf("plots/mu_density_plot.pdf", width = 12, height = 7)
ggplot(results, aes(x = mu)) + 
  geom_density(fill = "lightblue", color = "black", alpha = 0.7) +
  geom_vline(xintercept = 20, color = "red", linetype = "dashed", size = 1) +
  xlab(expression(hat(mu))) + ylab("Density") +
  xlim(18, 22)
dev.off()

# Save the sigma2 density plot as a PDF
pdf("plots/sigma2_density_plot.pdf", width = 12, height = 7)
ggplot(results, aes(x = sigma2)) + 
  geom_density(fill = "lightblue", color = "black", alpha = 0.7) +
  geom_vline(xintercept = 2, color = "red", linetype = "dashed", size = 1) +
  xlab(expression(hat(sigma)^2)) + ylab("Density") +
  xlim(1.4, 2.6)
dev.off()

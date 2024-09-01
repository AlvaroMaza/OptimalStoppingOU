# Load necessary library
library(ggplot2)

# Function to estimate parameters from saved paths
parameter_estimation_table_from_rds <- function(file_paths, delta = 0.01) {
  sample_size <- length(file_paths)
  results <- data.frame(alpha = numeric(sample_size),
                        mu = numeric(sample_size),
                        sigma2 = numeric(sample_size))
  
  for (i in 1:sample_size) {
    # Load the sample path from the .rds file
    X <- readRDS(file_paths[i])
    
    # Estimate the parameters using the est_OU function
    result <- est_OU(X, delta)

    # Store the estimates in the data frame
    results[i,] <- c(result$alpha, result$mu, result$sigma2)
  }
  
  return(results)
}

# Example usage: specify the path where your .rds files are stored
folder_name <- "rds_files" # Folder containing .rds files
file_paths <- list.files(path = folder_name, pattern = "X_path_.*\\.rds", full.names = TRUE)

# Estimate parameters using the saved paths
results <- parameter_estimation_table_from_rds(file_paths, delta = 0.01)

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


if (!dir.exists("scenario_plots")) {
  dir.create("scenario_plots")
}

# Save the alpha density plot as a PDF
pdf("scenario_plots/alpha_density_plot.pdf", width = 12, height = 7)
ggplot(results, aes(x = alpha)) + 
  geom_density(fill = "lightblue", color = "black", alpha = 0.7) +
  geom_vline(xintercept = 3, color = "red", linetype = "dashed", size = 1) +
  xlab(expression(hat(alpha))) + ylab("Density")
dev.off()

# Save the mu density plot as a PDF
pdf("scenario_plots/mu_density_plot.pdf", width = 12, height = 7)
ggplot(results, aes(x = mu)) + 
  geom_density(fill = "lightblue", color = "black", alpha = 0.7) +
  geom_vline(xintercept = 20, color = "red", linetype = "dashed", size = 1) +
  xlab(expression(hat(mu))) + ylab("Density")
dev.off()

# Save the sigma2 density plot as a PDF
pdf("scenario_plots/sigma2_density_plot.pdf", width = 12, height = 7)
ggplot(results, aes(x = sigma2)) + 
  geom_density(fill = "lightblue", color = "black", alpha = 0.7) +
  geom_vline(xintercept = 2, color = "red", linetype = "dashed", size = 1) +
  xlab(expression(hat(sigma)^2)) + ylab("Density")
dev.off()

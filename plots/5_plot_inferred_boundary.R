library(goffda)

# Function to plot the boundary
plot_boundary <- function(data, inferred_results, actual_boundary, split_index) {
  # Plot the data up to the split point
  plot(as.numeric(data[1:split_index]), type = "l", col = "black", xlim = c(0, length(data)), ylim = range(c(min(data), max(data))),
       xlab = "Index", ylab = "")
  # Add the data after the split as a dotted line
  lines((split_index + 1):length(data), data[(split_index + 1):length(data)], col = "black", lty = 2)
  
  # Add the inferred boundary estimates and confidence intervals
  lines(split_index:(split_index + length(inferred_results$boundary_est) - 1), inferred_results$boundary_est, col = "blue", lwd = 2)
  lines(split_index:(split_index + length(inferred_results$upper_bound) - 1), inferred_results$upper_bound, col = "blue", lty = 2)
  lines(split_index:(split_index + length(inferred_results$lower_bound) - 1), inferred_results$lower_bound, col = "blue", lty = 2)
  
  # Add a shaded area for the confidence interval (optional)
  polygon(c(split_index:(split_index + length(inferred_results$boundary_est) - 1), rev(split_index:(split_index + length(inferred_results$boundary_est) - 1))), 
          c(inferred_results$upper_bound, rev(inferred_results$lower_bound)), col = rgb(0, 0, 1, 0.2), border = NA)
  
  # Plot the actual boundary for comparison
  lines(split_index:(split_index + length(actual_boundary) - 1), actual_boundary, col = "red", lty = 1)
  
  # Add a vertical line at the split point
  abline(v = split_index, col = "black", lwd = 1)
}

# Generate data with len = 500
data <- r_ou(n = 1, t = seq(0, 5, len = 500), x0 = 18, mu = 20, sigma = sqrt(2), alpha = 3)$data

# Subset the first 400 data points
data_first_400 <- data[1:400]

# Infer boundary using the first 400 data points
inferred_results <- infer_boundary(data_first_400, delta = 0.01, z_alpha = 0.1)

# Define the split index for the boundary lines
split_index <- length(data_first_400)

# Actual boundary for comparison
actual_boundary <- boundary_wrapper(c(alpha = 3, mu = 20, sigma2 = 2))

# Plot the boundary
plot_boundary(data, inferred_results, actual_boundary, split_index)

# Optionally, save the plot as a PNG file
# png("plots/boundary_plot.png", width = 800, height = 600, res = 150)
# plot_boundary(data, inferred_results, actual_boundary, split_index)
# dev.off()



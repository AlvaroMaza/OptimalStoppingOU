# Load required packages
library(ggplot2)

# Function to plot the boundary with time as the x-axis
plot_boundary <- function(data, inferred_results, actual_boundary, split_index, timestep) {
  # Calculate time points
  time_points <- (0:(length(data) - 1)) * timestep
  split_time <- split_index * timestep
  
  # Create a plotting device with time as the x-axis
  plot(time_points[1:split_index], as.numeric(data[1:split_index]), type = "l", col = "black",
       xlim = c(0, max(time_points)), ylim = range(c(14, max(data))),
       xlab = "Time", ylab = "Value")
  
  # Add the data after the split as a dotted line
  lines(time_points[(split_index + 1):length(data)], data[(split_index + 1):length(data)], col = "black", lty = 2)
  
  # Add the inferred boundary estimates and confidence intervals
  lines(time_points[split_index:(split_index + length(inferred_results$boundary_est) - 1)], 
        inferred_results$boundary_est, col = "blue", lwd = 2)
  lines(time_points[split_index:(split_index + length(inferred_results$upper_bound) - 1)], 
        inferred_results$upper_bound, col = "blue", lty = 2)
  lines(time_points[split_index:(split_index + length(inferred_results$lower_bound) - 1)], 
        inferred_results$lower_bound, col = "blue", lty = 2)
  
  # Add a shaded area for the confidence interval
  polygon(c(time_points[split_index:(split_index + length(inferred_results$boundary_est) - 1)], 
            rev(time_points[split_index:(split_index + length(inferred_results$boundary_est) - 1)])), 
          c(inferred_results$upper_bound, rev(inferred_results$lower_bound)), col = rgb(0, 0, 1, 0.2), border = NA)
  
  # Plot the actual boundary for comparison
  lines(time_points[split_index:(split_index + length(actual_boundary) - 1)], actual_boundary, col = "red", lty = 1)
  
  # Add a vertical line at the split point
  abline(v = split_time, col = "black", lwd = 1)
}

# Load the specified file index
file_index <- 1  # Replace with the desired file index
folder_name <- "rds_files"

# Define file paths
x_path_file <- file.path(folder_name, paste0("X_path_", file_index, ".rds"))
inferred_boundary_file <- file.path(folder_name, paste0("inferred_boundary_", file_index, ".rds"))

# Load the data from the .rds files
data <- readRDS(x_path_file)
inferred_results <- readRDS(inferred_boundary_file)

# Define the split index and timestep
split_index <- 400  # Assuming the split is always at 400 points
timestep <- 0.01

# Define the actual boundary
actual_boundary <- boundary_wrapper(c(alpha = 3, mu = 20, sigma2 = 2), partition_length = 100, strike = 19.5, expiration = 1)


pdf_filename <- paste0("scenario_plots/boundary_plot_path_", file_index, ".pdf")
# Save the plot as a PDF
pdf(file = pdf_filename, width = 12, height = 7)
plot_boundary(data, inferred_results, actual_boundary, split_index, timestep)
dev.off()


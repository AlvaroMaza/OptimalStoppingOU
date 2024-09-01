# Load required packages
library(ggplot2)

# Function to plot the boundary with time as the x-axis
plot_boundary <- function(data, inferred_results, actual_boundary, split_index, timestep) {
  # Calculate time points
  time_points <- (0:(length(data) - 1)) * timestep
  split_time <- split_index * timestep
  
  # Create a plotting device with time as the x-axis
  plot(time_points[(split_index + 1):length(data)], as.numeric(data[(split_index + 1):length(data)]), type = "l", col = "black",
       xlim = c(4, max(time_points)), ylim = range(c(17, max(data[(split_index + 1):length(data)] + 1))),
       xlab = "Time", ylab = "Value")
  
  
  # Add the inferred boundary estimates and confidence intervals
  lines(time_points[(split_index+1):(split_index + length(inferred_results$boundary_est))], 
        inferred_results$boundary_est, col = "blue", lwd = 2)
  lines(time_points[(split_index+1):(split_index + length(inferred_results$upper_bound))], 
        inferred_results$upper_bound, col = "blue", lty = 2)
  lines(time_points[(split_index+1):(split_index + length(inferred_results$lower_bound))], 
        inferred_results$lower_bound, col = "blue", lty = 2)
  
  # Add a shaded area for the confidence interval
  polygon(c(time_points[(split_index+1):(split_index + length(inferred_results$boundary_est))], 
            rev(time_points[(split_index+1):(split_index + length(inferred_results$boundary_est))])), 
          c(inferred_results$upper_bound, rev(inferred_results$lower_bound)), col = rgb(0, 0, 1, 0.2), border = NA)
  
  # Find the first crossings
  first_cross_boundary <- which(!data[(split_index+1):length(data)] > inferred_results$boundary_est)[1] + split_index - 1
  first_cross_upper <- which(!data[(split_index+1):length(data)] > inferred_results$upper_bound)[1] + split_index - 1
  first_cross_lower <- which(!data[(split_index+1):length(data)] > inferred_results$lower_bound)[1] + split_index - 1
  
  # Add vertical lines at the first crossings
  if (!is.na(first_cross_boundary)) {
    segments(x0 = time_points[first_cross_boundary], y0 = 0, 
             x1 = time_points[first_cross_boundary], y1 = data[first_cross_boundary], 
             col = "green", lwd = 1, lty = 2)
  }
  if (!is.na(first_cross_upper)) {
    segments(x0 = time_points[first_cross_upper], y0 = 0, 
             x1 = time_points[first_cross_upper], y1 = data[first_cross_upper + 1], 
             col = "purple", lwd = 1, lty = 2)
  }
  if (!is.na(first_cross_lower)) {
    segments(x0 = time_points[first_cross_lower], y0 = 0, 
             x1 = time_points[first_cross_lower], y1 = data[first_cross_lower], 
             col = "orange", lwd = 1, lty = 2)
  }
  segments(x0 = time_points[length(time_points)], y0 = 0, 
           x1 = time_points[length(time_points)], y1 = data[length(data)], 
           col = "red", lwd = 1, lty = 2)
}

# Load the specified file index
file_index <- 4  # Replace with the desired file index
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


pdf_filename <- paste0("plots/boundary_plot_path_with_lines", ".pdf")
# Save the plot as a PDF
pdf(file = pdf_filename)
plot_boundary(data, inferred_results, actual_boundary, split_index, timestep)
dev.off()

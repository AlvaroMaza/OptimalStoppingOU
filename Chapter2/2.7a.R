library(ggplot2)

# Define the folder where the .rds files are stored
folder_name <- "rds_files"

# Initialize a list to store the inferred boundaries
boundary_estimates_list <- list()

# Get the list of all inferred boundary files in the folder
inferred_files <- list.files(folder_name, pattern = "inferred_boundary_.*\\.rds", full.names = TRUE)

# Loop through the files and read the inferred boundary estimates
for (file in inferred_files) {
  inferred_results <- readRDS(file)
  boundary_estimates_list[[file]] <- inferred_results$boundary_est
}

# Define the real boundary
real_boundary <- boundary_wrapper(c(alpha = 3, mu = 20, sigma2 = 2), partition_length = 100, strike = 19.5, expiration = 1)

# Define the timestep
timestep <- 0.01

# Calculate time points based on the timestep and length of the boundary estimates
time_points <- (0:(length(boundary_estimates_list[[1]]) - 1)) * timestep

# Define the output PDF filename
pdf_filename <- "scenario_plots/all_inferred_boundaries.pdf"

# Save the plot as a PDF with specified dimensions
pdf(file = pdf_filename, width = 10, height = 7)
plot(NULL, xlim = c(0, max(time_points)), ylim = range(sapply(boundary_estimates_list, range)),
     xlab = "Time", ylab = "Value")

# Add each boundary estimate to the plot
for (boundary_estimate in boundary_estimates_list) {
  lines(time_points, boundary_estimate, col = rgb(0, 0, 1, alpha = 0.1))
}

# Add the real boundary in red
lines(time_points, real_boundary, col = "red", lwd = 2)

# Close the PDF device
dev.off()


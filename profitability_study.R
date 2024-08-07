# Load necessary libraries
library(dplyr)

# Define the folder and number of files
folder_name <- "rds_files"
num_files <- 200  # Number of files to process

# Initialize vectors to store profits for each path
profits_expiration <- numeric(num_files)
profits_crossing <- numeric(num_files)
profits_lower_bound <- numeric(num_files)
profits_upper_bound <- numeric(num_files)

# Loop through each file and compute the profits
for (file_index in 1:num_files) {
  # Define file paths
  x_path_file <- file.path(folder_name, paste0("X_path_", file_index, ".rds"))
  inferred_boundary_file <- file.path(folder_name, paste0("inferred_boundary_", file_index, ".rds"))
  
  # Load the data from the .rds files
  data <- readRDS(x_path_file)
  inferred_results <- readRDS(inferred_boundary_file)
  
  # Extract the testing data and boundary estimates
  testing_data <- data[401:500]
  boundary_est <- inferred_results$boundary_est
  lower_bound <- inferred_results$lower_bound
  upper_bound <- inferred_results$upper_bound
  
  # Compute profits for waiting until expiration
  last_test_value <- testing_data[length(testing_data)]
  profits_expiration[file_index] <- max(0, 19.5 - last_test_value)
  
  # Check initial condition for crossing the estimated boundary
  initial_test_value <- testing_data[1]
  initial_boundary_value <- boundary_est[1]
  if (initial_test_value > initial_boundary_value) {
    for (i in 1:length(testing_data)) {
      if (testing_data[i] < boundary_est[i]) {
        profits_crossing[file_index] <- 19.5 - boundary_est[i]
        break
      }
    }
  }
  
  # Check initial condition for crossing the lower bound
  if (initial_test_value > lower_bound[1]) {
    for (i in 1:length(testing_data)) {
      if (testing_data[i] < lower_bound[i]) {
        profits_lower_bound[file_index] <- 19.5 - lower_bound[i]
        break
      }
    }
  }
  
  # Check initial condition for crossing the upper bound
  if (initial_test_value < upper_bound[1]) {
    for (i in 1:length(testing_data)) {
      if (testing_data[i] > upper_bound[i]) {
        profits_upper_bound[file_index] <- max(0, 19.5 - upper_bound[i])
        break
      }
    }
  }
}

# Summarize the results
# Compute mean profits for each strategy
mean_profit_expiration <- mean(profits_expiration)
mean_profit_crossing <- mean(profits_crossing)
mean_profit_lower_bound <- mean(profits_lower_bound)
mean_profit_upper_bound <- mean(profits_upper_bound)

# Summarize the results
cat("Mean Profit from Holding Until Expiration:", mean_profit_expiration, "\n")
cat("Mean Profit from Exiting on Crossing the Boundary:", mean_profit_crossing, "\n")
cat("Mean Profit from Exiting on Crossing the Lower Bound:", mean_profit_lower_bound, "\n")
cat("Mean Profit from Exiting on Crossing the Upper Bound:", mean_profit_upper_bound, "\n")

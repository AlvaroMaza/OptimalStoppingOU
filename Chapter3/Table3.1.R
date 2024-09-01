# Function to check if at least 90% of the real path is within the confidence intervals
check_path_inclusion <- function(M = 100) {
  folder_name <- "rds_files"
  
  # Initialize vector to store results
  inclusion_results <- vector("logical", M)
  real_boundary <- boundary_wrapper(c(alpha = 1, mu = 20, sigma2 = 2), partition_length = 100, strike = 19.5, expiration = 1)
  
  # Loop through each trial index
  for (trial_index in 1:M) {
    x_sample_file <- file.path(folder_name, paste0("X_path_", trial_index, ".rds"))
    inferred_boundary_file <- file.path(folder_name, paste0("inferred_boundary_", trial_index, ".rds"))
    
    # Load the data
    X_sample <- readRDS(x_sample_file)
    inferred_results <- readRDS(inferred_boundary_file)
    
    # Count points within the confidence intervals
    inside_count <- sum(real_boundary >= inferred_results$lower_bound & real_boundary <= inferred_results$upper_bound)
    # Calculate the proportion of points inside the confidence intervals
    proportion_inside <- inside_count / length(real_boundary)
    
    # Check if at least 90% of the real path is inside the confidence intervals
    inclusion_results[trial_index] <- proportion_inside >= 0.9
  }
  
  # Calculate the proportion of paths with at least 90% inclusion
  proportion_paths_with_90_inclusion <- mean(inclusion_results)
  
  return(proportion_paths_with_90_inclusion)
}

# Set parameters
M <- 200

# Check path inclusion from saved files
proportion_90_inclusion <- check_path_inclusion(M = M)
print(proportion_90_inclusion)

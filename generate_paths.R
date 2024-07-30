library(parallel)
library(goffda)

compute_proportion_parallel <- function(M = 100) {
  # Set up the cluster with the number of cores to use
  num_cores <- detectCores() - 1
  cl <- makeCluster(num_cores)
  
  # Export necessary functions and objects to the cluster
  clusterExport(cl, c("r_ou", "infer_boundary", "boundary_wrapper", "boundary",
                      "est_OU", "NegLogLikOU", "fisher_info_gradient",
                      "d_logL_dalpha", "d_logL_dmu", "d_logL_dsigma2", "qnorm"))
  
  # Define the function to run in parallel
  process_trial <- function(trial_index) {
    
    # Generate sample paths
    X_sample <- r_ou(n = 1, t = seq(0, 5, len = 500), x0 = 18, mu = 20, sigma = sqrt(2), alpha = 3)$data
    
    # Split the sample
    X_training <- X_sample[1:400]
    
    # Infer boundary using the first half of the data
    results_first_half <- infer_boundary(X_training, delta = 0.01, z_alpha = 0.1, strike = 19.5)
    
    # Define the folder and file names for saving
    folder_name <- "rds_files"
    dir.create(folder_name, showWarnings = FALSE)  # Create the folder if it doesn't exist
    x_sample_file <- file.path(folder_name, paste0("X_path_", trial_index, ".rds"))
    inferred_boundary_file <- file.path(folder_name, paste0("inferred_boundary_", trial_index, ".rds"))
    
    # Save data to .rds files
    saveRDS(X_sample, file = x_sample_file)
    saveRDS(results_first_half, file = inferred_boundary_file)
    
    return(NULL)
  }
  
  # Run the trials in parallel
  parLapply(cl, 1:M, process_trial)
  
  # Stop the cluster
  stopCluster(cl)
}



# Compute proportions in parallel
compute_proportion_parallel(M=1000)




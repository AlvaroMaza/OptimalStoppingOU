# Define a function to compute the proportion of non-inclusions
compute_proportion <- function(z_alpha, M = 20) {
  proportions <- numeric(length = 200/2)
  for (i in 1:M) {
    print(i)
    # Generate sample paths
    X_sample <- r_ou(n = 1, t = seq(0, 5, len = 500), x0 = 18, mu = 20, sigma = sqrt(2), alpha = 3)$data
    
    # Split the sample
    X_training <- X_sample[1:400]
    
    # Infer boundary using the first half of the data
    results_first_half <- infer_boundary(X_training, delta = 0.01, z_alpha = z_alpha)
    
    # Calculate real boundary
    real_boundary <- boundary_wrapper(c(alpha = 3, mu = 20, sigma2 = 2))
    
    # Check if true boundary belongs to the confidence interval pointwise
    for (j in 1:length(real_boundary)) {
      if (real_boundary[j] < results_first_half$lower_bound[j] || real_boundary[j] > results_first_half$upper_bound[j]) {
        proportions[j] <- proportions[j] + 1
      }
    }
    
    # Plot the data
    plot(as.numeric(X_sample[1:start_index]), type = "l", col = "black", xlim = c(0, length(X_sample)), ylim = range(c(results_first_half$lower_bound, max(X_sample))),
         xlab = "Index", ylab = "Value", main = "Barrier Boundary with Confidence Intervals")
    lines((start_index + 1):length(X_sample), X_sample[(start_index + 1):length(X_sample)], col = "black", lty = 2)
    lines(start_index:(start_index + length(results_first_half$boundary_est) - 1), results_first_half$boundary_est, col = "blue", lwd = 2)
    lines(start_index:(start_index + length(results_first_half$upper_bound) - 1), results_first_half$upper_bound, col = "blue", lty = 2)
    lines(start_index:(start_index + length(results_first_half$lower_bound) - 1), results_first_half$lower_bound, col = "blue", lty = 2)
    polygon(c(start_index:(start_index + length(results_first_half$boundary_est) - 1), rev(start_index:(start_index + length(results_first_half$boundary_est) - 1))),
            c(results_first_half$upper_bound, rev(results_first_half$lower_bound)), col = rgb(0, 0, 1, 0.2), border = NA)
    lines(start_index:(start_index + length(real_boundary) - 1), real_boundary, col = "red", lty = 1)
    abline(v = start_index, col = "black", lwd = 1)
  }
  return(proportions/M)
}

# Calculate alpha
alpha <- 0.1
M <- 50
# Compute proportions
proportions <- compute_proportion(z_alpha=alpha, M = M)

# Calculate q_alpha
q_alpha <- qnorm(1 - alpha/2)

# Calculate dotted lines values
dotted_lines <- alpha + c(-1, 1) * q_alpha * sqrt(alpha * (1 - alpha) / M)

# Plot proportions
plot(proportions, type = "l", col = "red", ylim = c(0, 0.5),
     xlab = "Time", ylab = "Proportion of Non-inclusions",
     main = "Proportion of Trials with True Boundary Outside Confidence Interval")

# Add dashed line for alpha
abline(h = alpha, col = "black", lty = 2)

# Add dotted lines
abline(h = dotted_lines, col = "black", lty = 3)


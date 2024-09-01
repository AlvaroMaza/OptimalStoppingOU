library(optimx)

# Define the objective function for optimization (negative log-likelihood)
objective_function <- function(B, gold, silver, delta = 1/252) {
  # Compute the OU path for a given B
  OU_path <- 1/gold[1] * gold - B/silver[1] * silver
  
  # Estimate the OU parameters
  est_params <- est_OU(X = OU_path, delta = delta)
  
  # Calculate the negative log-likelihood
  neg_log_likelihood <- NegLogLikOU(theta = c(log(est_params$alpha), est_params$mu, log(est_params$sigma2)), 
                                    delta = delta, X = OU_path)
  
  return(neg_log_likelihood)
}

# Load data
silver_data <- read.csv('SLV.csv')
gold_data <- read.csv('GLD.csv')

# Extract the 'Open' prices
silver <- silver_data$Open[1:1100]
gold <- gold_data$Open[1:1100]
date_line <- silver_data$Date[1:1100]
date_line <- as.Date(date_line)

# Set the initial guess for B and range for optimization
initial_B <- 0.5
B_range <- c(0, 1)  # Define a reasonable range for B

# Optimize B using the optimize function
optimal_B <- optimize(f = objective_function, interval = B_range, 
                      gold = gold, silver = silver, delta = 1/252)

# Display the optimal B and the corresponding minimum negative log-likelihood
optimal_B_value <- optimal_B$minimum
optimal_neg_log_likelihood <- optimal_B$objective

cat("Optimal B:", optimal_B_value, "\n")
cat("Minimum Negative Log-Likelihood:", optimal_neg_log_likelihood, "\n")

# Optimal B obtained from the previous optimization step
optimal_B <- 0.5944895

# Compute the optimal OU path
OU_path <- 1/gold[1] * gold - optimal_B/silver[1] * silver

pdf("silver_gold_OU.pdf", width = 12, height = 7)
plot(date_line, OU_path, type = "l",
     xlab = "Date", ylab = "Value")
dev.off()







## Plot of the likelihood

# Generate a sequence of B values over the range
B_values <- seq(0, 1, by = 0.01)

# Compute the negative log-likelihood for each B value
neg_log_likelihood_values <- sapply(B_values, function(B) {
  objective_function(B, gold, silver, delta = 1/252)
})


pdf("likelihood.pdf", width = 12, height = 7)
# Plot the negative log-likelihood as a function of B
plot(B_values, -neg_log_likelihood_values, type = "l", col = "black", lwd = 2,
     xlab = "B", ylab = "Total Log-Likelihood")

# Add a point for the optimal B
points(optimal_B_value, -optimal_neg_log_likelihood, col = "red", pch = 19)

# Add a dashed vertical red line from the x-axis to the optimal point
segments(x0 = optimal_B_value, y0 = par("usr")[3], x1 = optimal_B_value, y1 = -optimal_neg_log_likelihood, col = "red", lty = 2)
dev.off()








# Function to plot the inferred boundary
plot_boundary_gold_silver <- function(data, inferred_results, split_index, date_line) {

  
  # Plot the data up to the split point
  plot(date_line[1:split_index], as.numeric(data[1:split_index]), type = "l", col = "black", 
       xlim = range(date_line), ylim = range(c(min(data), max(data))),
       xlab = "Date", ylab = "Value")
  
  # Add the data after the split as a dotted line
  lines(date_line[(split_index + 1):length(data)], data[(split_index + 1):length(data)], col = "black", lty = 2)
  
  # Add the inferred boundary estimates and confidence intervals
  boundary_dates <- date_line[split_index:(split_index + length(inferred_results$boundary_est) - 1)]
  lines(boundary_dates, inferred_results$boundary_est, col = "blue", lwd = 2)
  lines(boundary_dates, inferred_results$upper_bound, col = "blue", lty = 2)
  lines(boundary_dates, inferred_results$lower_bound, col = "blue", lty = 2)
  
  # Add a shaded area for the confidence interval (optional)
  polygon(c(boundary_dates, rev(boundary_dates)), 
          c(inferred_results$upper_bound, rev(inferred_results$lower_bound)), 
          col = rgb(0, 0, 1, 0.2), border = NA)
  
  # Add a vertical line at the split point
  abline(v = date_line[split_index], col = "black", lwd = 1)
  
  # Add the shaded area from the 400th observation to the split point
  shade_start <- split_index - 200
  rect(xleft = date_line[shade_start + 1], ybottom = min(data), 
       xright = date_line[split_index], ytop = max(data), 
       col = rgb(0.6, 0.6, 1, 0.2), border = NA)
}



# Initialize variables to accumulate total profits for each strategy
total_profit_crossing <- 0
#total_profit_crossing_upper_bound <- 0
total_profit_crossing_lower_bound <- 0
total_profit_expiration <- 0

# Initialize vectors to store profits for each iteration
profit_crossing_values <- c()
profit_expiration_values <- c()

testing_length <- 60
# Loop to update start variable and plot results
for (start in seq(0, length(OU_path) - 260, by = testing_length)) {
  
  # Define training and testing length
  training_length <- start + 200
  
  # Extract training and testing data
  training_data <- OU_path[start + 1:(training_length - start)]
  testing_data <- OU_path[(training_length + 1):(training_length + testing_length)]
  
  # Infer the boundary using the training data
  inferred_results <- infer_boundary_gold_silver(training_data, delta = 1/252, z_alpha = 0.1, partition_length = testing_length)
  #print(inferred_results$boundary_est)
  # Define the split index
  split_index <- training_length
  
  # Create a filename for each plot
  plot_filename <- paste0("gold_silver_", start, ".pdf")
  
  # Open a PDF device
  pdf(plot_filename, width = 12, height = 7)
  
  # Plot the boundary
  plot_boundary_gold_silver(OU_path, inferred_results, split_index, date_line)
  
  # Close the PDF device
  dev.off()
  
  # Check the relationship between the first testing data point and the inferred boundary
  initial_test_value <- testing_data[1]
  initial_boundary_value <- inferred_results$boundary_est[1]
  
  if (initial_test_value > initial_boundary_value) {
    # Strategy 1: Exit on Crossing
    for (i in 1:length(testing_data)) {
      if (testing_data[i] < inferred_results$boundary_est[i]) {
        profit_crossing <- inferred_results$est_theta$mu - inferred_results$boundary_est[i]
        total_profit_crossing <- total_profit_crossing + profit_crossing
        break
      }
    }
    
    # Strategy 2: Exit on Crossing lower bound
    for (i in 1:length(testing_data)) {
      if (testing_data[i] < inferred_results$lower_bound[i]) {
        profit_crossing_lower <- inferred_results$est_theta$mu - inferred_results$lower_bound[i]
        total_profit_crossing_lower_bound <- total_profit_crossing_lower_bound + profit_crossing_lower
        break
      }
    }
    
    # Strategy 3: Hold Until Expiration
    last_test_value <- testing_data[length(testing_data)]
    profit_expiration <- max(0,inferred_results$est_theta$mu - last_test_value)
    total_profit_expiration <- total_profit_expiration + profit_expiration
    
  }
  # Append profits to the vectors
  profit_crossing_values <- c(profit_crossing_values, total_profit_crossing)
  profit_expiration_values <- c(profit_expiration_values, total_profit_expiration)
  cat("Profit crossing:", total_profit_crossing, ", Profit expiration:", total_profit_expiration, "\n")
}

cat("Total Profit from Exit on Crossing:", total_profit_crossing, "\n")
cat("Total Profit from Exit on Crossing Lower Bound:", total_profit_crossing_lower_bound, "\n")
cat("Total Profit from Holding Until Expiration:", total_profit_expiration, "\n")


profit_crossing_values <- c(rep(0, 5), profit_crossing_values)
profit_expiration_values <- c(rep(0, 5), profit_expiration_values)





pdf("profit_over_time.pdf", width = 12, height = 7)
# Plot OU_path over time with the primary y-axis
plot(date_line, OU_path, type = "l", col = "black", lwd = 2,
     xlab = "Date", ylab = "OU Path Value")

# Create a secondary y-axis for the profit
par(new = TRUE)

# Plot the profit data using the secondary y-axis
plot(profit_crossing_values, type = "l", col = "blue", lwd = 2,
     xaxt = "n", yaxt = "n", xlab = "", ylab = "", ylim = range(c(profit_crossing_values, profit_expiration_values)))

# Add the profit expiration line
lines(profit_expiration_values, type = "l", col = "red", lwd = 2)

# Add secondary y-axis labels
axis(4)

dev.off()

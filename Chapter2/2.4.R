library(goffda)
library(ggplot2)
library(plotly)

check_convergence <- function(mu, sigma2, alpha, T = 4, delta = 0.01, M = 1000, alpha_level = 0.05) {
  n <- T / delta
  theta <- c(alpha, mu, sigma2)
  K <- length(theta)
  
  chi_squared_stats <- numeric(M)
  vectors <- matrix(NA, nrow = M, ncol = K)
  
  for (i in 1:M) {
    X <- r_ou(n = 1, t = seq(0, T, len = n), x0 = 28, mu = mu, sigma = sqrt(sigma2), alpha = alpha)$data
    
    result <- est_OU(X, delta)
    
    fisher_info <- fisher_info_second_derivatives(result$alpha, result$mu, result$sigma2, X, delta,
                                                  d2_logL_dalpha2, d2_logL_dmu2, d2_logL_dsigma2,
                                                  d2_logL_dalphadmu, d2_logL_dalphadsigma, d2_logL_dmudsigma)
    
    R <- chol(fisher_info)
    
    z <- sqrt(n) * (R %*% (c(result$alpha, result$mu, result$sigma2) - theta))
    chi_squared_stats[i] <- sum(z^2)
    vectors[i, ] <- z
  }
  
  chi_squared_critical_value <- qchisq(1 - alpha_level, df = K)
  coverage_probability <- mean(chi_squared_stats <= chi_squared_critical_value)
  
  return(list(vectors = vectors, chi_squared_stats = chi_squared_stats, 
              critical_value = chi_squared_critical_value, 
              coverage_probability = coverage_probability))
}

# Parameters
mu <- 20
sigma2 <- 2
alpha <- 3
M <- 1000
alpha_level <- 0.1



# Check convergence for n=1000
convergence_results <- check_convergence(mu, sigma2, alpha, M = M, alpha_level = alpha_level)

coverage_probability <- convergence_results$coverage_probability
vectors <- convergence_results$vectors
critical_value <- convergence_results$critical_value
chi_squared_stats <- convergence_results$chi_squared_stats

# Calculate the critical radius
critical_radius <- sqrt(critical_value)

# Function to generate sphere coordinates
generate_sphere <- function(radius, n = 100) {
  theta <- seq(0, 2 * pi, length.out = n)
  phi <- seq(0, pi, length.out = n)
  theta_grid <- rep(theta, each = length(phi))
  phi_grid <- rep(phi, length(theta))
  
  x <- radius * sin(phi_grid) * cos(theta_grid)
  y <- radius * sin(phi_grid) * sin(theta_grid)
  z <- radius * cos(phi_grid)
  
  return(list(x = x, y = y, z = z))
}

# Generate sphere coordinates
sphere_coords <- generate_sphere(critical_radius)

# Create a 3D scatter plot with Plotly
fig <- plot_ly(x = vectors[,1], y = vectors[,2], z = vectors[,3], 
               type = 'scatter3d', mode = 'markers',
               marker = list(color = ifelse(chi_squared_stats > critical_value, 'red', 'blue'),
                             size = 3)) %>%
  add_trace(x = sphere_coords$x, y = sphere_coords$y, z = sphere_coords$z,
            type = 'scatter3d', mode = 'markers', 
            marker = list(color = 'lightblue', size = 2, opacity = 0.3),
            showlegend = FALSE) %>%
  layout(scene = list(xaxis = list(title = 'Z1'),
                      yaxis = list(title = 'Z2'),
                      zaxis = list(title = 'Z3')))

# Show the plot
fig
#kaleido(fig, "C:/Users/alvar/Desktop/OptimalStoppingOUz_j_plot.pdf")


# Calculate the critical value for the given confidence level
z_alpha <- qnorm(1 - alpha_level / 2)

# Calculate the standard error
standard_error <- sqrt(alpha_level * (1 - alpha_level) / M)

# Compute the confidence interval
lower_bound <- (1-alpha_level) - z_alpha * standard_error
upper_bound <- (1-alpha_level) + z_alpha * standard_error

cat("Observed coverage probability:", coverage_probability, "\n")
cat("90% Confidence interval:", lower_bound, "-", upper_bound, "\n")



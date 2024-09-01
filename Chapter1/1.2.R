library(ggplot2)
library(manipulate)

# Function to calculate transition probability density for Ornstein-Uhlenbeck process
transition_prob_density_OU <- function(x, x0, t, mu, sigma, theta) {
  density <- 1 / sqrt(2 * pi * (sigma^2 / (2 * theta)) * (1 - exp(-2 * theta * t))) * 
    exp(-(x - mu - (x0 - mu) * exp(- theta * t))^2 / ((sigma^2 / theta) * (1 - exp(-2 * theta * t))))
  return(density)
}

# Function to generate data for heatmap for Ornstein-Uhlenbeck t.p.d.
generate_heatmap_data_OU <- function(x0, sigma, theta, mu) {
  x_values <- seq(4.5, 11.5, length.out = 1000)
  t_values <- seq(0.1, 1, by = 0.001)
  
  data <- expand.grid(x = x_values, t = t_values)
  data$density <- with(data, transition_prob_density_OU(x, x0, t, mu, sigma, theta))
  
  return(data)
}

# Uncomment the following code to use manipulate for interactive plotting
# manipulate(
#   {
#     data <- generate_heatmap_data_OU(x0, sigma, theta, mu)
#     ggplot(data, aes(y = x, x = t, fill = density)) +
#       geom_tile() +
#       scale_fill_gradient(low = "blue", high = "yellow") +
#       labs(y = expression(x[t]), x = "t", fill = "Density") +
#       theme_minimal()
#   },
#   x0 = slider(-10, 10, initial = 5, step = 0.1, label = "X_0"),
#   sigma = slider(0.1, 2, initial = 2, step = 0.1, label = "Sigma"),
#   theta = slider(0.1, 3, initial = 3, step = 0.1, label = "Theta"),
#   mu = slider(-10, 10, initial = 10, step = 0.1, label = "Mu")
# )

# Set parameters for saving plots
x0 <- 5
sigma <- 2
theta <- 3
mu <- 10

# Generate data for heatmap
data <- generate_heatmap_data_OU(x0, sigma, theta, mu)

# Save the plot as a PDF
pdf("ou_tpd.pdf", width = 12, height = 7)
ggplot(data, aes(y = x, x = t, fill = density)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "yellow") +
  labs(y = "X(t)", x = "t", fill = "Density") +
  theme_minimal()
dev.off()

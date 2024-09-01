# Load package
library(goffda)
library(manipulate)

# Custom color palette
custom_palette <- colorRampPalette(c("lightblue", "blue"))

# Uncomment the following code to use manipulate for interactive plotting
# manipulate(
#   {
#     BM_data <- r_ou(n = 10, x0 = x0, mu = mu, sigma = sigma, alpha = alpha)
#     plot(BM_data, main = "Ornstein-Uhlenbeck Paths", col = custom_palette(5))
#   },
#   x0 = slider(0, 100, initial = 5, step = 0.1, label = "X_0"),
#   sigma = slider(0.1, 10, initial = 2, step = 0.1, label = "Sigma"),
#   alpha = slider(0.1, 10, initial = 3, step = 0.1, label = "Alpha"),
#   mu = slider(-10, 10, initial = 10, step = 0.1, label = "Mu")
# )

# Set parameters for saving plots
x0 <- 5
sigma <- 2
alpha <- 3
mu <- 10

# Generate Ornstein-Uhlenbeck Paths
BM_data <- r_ou(n = 10, x0 = x0, mu = mu, sigma = sigma, alpha = alpha)

# Save the plot as a PDF
pdf("ou_paths.pdf", width = 12, height = 7)
plot(BM_data, main = "", col = custom_palette(5))
dev.off()

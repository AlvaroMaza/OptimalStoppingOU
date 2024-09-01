# Define parameters
K <- 100  # Strike price of the put option

# Create a sequence of possible stock prices at expiration
stock_prices <- seq(50, 150, by = 1)

# Calculate the payoff of the put option at expiration
put_payoff <- pmax(K - stock_prices, 0) - 10

pdf("put.pdf", width=12, height = 7)
# Plot the payoff diagram
plot(stock_prices, put_payoff, type = "l", col = "blue", lwd = 2,
     xlab = "Stock Price at Expiration", ylab = "Profit")

# Add a line for the x-axis
abline(h = 0, col = "black")
dev.off()




# Calculate the payoff of the call option at expiration
call_payoff <- pmax(stock_prices - K, 0) - 10

pdf("call.pdf", width=12, height = 7)
# Plot the payoff diagram
plot(stock_prices, call_payoff, type = "l", col = "blue", lwd = 2,
     xlab = "Stock Price at Expiration", ylab = "Profit")

# Add a line for the x-axis
abline(h = 0, col = "black")
dev.off()

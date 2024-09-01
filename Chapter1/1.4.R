library(ggplot2)
library(readr)

nvidia_data <- read.csv("NVDA.csv")

# Step 3: Convert the Date column to Date type if it's not already
nvidia_data$Date <- as.Date(nvidia_data$Date, format = "%Y-%m-%d")

# Define the cutoff date
august_01 <- as.Date("2024-08-01")

pdf("nvda.pdf", width = 12, height = 7)
# Step 4: Plot the "Open" prices using base R
plot(nvidia_data$Date, nvidia_data$Open, type = "n", col = "black", lwd = 2,
     xlab = "Date", ylab = "Open Price")

# Add a solid line for dates before August 1
lines(nvidia_data$Date[nvidia_data$Date <= august_01], 
      nvidia_data$Open[nvidia_data$Date <= august_01], 
      col = "black", lwd = 2)

# Add a dashed line for dates from August 1 onwards
lines(nvidia_data$Date[nvidia_data$Date >= august_01], 
      nvidia_data$Open[nvidia_data$Date >= august_01], 
      col = "black", lwd = 2, lty = 2)

# Step 5: Add a vertical line at August 1st
abline(v = august_01, col = "red", lwd = 1, lty = 2)
dev.off()

library(goffda)
library(GGally)
library(ggplot2)
library(gridExtra)
library(grid)

############### PAIRWISE PLOT FUNCTION ###############

pairwise_plot <- function(mu, sigma2, alpha, T = 1, delta = 0.01) {
  n <- T/delta
  theta <- c(alpha, mu, sigma2)
  
  # Generate data for sample_size = 100
  data_100 <- data.frame(alpha = numeric(length = 100),
                         mu = numeric(length = 100),
                         sigma2 = numeric(length = 100))
  for (i in 1:100) {
    X <- r_ou(n = 1, t = seq(0, T, len = n), x0 = 28, mu = mu, sigma = sqrt(sigma2), alpha = alpha)$data
    
    result <- est_OU(X, delta)
    
    fisher_info <- fisher_info_second_derivatives(result$alpha, result$mu, result$sigma2, X, delta,
                                                  d2_logL_dalpha2, d2_logL_dmu2, d2_logL_dsigma2,
                                                  d2_logL_dalphadmu, d2_logL_dalphadsigma, d2_logL_dmudsigma)
    
    # fisher_info_gradient <- fisher_info_gradient(result$alpha, result$mu, result$sigma2, X, delta,
    #                                              d_logL_dalpha, d_logL_dmu, d_logL_dsigma2)
    
    R <- chol(fisher_info)
    
    point <- sqrt(n) * (R %*% (c(result$alpha, result$mu, result$sigma2) - theta))
    data_100[i,] <- point
  }
  
  # Generate data for sample_size = 1000
  data_1000 <- data.frame(alpha = numeric(length = 1000),
                          mu = numeric(length = 1000),
                          sigma2 = numeric(length = 1000))
  for (i in 1:1000) {
    X <- r_ou(n = 1, t = seq(0, T, len = n), x0 = 28, mu = mu, sigma = sqrt(sigma2), alpha = alpha)$data
    
    result <- est_OU(X, delta)
    
    fisher_info <- fisher_info_second_derivatives(result$alpha, result$mu, result$sigma2, X, delta,
                                                  d2_logL_dalpha2, d2_logL_dmu2, d2_logL_dsigma2,
                                                  d2_logL_dalphadmu, d2_logL_dalphadsigma, d2_logL_dmudsigma)
    
    # fisher_info_gradient <- fisher_info_gradient(result$alpha, result$mu, result$sigma2, X, delta,
    #                                              d_logL_dalpha, d_logL_dmu, d_logL_dsigma2)
    
    R <- chol(fisher_info)
    
    point <- sqrt(n) * (R %*% (c(result$alpha, result$mu, result$sigma2) - theta))
    data_1000[i,] <- point
  }
  
  # Combine data
  data_100$sample_size <- "100"
  data_1000$sample_size <- "1000"
  combined_data <- rbind(data_100, data_1000)
  
  # Compute correlations
  cor_100 <- cor(data_100[,1:3])
  cor_1000 <- cor(data_1000[,1:3])
  
  #print(cor_100)
  
  # Plot functions
  density_with_normal <- function(data, mapping, ...) {
    ggplot(data = data, mapping = mapping) +
      geom_density(aes(color = sample_size), alpha = 0.7) + 
      stat_function(fun = dnorm, args = list(mean = 0, sd = 1)) +
      theme_minimal() +
      xlim(c(-6, 6))
  }
  
  scatterplots <- function(data, mapping, ...) {
    ggplot(data = data, mapping = mapping) +
      geom_point(aes(color = sample_size), alpha = 0.7) +
      theme_minimal() +
      xlim(c(-6, 6)) +
      ylim(c(-6, 6))
  }
  
  upper_corr <- function(data, mapping, ...) {
    x_var <- as.character(mapping$x[2])
    y_var <- as.character(mapping$y[2])
    
    corr_100 <- round(cor_100[x_var, y_var], 4)
    corr_1000 <- round(cor_1000[x_var, y_var], 4)
    
    ggplot() + 
      annotate("text", x = 0.5, y = c(0.55, 0.5, 0.49, 0.45), 
               label = c("", paste(corr_100), paste(corr_1000), ""), 
               color = c("white", "red", "turquoise", "white"), size = 5, hjust = 0.5) +
      theme_minimal()
  }
  
  ggpairs(combined_data, columns = 1:3, 
          upper = list(continuous = wrap(upper_corr)),
          lower = list(continuous = scatterplots),
          diag = list(continuous = density_with_normal)
  )
}

# Generate and save the combined plot
#pdf("plots/3_pairwise_plot.pdf", width = 12, height = 12)
plot_combined <- pairwise_plot(mu = 20, sigma2 = 2, alpha = 3)
print(plot_combined)
#dev.off()

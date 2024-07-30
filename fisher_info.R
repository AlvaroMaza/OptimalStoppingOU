######### THEORETICAL FISHER INFORMATION #########
mu_fisher <- function (alpha, sigma2, delta) {
  return( 2*alpha/sigma2 * (((1-exp(-alpha*delta))^2)/(1-exp(-2*alpha*delta))) )
}

sigma2_fisher <- function(alpha, sigma2, delta) {
  return( 1/(2*(sigma2^2)) )
}


alpha_fisher <- function(alpha, sigma2, mu, x0, delta) {
  gamma <- exp(-alpha * delta)
  
  result <- - 1/(2*alpha^2) + 
    2 * delta^2*gamma^2/(1-gamma^2)^2 -
    2 * delta^2*gamma^2/(1-gamma^2) -
    4 * delta^2*gamma^2/(1-gamma^2)^2 -
    2 * alpha/sigma2 * delta^2*gamma^2/(1-gamma^2) *
    (sigma2/(2*alpha) * (1-gamma^2) + (x0 - mu)^2*gamma^2) +
    2 * delta*gamma^2/(1-gamma^2)
  
  return(result)
}
  


mu_fisher(3, 2, 0.01)
sigma2_fisher(3, 2, 0.01)
alpha_fisher(3, 2, 20, 14, 0.01)



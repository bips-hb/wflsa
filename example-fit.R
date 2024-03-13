library(wflsa)
library(flsa)
library(ggplot2)

# number of parameters. NOTE: is fixed here
p <- 100

# the average
mu <- as.vector(sapply(1:10*10, function(mu) rep(mu, 10)))
mu <- mu - mean(mu)

# generating the data
y <- rnorm(100, mean = mu)

# band_matrix function creates a square matrix with a band around the diagonal.
# Parameters:
#   - p: Size of the square matrix.
# Usage: band_matrix(p)
# Default band width is set to 1, but it can be adjusted as needed.

band_matrix <- function(p) {
  band_width <- 1  # Adjust the band width as needed
  
  # Create a matrix with a band around the diagonal
  my_matrix <- matrix(0, nrow = p, ncol = p)
  for (i in 1:p) {
    lower <- max(1, i - band_width)
    upper <- min(p, i + band_width)
    my_matrix[i, lower:upper] <- 1
    
  }
  diag(my_matrix) <- 0  # Set diagonal elements to 0 (optional, as they are already 1)
  return(my_matrix)
}

# creating the weight matrix commonly used for the 1-dimensional fused lasso signal approximator
W <- band_matrix(p)

fit <- wflsa::wflsa(y, W, lambda1 = .1, lambda2 = 2)
fit_flsa <- as.vector(flsa(y, lambda1 = .1, lambda2 = 2))

data <- dplyr::tibble(
  index = 1:p, 
  mu = mu,
  y = y, 
  beta_fit = fit$betas[[1]], 
  beta_fit_flsa = fit_flsa
)

ggplot2::ggplot(data) + 
  geom_point(mapping = aes(index, y), color = 'blue') + 
  geom_point(mapping = aes(index, beta_fit), color = 'red') + 
  geom_point(mapping = aes(index, beta_fit_flsa), color = 'purple') + 
  # geom_point(mapping = aes(index, beta_fit_flsa), color = 'green') + 
  ylab("value") + 
  xlab("variable") + 
  ggtitle("Raw data (blue) and fitted values (red)")


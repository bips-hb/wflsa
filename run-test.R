library(wfla)
library(genlasso)

m <- 9

lambda1 <- .4
lambda2 <- .01
rho <- 1

eta1 <- lambda1 / rho
eta2 <- lambda2 / rho


a <- 10

genlasso_implementation <-  function(){
  out <- genlasso::genlasso(y, diag(1, m), create_matrix_D(W, lambda1, lambda2, rho = rho, remove_zero_row = TRUE), minlam = 1)
  b   <- coef(out, lambda = 1)$beta 
  sapply(b, function(x) 
    if (abs(x) < 10^-7) {
      return(0)
    } else {
      return(x)
    })}

# draw random data for y

set.seed(2)
y <- rnorm(m)



W <- matrix(rep(.5, m*m), ncol = m)
genlasso_implementation()
wfla::genlasso_wrapper(y, W, m, m + (m)*(m-1)/2, eta1, eta2, a, rho = rho, max_iter = 1000, eps = 10^-10, truncate = 10^-5)

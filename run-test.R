library(wfla)
library(genlasso)

m <- 9

lambda1 <- .4
lambda2 <- .01
rho <- 1

eta1 <- lambda1 / rho
eta2 <- lambda2 / rho


a <- .4

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

set.seed(5)
y <- rnorm(m)


#eta1 = 0
#eta2 = 0
W <- matrix(rep(1, m*m), ncol = m)
W <- CVN::create_weight_matrix(type = "uniform-random", m = m)

(g = genlasso_implementation())
(w = wfla::genlasso_wrapper(y, W, m, eta1, eta2, a, rho = rho, 
                       max_iter = 1e5, eps = 10^-10, truncate = 10^-4))

microbenchmark::microbenchmark(genlasso_implementation(), 
                               wfla::genlasso_wrapper(y, W, m, eta1, eta2, a, rho = rho, 
                                                      max_iter = 1e5, eps = 10^-10, truncate = 10^-4), times = 4)

sum(abs(g - w))

sum(abs(w - y))

y / w


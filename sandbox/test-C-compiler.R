library(wfla)

m <- 9

lambda1 <- .4
lambda2 <- .8
rho <- 1

eta1 <- lambda1 / rho
eta2 <- lambda2 / rho

W <- matrix(rep(1, m*m), ncol = m)

a <- 10
rho_genlasso <- 1
maxiter_genlasso <- 1000
eps_genlasso <- 10^-10
truncate_genlasso <- 10^-5

# draw random data for y
y <- rnorm(m)

wfla::genlasso_wrapper(y, W, m, m + (m*(m-1)/2), eta1, eta2, a, rho = 1, max_iter = 1000, eps = 10^-10, truncate = 10^-5)


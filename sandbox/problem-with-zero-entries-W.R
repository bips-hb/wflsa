library(CVN)
library(CVNSim)

m <- 9

y <- rnorm(m)

lambda1 <- .4
lambda2 <- .8
rho <- 1

eta1 <- lambda1 / rho
eta2 <- lambda2 / rho

W <- CVNSim::create_weight_matrix("glasso")

D <- CVN::create_matrix_D(W, lambda1, lambda2,remove_zero_row = FALSE)

a <- CVN::matrix_A_inner_ADMM(m, D) + 1
rho_genlasso <- 1
maxiter_genlasso <- 1000
eps_genlasso <- 10^-10
truncate_genlasso <- 10^-5

CVN::genlasso_wrapper(y, W, m, nrow(D), eta1, eta2, a, rho = 1, max_iter = 1000, eps = 10^-10, truncate = 10^-5)


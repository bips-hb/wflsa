library(genlasso)
library(CVN)
library(matrixcalc)
library(microbenchmark)
library(glmnet)

# trial for new generalized LASSO estimation
m <- 4 # number of graphs

global_rho <- 1
lambda1 <- .2
lambda2 <- .7

n1 <- global_rho * lambda1
n2 <- global_rho * lambda2


W <- matrix(1, m, m)
e = 0
W <- matrix(c(0, 1, e, e,
              1, 0, 1, e,
              e, 1, 0, 1,
              e, e, 1, 0), ncol = 4 )
# W <- matrix(runif(m*m), ncol = m)
# W <- W %*% t(W) 
# W <- W / max(W) 
# W <- W 
# diag(W) <- 0
#W <- matrix(rbinom(m*m, 1, .5), ncol = m) * matrix(runif(m*m), ncol = m)
#W <- matrix(rbinom(m*m, 1, .5), ncol = m)
#W <- W %*% t(W) 
#diag(W) <- 0

D <- create_matrix_D(W, lambda1, lambda2, global_rho)
c <- nrow(D)

-t(D) %*% D

a = CVN::matrix_A_inner_ADMM(m, D)-.0001

is.positive.semi.definite(diag(a, m) - t(D) %*% D)


a = n1^2 + 3*n2^2

y <- rnorm(m, mean = 0, sd = 1)

# apply the generalized LASSO 
out <- genlasso::genlasso(y, diag(1, m), D, minlam = 1)
coef(out, lambda = 1)$beta

altZ(y, D, W, n1, n2, diagA = 20, rho = 1, max_iter = 1000)
  

CVN::aug_genlasso(y, W, m, c, lambda1, lambda2, global_rho, a, rho = 1, max_iter = 1000, eps = 10^-10)
  
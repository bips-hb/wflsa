library(genlasso)
library(CVN)
library(matrixcalc)
library(microbenchmark)
library(glmnet)
library(tictoc) 

source("R-implementation-genlasso.R")

# trial for new generalized LASSO estimation
m <- 4 # number of graphs

repetitions <- 20  

est_script   <- matrix(rep(0, repetitions*m), nrow = repetitions)
est_cpp      <- matrix(rep(0, repetitions*m), nrow = repetitions)
est_genlasso <- matrix(rep(0, repetitions*m), nrow = repetitions)

# Choose parameters and weight matrix
lambda1    <- .2
lambda2    <- 1
global_rho <- 1
rho        <- 1
eta1       <- lambda1 / global_rho
eta2       <- lambda2 / global_rho
truncate   <- 1e-5
eps        <- 1e-12

# Choose of W
uniform_random <- T
predefined <- F
disconnected <- F

W <- matrix(1, m, m) # standard fully-connected 

if (disconnected) { 
  W <- matrix(0, m, m)
}

if (uniform_random) { 
  W <- matrix(runif(m*m), ncol = m)
  W <- W %*% t(W)
  W <- W / max(W)
  W <- W
  diag(W) <- 0
}

if (predefined) { 
  e <- 0
  W <- matrix(c(0, 1, e, 1,
                1, 0, 1, e,
                e, 1, 0, 1,
                1, e, 1, 0), ncol = 4)
}

# create D and determine a
D <- create_matrix_D(W, lambda1, lambda2, rho = global_rho, remove_zero_row = FALSE)
a <- CVN::matrix_A_inner_ADMM(m, D) + 1

# Approaches ---------------
# R-implementation 
r_implementation <- function(){altZ(y, D, W, lambda1, lambda2, global_rho, diagA = a, max_iter = 1000, old = FALSE, eps = eps, truncate = truncate)}

genlasso_implementation <-  function(){
  out <- genlasso::genlasso(y, diag(1, m), create_matrix_D(W, lambda1, lambda2, rho = global_rho, remove_zero_row = TRUE), minlam = 1)
  b   <- coef(out, lambda = 1)$beta 
  sapply(b, function(x) 
    if (abs(x) < 10^-7) {
      return(0)
    } else {
      return(x)
  })}

cpp_implementation <- function() { 
  CVN::aug_genlasso(y, W, as.integer(m), nrow(D), eta1, eta2, a,
                   1, 1000, eps, truncate)
}

# perform tests
tic()
for (r in 1:repetitions) { 

  y <- rnorm(m, mean = 0, sd = 1)

  est_script[r,]   <- r_implementation()
  est_cpp[r,]      <- cpp_implementation()
  est_genlasso[r,] <- genlasso_implementation()
}
toc()

cat(sprintf("max. abs. difference ----------\n"))
cat(sprintf("\tscript vs. cpp:\t\t%f\n",    max(abs(est_script - est_cpp))))
cat(sprintf("\tgenlasso vs. script:\t%f\n", max(abs(est_genlasso - est_script))))
cat(sprintf("\tgenlasso vs. cpp:\t%f\n",    max(abs(est_genlasso - est_cpp))))

hist(est_genlasso - est_script, breaks = 20)

#print(summary(abs(est_genlasso - est_script)))

#f()
#as.vector(g())
#h()

#microbenchmark(f(), g(), h(), times = 4)


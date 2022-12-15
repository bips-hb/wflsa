library(microbenchmark)


# generate a weight matrix 
toy_weightMatrix <- function(m = 5) { 
  matrix(1, nrow = m, ncol = m) 
}


# The sandbox to test ideas and code
toy_sym_matrix <- function(p = 20, mean = 0, sd = 1) { 
  s <- matrix(rnorm(p*p, mean = mean, sd = sd), p)
  s[lower.tri(s)] = t(s)[lower.tri(s)]
  s
}

toy_list_sym_matrices <- function(m = 10, p = 20, mean = 0, sd = 1) { 
  lapply(1:m, function(i) {
    toy_sym_matrix(p = p, mean = mean, sd = sd)
  }) 
}

# generate a simple dataset X
toy_generateX <- function(n = 20, p = 30, mean = 0, sd = 1) { 
  matrix(rnorm(n*p, mean = mean, sd = sd), ncol = p) 
}

# generate list with square matrices


toy_generateRawDataset <- function(n = rep(20,10), p = 30, mean = 0, sd = 1) { 
  lapply(n, function(ni) { 
    toy_generateX(ni, p = p, mean = mean, sd = sd)
  }) 
}


seed(1)
data = toy_generateRawDataset(rep(50, 10), p = 100)
W = toy_weightMatrix(m = 10)
lambda1 = 1
lambda2 = 1
rho = 1
# Check correctness input -------------------------------
check_correctness_input(data, W, lambda1, lambda2, rho)

system.time(CVN(data, 
    W = W, 
    lambda1 = lambda1, 
    lambda2 = lambda2, 
    rho = rho,
    epsilon = 10^(-2),
    maxiter = 100, 
    n_cores = 1, 
    normalized = FALSE, 
    verbose = TRUE))





W <- matrix(0, 4, 4)
D <- create_matrix_D(W, 1, .3)
D
row_with_only_zeros <- apply(D, 1, function(row) all(row == 0))
D[!row_with_only_zeros, ]

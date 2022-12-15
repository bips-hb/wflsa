library(CVN)
library(CVNSim)

p <-  100
m <- 9
lambda1 <- .7
lambda2 <- .1
rho <- 1
normalized = T
maxiter = 100
n = 100

i = 1

rho = 1
eps = 1e-4
maxiter = 100 
truncate = 1e-5
rho_genlasso = 1
eps_genlasso = 1e-10
maxiter_genlasso = 100
truncate_genlasso = 1e-6
n_cores = min(length(lambda1)*length(lambda2), parallel::detectCores() - 1)
normalized = FALSE
warmstart = FALSE
use_genlasso_package = FALSE
minimal = FALSE
verbose = TRUE
use_new_updateZ = TRUE

eta1 <- lambda1 / rho
eta2 <- lambda2 / rho

set.seed(1)
starting_graph <- CVNSim::generate_graph(p = p, type = "random", probability = .5)
grid_of_graph <- CVNSim::create_grid_of_graphs(starting_graph = starting_graph, 
                                               n_edges_added_x = 2, 
                                               n_edges_removed_x = 1, 
                                               n_edges_added_y = 3, 
                                               n_edges_removed_y = 4)

# generate a single raw dataset for each graph in the grid
#set.seed(2)
data <- CVNSim::generate_raw_data_grid(n = n, grid_of_graph) 

W <- CVNSim::create_weight_matrix("grid")
D <- CVN::create_matrix_D(W, lambda1 = lambda1, lambda2 = lambda1, rho = rho)

a <- CVN::matrix_A_inner_ADMM(m, D) + 1

#cvn <- CVN(data, W, use_new_updateZ = F)

# Extract variables -------------------------------------
m <- length(data)       # total number of graphs  
p <- ncol(data[[1]])    # total number of variables
n_obs <- sapply(data, function(X) nrow(X))    # no. of observations per graph 


# Center or normalize data -------------------------------
data <- lapply(data, function(X) scale(X, scale = normalized))

# Compute the empirical covariance matrices --------------
Sigma <- lapply(1:m, function(i) cov(data[[i]])*(n_obs[i] - 1) / n_obs[i])

# initialize results list ------------
global_res <- list(
  Theta    = list(),          
  adj_matrices = list(),   
  Sigma    = Sigma,
  m        = m, 
  p        = p, 
  n_obs    = n_obs,
  data     = data,
  W        = W, 
  maxiter  = maxiter, 
  rho      = rho, 
  eps      = eps,
  truncate = truncate, 
  rho_genlasso      = rho_genlasso, 
  eps_genlasso      = eps_genlasso, 
  maxiter_genlasso  = maxiter_genlasso, 
  truncate_genlasso = truncate_genlasso, 
  n_lambda_values   = length(lambda1) * length(lambda2), 
  normalized = normalized,
  warmstart  = warmstart, 
  minimal = minimal, 
  use_genlasso_package  = use_genlasso_package
)

# data frame with the results for each unique (lambda1,lambda2) pair
res <- data.frame(expand.grid(lambda1 = lambda1, 
                              lambda2 = lambda2, 
                              converged = FALSE, 
                              value = NA, 
                              n_iterations = NA, 
                              aic = NA))
res$id <- 1:nrow(res)

# Create initial values for the ADMM algorithm ----------
Z <- rep(list(matrix(0, nrow = p, ncol = p)), m) # m (p x p)-dimensional zero matrices
Y <- rep(list(matrix(0, nrow = p, ncol = p)), m)

# if warmstart, the individual graphs are first estimated using the GLASSO.
if (warmstart) { 
  # We use the lambda1 value 
  Theta <- lapply(Sigma, function(S) { 
    est <- glasso::glasso(s = S, rho = res$lambda1[i])
    est$w
  })
}  else { 
  Theta <- lapply(Sigma, function(S) diag(1/diag(S)))
}

# Initialize variables for the algorithm -----------------
# Generate matrix D for the generalized LASSO 
D <- CVN::create_matrix_D(W, res$lambda1[i], res$lambda2[i], rho, 
                          remove_zero_row = use_genlasso_package)

a <- CVN::matrix_A_inner_ADMM(m, D) + 1

# Estimate the graphs -------------------------------------
eta1 <- res$lambda1[i] / rho 
eta2 <- res$lambda2[i] / rho 

# keep track whether the algorithm finished, either by 
# whether the stopping condition was met, or the maximum
# number of iterations was reached
converged  <- FALSE 
iter       <- 1 
difference <- 1 # stores the discrepancy 

# Create initial values for the ADMM algorithm ----------
Z <- rep(list(matrix(0, nrow = p, ncol = p)), m) # m (p x p)-dimensional zero matrices
Y <- rep(list(matrix(0, nrow = p, ncol = p)), m)

# if warmstart, the individual graphs are first estimated using the GLASSO.
if (warmstart) { 
  # We use the lambda1 value 
  Theta <- lapply(Sigma, function(S) { 
    est <- glasso::glasso(s = S, rho = lambda1)
    est$w
  })
}  else { 
  Theta <- Sigma
}

Theta <- Sigma

# Initialize a Temp variable for Theta
Temp  <- Theta
Theta <- Theta
Z     <- Z
Y     <- Theta

Z_new <- CVN::updateZ_wrapper(m, p, m*(m-1)/2 + m, 
                   as.matrix(Theta), as.matrix(Y), W, eta1, eta2, a, 
                   rho_genlasso, maxiter_genlasso, eps_genlasso, 
                   truncate_genlasso)  
Z_new <- as.list(Z_new[,1])

Z_old <- CVN::updateZ(m, p, m*(m-1)/2 + m, Theta, Y, W, 
                    eta1, eta2, a, 
                    rho_genlasso, maxiter_genlasso, eps_genlasso, 
                    truncate_genlasso, use_genlasso_package) 

sapply(1:9, function(i) sum(Z_new[[i]] - Z_old[[i]]))

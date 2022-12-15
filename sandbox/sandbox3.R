library(CVNSim)
library(tictoc)
library(microbenchmark)
library(profvis)
library(CVN)

p <- 10
m <- 9

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



starting_graph <- CVNSim::generate_graph(p, type = "random", probability = .5)
grid3 <- create_grid_of_graphs(starting_graph = starting_graph, 
                      n_edges_added_x = 2, 
                      n_edges_removed_x = 2, 
                      n_edges_added_y = 2,
                      n_edges_removed_y = 2,
                      verbose = TRUE)

data <- CVNSim::generate_raw_data_grid(100, grid3)

lambda1 = seq(.1,20, length.out = 5)
lambda2 = seq(.1,20, length.out = 5)

microbenchmark(f(), g())
tic()
cvn <- CVN::CVN(data = data, W, lambda1 = lambda1, lambda2 = lambda2, 
                eps = 10^-2, maxiter = 1000, 
                  verbose = TRUE, warmstart = T, use_previous_estimate = F, 
                  use_genlasso_package = T)
toc()

genl <- cvn$adj_matrices[[1]][[9]]

tic()
cvn <- CVN::CVN(data = data, W, lambda1 = lambda1, lambda2 = lambda2, 
                eps = 10^-3, maxiter = 1000, 
                verbose = T, warmstart = T, use_previous_estimate = F, 
                use_genlasso_package = F)
toc()

cpp <- cvn$adj_matrices[[1]][[9]]

g = grid3[[1]]$adj_matrix
cvn$adj_matrices[[2]][[1]]



cvn1 <- CVN::CVN(data = data, W, lambda1 = lambda1, lambda2 = lambda2, 
                 epsilon = 10^-3, maxiter = 1000, 
                 verbose = TRUE, warmstart = T, use_previous_version = TRUE)

tic()
profvis(cvn2 <- CVN::CVN(data = data, W, lambda1 = lambda1, lambda2 = lambda2, 
                 epsilon = 10^-2, maxiter = 1000, 
                 verbose = TRUE, warmstart = T, use_previous_version = FALSE))
toc()

grid3$`(1,1)`
cvn2$adj_matrices[[3]][1]

grid3$`(1,2)`
cvn2$adj_matrices[[3]][2]



cvn2$results

cvn <- CVN::CVN(data = data, W, lambda1 = c(.1), lambda2 = c(.1), 
                 epsilon = 10^-3, maxiter = 1000, 
                 verbose = TRUE, warmstart = FALSE) 




tic()
cvn2 <- CVN::CVN(data = data, W, lambda1 = lambda1, lambda2 = lambda2, verbose = T, epsilon = 10^-2, maxiter = 1000, use_previous_version = TRUE)
toc()

tic()
m1 <- microbenchmark(
  cvn1 <- CVN::CVN(data = data, W, lambda1 = lambda1, lambda2 = lambda2, epsilon = 10^-2, maxiter = 1000, verbose = TRUE, use_previous_version = TRUE), 
  times = 1
)

m2 <- microbenchmark(
  cvn2 <- CVN::CVN(data = data, W, lambda1 = lambda1, lambda2 = lambda2, epsilon = 10^-2, maxiter = 1000, verbose = FALSE), 
  times = 1
)

m3 <- microbenchmark(
  cvn3 <- CVN::CVN(data = data, W, lambda1 = lambda1, lambda2 = lambda2, epsilon = 10^-2, maxiter = 1000, verbose = TRUE, warmstart = TRUE), 
  times = 1
)

m2 <- microbenchmark(
  cvn1 <- CVN::CVN(data = data, W, lambda1 = lambda1, lambda2 = lambda2, epsilon = 10^-2, maxiter = 1000, verbose = FALSE, use_previous_version = T, warmstart = F), 
  times = 3
)
toc()

cvn2$value
cvn2$n_iterations

cvn1$value
cvn1$n_iterations


for (i in 1:4) { 
  for (j in 1:9) { 
    print(hamming_distance(cvn1$adj_matrices[[i]][[j]], cvn2$adj_matrices[[i]][[j]]))
  }
}

i = 1
j = 3
cvn1$adj_matrices[[i]][[j]]
cvn2$adj_matrices[[i]][[j]]

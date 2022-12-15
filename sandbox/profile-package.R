library(CVN)
library(CVNSim)
library(profvis)


p <-  10
m <- 9
lambda1 <- 1
lambda2 <- 1
rho <- 1
n = 100

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
set.seed(2)
data <- CVNSim::generate_raw_data_grid(n = n, grid_of_graph) 

W <- CVNSim::create_weight_matrix("full", m = m)
#D <- CVN::create_matrix_D(W, lambda1 = 1, lambda2 = 1, rho = 1)

#a <- CVN::matrix_A_inner_ADMM(m, D) + 1

#cvn <- CVN(data, W, use_new_updateZ = F)
#cvn_new <- CVN(data, W, use_new_updateZ = T)


#Z = updateZ_wrapper(m, p, p*(p-1)/2, as.matrix(cvn$Theta[[1]]), as.matrix(cvn$Sigma), W, eta1, eta2, a, 
#                rho = 1,  maxiter_genlasso = 1000, eps_genlasso = 1e-4, truncate_genlasso = 1e-5)


#library(microbenchmark)
#microbenchmark({CVN(data, W, use_new_updateZ = F)}, {CVN(data, W, use_new_updateZ = T)}, times = 1)


library(tictoc)
tic()

cvn_new <- CVN(data, W, n_cores = 1, maxiter = 1000, lambda1 = lambda1, lambda2 = lambda2, verbose = F)
toc()

 
tic()
cvn <- CVN(data, W, use_new_updateZ = F, n_cores = 2, lambda1 = lambda1, lambda2 = lambda2, verbose = T)
toc()

#as.list(Z[,1])


# 
# cvn$Theta[[1]][[1]] + cvn$Sigma[[1]]
# Z[[1]] - cvn$Theta[[1]][[1]]
# 
sapply(1:cvn$n_lambda_values, function(i) { 
   sapply(1:cvn$m, function(j) { 
     sum(cvn$Theta[[i]][[j]] - cvn_new$Theta[[i]][[j]])
   })
 })
# sum(cvn$adj_matrices[[1]][[1]] - cvn_new$adj_matrices[[1]][[1]])

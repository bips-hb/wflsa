library(CVN) 
library(CVNSim)
library(parallel)
library(profvis)

set.seed(1)

m = 2
n = 100
p = 5

graphs <- lapply(1:m, function(i) CVNSim::generate_graph(p = p))
data <- lapply(graphs, function(g) generate_data_given_adjaceny_matrix(n = n, g$adj_matrix))

W <- matrix(.5, m, m)

res = CVN::CVN(data, W = W, maxiter = 100, epsilon = 1e-3, lambda1 = 20, lambda2 = 2, rho = 1, n_cores = 8, verbose = TRUE)


profvis({
res = CVN::CVN(data, W = W, maxiter = 100, epsilon = 1e-3, lambda1 = 20, lambda2 = 2, rho = 1, n_cores = 8, verbose = TRUE)
})
Z = res$Z
Z[[1]]
graphs[[1]]$adj_matrix

mapply('+', res$Theta, res$Y, SIMPLIFY = FALSE)
rho = 1
n_obs <- rep(n, m)
Theta <- res$Theta
Z <- res$Z
Y <- res$Y
Sigma <- res$Sigma

t1 = res$Theta[[1]]

gl = glasso(cov(data[[1]]), 0.1, nobs = n)
t2 = gl$wi

t1
t2

a1 <- as.numeric(abs(t1) >= 1e-4)
a2  <- as.numeric(t2 != 0)

sum(abs(t1 - t2))

heatmap(t1)
hist(a1)


adj_matrices <- lapply(res$Theta, function(A) abs(A) >= 1e-5)

mapply(sum, adj_matrices)
sapply(graphs, function(g) sum(g$adj_matrix))

gl = glasso(cov(data[[1]]), 0.1, nobs = n)
gl$wi
View(gl$wi)
load("../MdS/Tests/alladjacs.RData")
load("../MdS/Tests/allsamples.RData")

data <- BigListof_samples[[1]][1:2]
true_adj_matrices <- BigListof_adjencies[[1]][1:2]

W <- matrix(1, 2, 2)

res <- CVN(data, W = W, maxiter = 1000, lambda1 = 10, lambda2 = 10, rho = 1, n_cores = 8, verbose = TRUE)

which( abs(res$Theta[[1]]) > 0.1)

which(true_adj_matrices[[1]] == 1)

W

D <- create_matrix_D(W, lambda1 = 1, lambda2 = 1)

y = rep(1,9)
r = genlasso(y, diag(1, 9), D = D, svd = TRUE)

beta = coef(r, lambda = 1)
beta$beta
library(R.utils)
isZero(beta$beta)
beta = beta$beta

beta[which(abs(beta) <= 10^-12)] <- 0

abs(beta$beta) <= 10^-12



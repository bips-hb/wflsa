library(microbenchmark)

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


toy_generateRawDataset <- function(m = 5, n = rep(20,m), p = 30, mean = 0, sd = 1) { 
  lapply(n, function(ni) { 
    toy_generateX(ni, p = p, mean = mean, sd = sd)
  }) 
}

# generate a weight matrix 
toy_weightMatrix <- function(m = 5) { 
  matrix(1, nrow = m, ncol = m) 
}



normalized = FALSE

X <- toy_generateRawDataset(mean = 2)

Y_old <- toy_generateRawDataset(mean = 0)
Theta_new <- toy_generateRawDataset()
Theta_old <- toy_generateRawDataset()
Z_new <- toy_generateRawDataset(mean = 0)
  
norm(Y_old[[2]])

lapply(X, function(X) X == 0)

# eigen decomposition ------
# m, Theta, Z, Y, Sigma, n_obs, rho = 1, n_cores = 1


m = 10
Sigma <- toy_list_sym_matrices()
Theta <- toy_list_sym_matrices()
Z <- toy_list_sym_matrices()
Y <- toy_list_sym_matrices()
rho = 1
n_cores = 1
n_obs = rep(20, m)

microbenchmark(updateTheta(m, Theta, Z, Y, Sigma, n_obs = n_obs), 
               updateTheta(m, Theta, Z, Y, Sigma, n_obs = n_obs, n_cores = 3))

mclapply(1:10, function(i) {
  
  # perform an eigendecomposition 
  eigen_decomposition <- eigen(Sigma[[i]] - (rho / n_obs[i]) * Z[[i]] + (rho / n_obs[i]) * Y[[i]])
  
  # obtain matrices Q and Lambda (Lambda is a diagonal matrix, so first only the diagonal is stored)
  Q <- eigen_decomposition$vectors 
  Lambda <- eigen_decomposition$values
  
  # Lambda is updated 
  Lambda <- n_obs[i]/(2 * rho) * (-Lambda + sqrt(Lambda^2 + 4*rho/n_obs[i]))
  
  # The new Theta matrix
  Q %*% diag(Lambda) %*% t(Q)}, mc.cores = n_cores # number of cores used
)

A = Sigma - (rho / n_obs[1]) * Z + (rho / n_obs[1]) * Y

eigen_decomposition <- eigen(A)
Q <- eigen_decomposition$vectors 
Lambda <- eigen_decomposition$values

Lambda <- n_obs[1]/(2 * rho) * (-Lambda + sqrt(Lambda^2 + 4*rho / n_obs[1]))
return(Q %*% diag(Lambda) %*% t(Q))






microbenchmark(sum(mcmapply(function(theta_new, theta_old) {norm(theta_new - theta_old)}, 
                            Theta_new, Theta_old)) / sum(mcmapply(norm, Theta_old)), 
               CVN::relative_difference_precision_matrices(Theta_new, Theta_old, n_cores = 1), 
               CVN::relative_difference_precision_matrices(Theta_new, Theta_old, n_cores = 8))

system.time(sum(mcmapply(function(theta_new, theta_old) {norm(theta_new - theta_old)}, 
       Theta_new, Theta_old)) / sum(mcmapply(norm, Theta_old)))

system.time(
  CVN::relative_difference_precision_matrices(Theta_new, Theta_old, n_cores = 1)
)

system.time(
CVN::relative_difference_precision_matrices(Theta_new, Theta_old, n_cores = 8)
)



mapply(norm, Theta_new - Theta_old)

mapply('-', Theta_new, Theta_old, SIMPLIFY = FALSE)
Theta_new


lapply()

y = mapply(function(y, theta, z) {y + theta - z}, Y_old, Theta_new, Z_new)

mclapply(1:5, function(i) Y_old[[i]] + Theta_new[[i]] - Z_new[[i]], mc.cores = 8)
b = mapply(function(y, theta, z) {y + theta + z}, Y_old, Theta_new, Z_new, SIMPLIFY = FALSE)

x = lapply(X, function(X) scale(X, center = TRUE, scale = normalized))

Sigma <- lapply(X, cov) 

lapply(X, mean)
lapply(x, mean)

lapply(X, mean)
lapply(X, scale)

r = scale(X[[1]])
mean(scale(X[[1]]))

mapply(function, ...)

X
mu <- lapply(X, mean)

X1 <- mapply('-', X, mu, SIMPLIFY = FALSE)
lapply(X1, mean)
mean(X)
lapply(X, function(x) mean(x))





library(genlasso)

set.seed(1)
n = 100
p = 10
X = matrix(rnorm(n*p),nrow=n)
y = 3*X[,1] + rnorm(n)
D = diag(1,p)
out = genlasso(y,X,D, minlam = 1)
coef(out, lambda=1)




m = 4
p = 4

Sigma <- toy_list_sym_matrices(m = m, p = p)
Theta <- toy_list_sym_matrices(m = m, p = p)
Z <- toy_list_sym_matrices(m = m, p = p)
Y <- toy_list_sym_matrices(m = m, p = p)

mapply[]
get_all_entries <- function(list_matrices, s, t) { 
   sapply(list_matrices, function(M) M[s,t])
}

set_entries_in_Z <- function(m, values, s, t) { 
  for (i in 1:m) { 
    Z[[i]][s,t] <<- values[i] 
    Z[[i]][t,s] <<- values[i] 
  }
}


get_all_entries(Z, 1,3)
set_entries_in_Z(4, c('a', 'b', 'c', 'd'), 1,3)

sapply(1:m, function(i) Sigma[[1]][i,i])


b <- mapply(rep, 1:m, 2, SIMPLIFY = FALSE)
lapply(b, function(i) get_all_entries(Sigma, i, i))
C = combn(1:m, 2, simplify = FALSE)
lapply(C, function(co) co[1] * co[2])

Theta 
Y
mapply()

mapply(function(theta, y) {diag( diag(theta) + diag(y) ) }, 
       Theta, Y, SIMPLIFY = FALSE)

diag(Theta[[1]])
Z <- rep(list(matrix(0, nrow = p, ncol = p)), m)

data = toy_generateRawDataset(m = 10, n = rep(40, m), p = 20)
W = toy_weightMatrix(m = 10)
lambda1 = 1
lambda2 = 1
rho = 1
  # Check correctness input -------------------------------
check_correctness_input(data, W, lambda1, lambda2, rho)

CVN(toy_generateRawDataset(m = 10, n = rep(40, m), p = 20), 
    W = toy_weightMatrix(m = 10), 
    lambda1 = 1, 
    lambda2 = 1, 
                rho = 1,
                epsilon = 10^(-5),
                maxiter = 1000, 
                n_cores = 1, 
                normalized = FALSE, 
                verbose = FALSE)

D <- create_matrix_D(W, 1, 1, 1)
qr(D)




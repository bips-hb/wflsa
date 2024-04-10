library(genlasso)
library(CVN)
library(wflsa)

p <- 9
y <- rnorm(p, mean = 1:p)
y <- y - mean(y)

W <- CVN::create_weight_matrix(type = 'uniform-random', m = p)

D <- CVN::create_matrix_D(W, 0, lambda2 = 1, rho = 1, remove_zero_row = T)

I <- diag(rep(1,p))

gf = genlasso(y, X = I, D = D)
gf$lambda[1]
gf$beta[,1]
gf$beta

f = wflsa(y = y, W = W, lambda2 = gf$lambda[32])

f$betas[[1]][[1]]


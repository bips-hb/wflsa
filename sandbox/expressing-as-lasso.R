library(CVN)

p <- 3
W <- CVN::create_weight_matrix(type = 'uniform-random', m = p)

# Creating matrix A such that beta = A %*% theta


ones_matrix <- function(p, k) {
  # Create an empty matrix filled with zeros
  result <- matrix(0, nrow = p, ncol = p)
  result[k, ] <- 1
  result[, k] <- -1
  result[k,k] <- 0 
  return(result)
}

create_row_matrix_A <- function(k, W) {
  p <- nrow(W)
  Ak <- ones_matrix(p, k)
  W1 <- 1/W
  W1[which(is.infinite(W1))] <- 0 
  W1 <- Ak * W1 / (p*(p-1)/2) #1/(p*(p-1)/2 - 1)
  diag(W1) <- 0
  return(t(vec(t(W1))))
}



create_row_matrix_A(1, W)
create_row_matrix_A(2, W)
create_row_matrix_A(3, W)

create_matrix_A <- function(W) {
  p <- nrow(W)
  
  res <- lapply(1:p, function(k) create_row_matrix_A(k, W))
  
  do.call(rbind, res)
}

A <- create_matrix_A(W)

# Identify columns with all zeros
# nonzero_cols <- apply(A, 2, function(x) any(x != 0))

# Subset the matrix to remove columns with all zeros
# A <- A[, nonzero_cols]

rowSums(A)

y <- rnorm(1:p)

y <- y - mean(y)

fit <- glmnet::glmnet(A, y, intercept = FALSE, offset = FALSE)
fit$beta
theta <- coef(fit, s = .1)[-1]
beta0 <- as.vector(A %*% theta)

theta <- predict(fit, s = .1, type = "coefficients")[-1]
beta0 <- as.vector(A %*% theta)

f <- wflsa::wflsa(y, W, lambda2 = .1, rho = 1, truncate = 1e-8)
sum(abs(f$betas[[1]][[1]] - beta0))
#theta


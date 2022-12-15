#' Create matrix \eqn{D} to be used for the Generalized LASSO
#'
#' Generates a matrix \eqn{D} to be used for the generalized LASSO. 
#' We solve a generalized LASSO problem for each edge \eqn{(s,t)} 
#' for each update step for \eqn{Z}. 
#' 
#' @param W The \eqn{(m \times m)}-dimensional upper-triangular 
#'          weight matrix \eqn{W}
#' @param lambda1 The \eqn{\lambda_1} LASSO penalty term 
#' @param lambda2 The \eqn{\lambda_2} global smoothing parameter 
#' @param rho The \eqn{\rho} ADMM's penalty parameter (Default: \code{1})
#' @param remove_zero_row If \code{TRUE}, rows with zeros are removed. 
#'                        (Default: \code{TRUE})
#' 
#' @return A \eqn{((m \cdot (m+1)/2) \times m)}-dimensional matrix 
#'         
#' @references Tibshirani, R. J., & Taylor, J. (2011). 
#'             The solution path of the generalized lasso. 
#'             Annals of Statistics, 39(3), 1335–1371. 
#'             https://doi.org/10.1214/11-AOS878
#' @examples 
#' m <- 4 # number of graphs
#' W <- matrix(1, nrow = m, ncol = m) 
#' 
#' # penalty terms:
#' lambda1 <- .2
#' lambda2 <- .4
#' rho <- 1
#' 
#' CVN::create_matrix_D(W, lambda1, lambda2, rho) 
#' #      [,1] [,2] [,3] [,4]
#' # [1,]  0.2  0.0  0.0  0.0
#' # [2,]  0.0  0.2  0.0  0.0
#' # [3,]  0.0  0.0  0.2  0.0
#' # [4,]  0.0  0.0  0.0  0.2
#' # [5,]  0.4 -0.4  0.0  0.0
#' # [6,]  0.4  0.0 -0.4  0.0
#' # [7,]  0.4  0.0  0.0 -0.4
#' # [8,]  0.0  0.4 -0.4  0.0
#' # [9,]  0.0  0.4  0.0 -0.4
#' # [10,]  0.0  0.0  0.4 -0.4
#' @export
create_matrix_D <- function(W, lambda1, lambda2, rho = 1, remove_zero_row = TRUE) { 
  
  eta1 <- lambda1 / rho 
  eta2 <- lambda2 / rho 
    
  # number of matrices that need to be estimated
  m <- nrow(W)
    
  # get all unique pairs, i.e., (1,2), (1,3) ... (m - 2, m), (m - 1, m)
  all_unique_pairs <- combn(1:m, 2) 
  
  # number of unique pairs 
  n_unique_pairs <- m*(m-1)/2 
  
  # matrix that will contain all pairs
  C <- matrix(rep(0, n_unique_pairs*m), ncol = m)
  
  # fill in each row of C independently
  sapply(1:n_unique_pairs, function(i) { 
      indices <- all_unique_pairs[,i] 
      combination <- rep(0, m) 
      combination[indices] <- c(1, -1) * W[indices[1], indices[2]]
      C[i, ] <<- combination
    })
  
  # identity matrix 
  I <- diag(m)
  
  # combinating the identity matrix with C
  D <- rbind(eta1 * I, eta2 * C)
  
  # remove any unnecessary rows from matrix D, i.e., 
  # rows with only zeros
  if (remove_zero_row) { 
    row_with_only_zeros <- apply(D, 1, function(row) all(row == 0))
    D <- D[!row_with_only_zeros, ]
  }
  
  return(D)
}
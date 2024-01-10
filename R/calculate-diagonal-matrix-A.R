#' Determine the Diagonal of Matrix \eqn{A} 
#' 
#' The wFLSA algorithm requires the choice of a matrix \eqn{A} such that 
#' \eqn{A - D'D} is positive semidefinite.  
#' We choose the matrix \eqn{A} to be diagonal with a fixed value \eqn{a}. 
#' This function determines the smallest value of \eqn{a} such that 
#' \eqn{A - D'D} is indeed positive semidefinite. We do this be determining 
#' the largest eigenvalue
#' 
#' @param W Weight matrix \eqn{W}
#' @param eta1,eta2 The values \eqn{\eta_1 = \lambda_1 / \rho} and 
#'                   \eqn{\eta_2 = \lambda_2 / \rho}
#' 
#' @return Value of \eqn{a}
#' @references 
#' Zhu, Y. (2017). An Augmented ADMM Algorithm With Application to the 
#' Generalized Lasso Problem. Journal of Computational and Graphical Statistics, 
#' 26(1), 195–204. https://doi.org/10.1080/10618600.2015.1114491
#' @export
calculate_diagonal_matrix_A <- function(W, eta1, eta2) { 
  
  # determine number of variables
  p <- ncol(W) 
  
  # determine D'D 
  DtD <- -eta2^2 * (W*W) + eta1^2 * diag(p) + eta2^2*diag(rowSums(W*W))
  
  # determine the eigenvalues
  eigenvectors <- eigen(DtD)
  
  # return the value of a
  return(max(eigenvectors$values))
}
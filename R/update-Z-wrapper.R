#' Wrapper for the \eqn{Z}-update Step for the ADMM
#' 
#' A wrapper for the \code{C} function that returns the 
#' updated value of \eqn{Z} for the ADMM given the previously 
#' updated values of \eqn{\Theta} and \eqn{Y}
#' 
#' @param m Number of graphs
#' @param p Number of variables
#' @param nrow_D Number of rows of the \eqn{D}-matrix
#' @param Theta A list with matrices with the current values of \eqn{\Theta}
#' @param Y A list with matrices with the current values of \eqn{Y} 
#' @param W Weight matrix 
#' @param eta1 
#' @param rho_genlasso The \eqn{\rho} penalty parameter for the ADMM algorithm 
#' @param eps_genlasso If the relative difference between two update steps is 
#'                smaller than \eqn{\epsilon}, the algorithm stops
#' @param maxiter_genlasso Maximum number of iterations for solving 
#'                the generalized LASSO problem 
#' @param truncate_genlasso All values of the final \eqn{\hat{\beta}} below 
#'                 \code{truncate_genlasso} will be set to \code{0}. 
#' 
#' @return A list with matrices with the new values of \eqn{Z}
#'
#' @seealso \code{\link{create_matrix_D}} 
#' 
#' @export
updateZ_wrapper <- function(m, p, nrow_D, 
                            Theta, Y, W, eta1, eta2, a, 
                            rho_genlasso, maxiter_genlasso, eps_genlasso, 
                            truncate_genlasso) { 
  
  updateZRcpp(m, p, nrow_D, 
               Theta, Y, W, eta1, eta2, a, 
               rho_genlasso, maxiter_genlasso, eps_genlasso, 
               truncate_genlasso)  
  
}
  
#' The Weighted Fused Lasso Signal Approximator (wFLSA) Algorithm
#'
#' Solves the weighted Fused LASSO Signal Approximator optimization problem 
#' using an ADMM-based approach. The problem is formulated as follows: 
#' \deqn{
#' \hat{\beta} = \argmin \frac{1}{2} || y - \beta ||_2^2 + \lambda_1 ||\beta ||_1 + \lambda_2 \sum_{i < j} w_{ij} | \beta_i - \beta_j |  
#' }
#' where:
#' \itemize{
#'  \item \eqn{y} is the response with mean \eqn{0}.
#'  \item \eqn{\beta} is the vector of coefficients to be estimated.
#'  \item \eqn{|| \cdot ||_1} and \eqn{|| \cdot ||_2} are the \eqn{L_1}- and \eqn{L_2}-norms, respectively.
#'  \item \eqn{\lambda_1 > 0} is the regularization parameter controlling the strength of the sparsity penalty. 
#'  \item \eqn{\lambda_2 > 0} is the regularization parameter controlling the smoothness. 
#'  \item \eqn{w_{ij} \in [0,1]} is the weight between the \eqn{i}-th and \eqn{j}-th coefficient.
#' }
#' This implementation solves the entire path for a fixed value of \eqn{\lambda_2}. 
#' 
#' @param y Vector of length \eqn{p} representing the response variable (assumed to be centered).
#' @param W Weight matrix of dimensions \eqn{p \times p}.
#' @param lambda2 Vector of positive regularization parameters for smoothness penalty (Default: \code{(.1, .2)})
#' @param rho ADMM's parameter (Default: \code{1}).
#' @param max_iter Maximum number of iterations (Default: \code{1e5}).
#' @param eps Stopping criterion. If differences are smaller than \code{eps}, 
#'            the algorithm halts (Default: \code{1e-10}).
#' @param truncate Values below \code{truncate} are set to \code{0} (Default: \code{1e-4}).
#' @param offset Logical indicating whether the data is centered (Default: \code{TRUE}).
#'
#' @return A list containing:
#' \itemize{
#'  \item \code{betas}: A list with estimates. Each entry are the \eqn{\beta}-estimates
#'        for a specific value of \eqn{\lambda_2}. 
#'  \item all input variables.
#' }
#' 
#' @note
#' \strong{Important Note:} The algorithm assumes \eqn{y} to be centered, i.e., its mean is 0. 
#' 
#' @seealso \code{\link{genlassoRcpp}}
#' @examples
#' # Example usage of the wflsa function
#' y <- c(1, 2, 3)
#' W <- matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1), ncol = 3)
#' lambda2 <- c(0.1, 0.2)
#' result <- wflsa::wflsa(y, W, lambda2)
#' wflsa::get_estimate_lambda1(result, lambda1 = .5)
#' @export
wflsa <- function(y, W, lambda2 = c(.1, .2), rho = 1, 
                  max_iter = 1e5, eps = 1e-10, truncate = 1e-4, offset = TRUE) {
  
  # number of variables
  p <- length(y)
  
  # Check if dimensions of matrix W match the length of vector y
  if (ncol(W) != p || nrow(W) != p) {
    stop("Error: Dimensions of matrix W must be equal to the length of vector y.")
  }
  
  # Check if lambda2 is a non-empty vector with all positive values
  if (!is.vector(lambda2) || length(lambda2) == 0 || any(lambda2 <= 0)) {
    stop("Error: lambda2 must be a non-empty vector with all positive values.")
  }
  
  # Check if rho is positive
  if (rho <= 0) {
    stop("Error: rho must be a positive value.")
  }
  
  # Check if max_iter is a positive integer
  if (max_iter <= 0) {
    stop("Error: max_iter must be a positive integer.")
  }
  
  # Check if eps is positive
  if (eps <= 0) {
    stop("Error: eps must be a positive value.")
  }
  
  # Check if truncate is positive
  if (truncate <= 0) {
    stop("Error: truncate must be a positive value.")
  }
  
  # Check if offset is logical
  if (!is.logical(offset)) {
    stop("Error: offset must be a logical value.")
  }
  
  if (offset) {
    y <- y - mean(y)
  }
  
  # Determine the eta2 values
  eta2 <- lambda2 / rho 
  
  # Determine the value of a such that aI - D'D is positive definite (see paper)
  a <- wflsa::calculate_diagonal_matrix_A(W, 0, max(eta2)) + 1
  
  # Calculate the regression coefficients for the different lambda2
  # Note that lambda1 = 0, since we can determine the values for lambda1 later 
  # through soft thresholding
  betas <- lapply(eta2, function(eta2) {
    genlassoRcpp(y, W, p, 0, eta2, a, rho, max_iter, eps, truncate)
  })
  
  # Determine the values of lambda1 for which at least one coefficient switches to non-zero
  lambda1 <- lapply(betas, function(beta) c(0, wflsa::get_unique_values(beta, digits = 7)))
  
  # Determine the number of breakpoints
  n_lambda1 <- sapply(lambda1, function(breaks) length(breaks))
  
  # Determine the beta estimates for each breakpoint for each value of lambda2
  betas <- lapply(1:length(betas), function(i) {
    beta <- betas[[i]]
    
    # apply the soft thresholding to determine the estimate for a given lambda1 value
    lapply(lambda1[[i]], function(breakpoint) {
      wflsa::soft_threshold(beta, breakpoint + 1e-7)
    })
  })
  
  result <- list(
    betas = betas,
    lambda1 = lambda1, 
    n_lambda1 = n_lambda1, 
    y = y,
    W = W,
    lambda2 = lambda2,
    rho = rho,
    max_iter = max_iter,
    eps = eps,
    truncate = truncate,
    offset = offset
  )
  
  class(result) <- c('wflsa.fit', class(result))
  
  return(result)
}

#' Print Function for the Weighted Fused LASSO Signal Approximator (wFLSA) fit object
#'
#' This function prints information about the fitted wFLSA object, including the 
#' number of variables, the number of lambda pairs, and the estimated beta coefficients.
#'
#' @param fit An object of class wflsa.fit.
#' @param ... Additional arguments to be passed to the print function.
#'
#' @export
print.wflsa.fit <- function(fit, ...) { 
  
  p <- length(fit$y) 
  
  cat(sprintf("Fit for the Weighted Fused LASSO Signal Approximator\n\n"))
  
  # Print the number of variables (p) and the number of lambda pairs
  cat(sprintf("Number of variables (p)  : %d\n", p))
  cat(sprintf("Number of lambda2 values : %d\n", length(fit$lambda2)))
}


#' Get Estimates for Beta Coefficients for a Given Lambda1 Value
#'
#' This function takes the resulting fit from the wflsa function and returns the
#' estimates for the beta coefficients for all lambda2 values that were considered
#' before.
#'
#' @param fit The result object obtained from wflsa function
#' @param lambda1 The lambda1 value for which beta coefficients are needed
#' @return A list containing:
#'   \item{betas}{Matrix of beta coefficients for each lambda2 value. Each row corresponds to a single lambda2 value}
#'   \item{lambda2}{Vector of lambda2 values considered in the fit.}
#' @export
#' @examples
#' # Example usage:
#' # fit <- wflsa(data, lambda1_values, lambda2_values)
#' # estimate <- get_estimate_lambda1(fit, lambda1_value)
#' # Now you can access estimated beta coefficients for different lambda2 values
get_estimate_lambda1 <- function(fit, lambda1) {
  
  betas <- lapply(fit$betas, function(est) wflsa::soft_threshold(est[[1]], lambda1))
  
  list(
    betas = t(sapply(betas, function(x) x)), # turn into a matrix
    lambda2 = fit$lambda2
  )
}


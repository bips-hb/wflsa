#' The Weighted Fused Lasso Signal Approximator (wFLSA) Algorithm
#'
#' Solves the weighted Fused LASSO Signal Approximator optimization problem
#' using an ADMM-based approach. The problem is formulated as follows:
#' \deqn{
#' \hat{\beta} = \operatorname{argmin} \frac{1}{2} || y - \beta ||_2^2 + \lambda_1 ||\beta ||_1 + \lambda_2 \sum_{i < j} w_{ij} | \beta_i - \beta_j |
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
#'
#' @param y Vector of length \eqn{p} representing the response variable (assumed to be centered).
#' @param W Weight matrix of dimensions \eqn{p \times p}.
#' @param lambda1 Vector of positive regularization parameters for \eqn{L_1} penalty.
#' @param lambda2 Vector of positive regularization parameters for smoothness penalty.
#' @param rho ADMM's parameter (Default: `1`).
#' @param max_iter Maximum number of iterations (Default: `1e5`).
#' @param eps Stopping criterion. If differences are smaller than `eps`,
#'            the algorithm halts (Default: `1e-10`).
#' @param truncate Values below `truncate` are set to `0` (Default: `1e-4`).
#' @param offset Logical indicating whether to include an intercept term (Default: `TRUE`).
#'
#' @return A list containing:
#' \itemize{
#'  \item `betas`: Estimated vector \eqn{\hat{\beta}} from the Weighted Fused LASSO.
#'  \item `tuning_parameters`: Data frame with tuning parameters. The column `df`
#'        contains the number of non-zero coefficients for the different lambda-values
#'  \item all input variables.
#' }
#'
#' @note
#' **Important Note:** The algorithm assumes \eqn{y} to be centered, i.e., its mean is 0.
#'
#' @seealso [genlassoRcpp()]
#' @examples
#' # Example usage of the wflsa function
#' y <- c(1, 2, 3)
#' W <- matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1), ncol = 3)
#' lambda1 <- c(0.1, 0.2)
#' lambda2 <- c(0.1, 0.2)
#' result <- wflsa(y, W, lambda1, lambda2)
#'
#' @export
wflsa <- function(y, W, lambda1 = c(0.1), lambda2 = c(0.1), rho = 1,
                  max_iter = 1e5, eps = 1e-10, truncate = 1e-4, offset = TRUE) {

  # number of variables
  p <- length(y)

  # Check if dimensions of matrix W match the length of vector y
  if (ncol(W) != p || nrow(W) != p) {
    stop("Error: Dimensions of matrix W must be equal to the length of vector y.")
  }

  # Check if lambda1 is a non-empty vector with all positive values
  if (!is.vector(lambda1) || length(lambda1) == 0 || any(lambda1 <= 0)) {
    stop("Error: lambda1 must be a non-empty vector with all positive values.")
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

  # Get all possible combinations and determine the eta1 and eta2 values
  tuning_parameters <- data.frame(expand.grid(lambda1 = lambda1, lambda2 = lambda2))
  tuning_parameters$eta1 <- tuning_parameters$lambda1 / rho
  tuning_parameters$eta2 <- tuning_parameters$lambda2 / rho

  # Determine the value of a such that aI - D'D is positive definite (see paper)
  a <- wflsa::calculate_diagonal_matrix_A(W, max(tuning_parameters$eta1), max(tuning_parameters$eta2)) + 10

  # Calculate the regression coefficients for all combinations of lambda1 and lambda2
  betas <- lapply(1:nrow(tuning_parameters), function(i) {
    eta1 <- tuning_parameters$eta1[i]
    eta2 <- tuning_parameters$eta2[i]
    genlassoRcpp(y, W, p, eta1, eta2, a, rho, max_iter, eps, truncate)
  })

  # determine number of covariates not equal to zero
  df <- sapply(betas, function(beta) sum(beta != 0))

  result <- list(
    betas = betas,
    tuning_parameters = tuning_parameters,
    df = df,
    y = y,
    W = W,
    lambda1 = lambda1,
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
#' @param x An object of class wflsa.fit.
#' @param ... Additional arguments to be passed to the print function.
#'
#' @export
print.wflsa.fit <- function(x, ...) {

  # Number of lambda pairs:
  n_lambda_pairs <- length(x$lambda1) * length(x$lambda2)
  p <- length(x$y)

  cat(sprintf("Weighted Fused LASSO Signal Approximator\n\n"))

  # Print the number of variables (p) and the number of lambda pairs
  cat(sprintf("Number of variables (p) : %d\n", p))
  cat(sprintf("Number of lambda pairs  : %d\n\n", n_lambda_pairs))

  # Print the header for the estimated beta coefficients
  cat("Estimated beta coefficients \n")

  # Loop over the lambda pairs (n_lambda_pairs) - stopping after the first 5 pairs
  for (i in 1:min(n_lambda_pairs, 5)) {

    # Print the lambda pair values
    cat(sprintf("(%.2f, %.2f):\t", x$tuning_parameters[i, 'lambda1'],
                x$tuning_parameters[i, 'lambda2']))

    # Loop over the entries in each beta vector (p) - stopping after the first 8 entries
    for (j in 1:min(p, 8)) {
      # Print each beta value with two digits after the comma, separated by tabs
      cat(sprintf("%.2f\t", x$betas[[i]][j]))
    }

    # If there are more than 8 entries, print dots (...) to indicate additional entries
    if (p > 8) {
      cat("...")
    }

    # Move to the next line after printing each beta vector
    cat("\n")
  }

  # If there are more than 5 lambda pairs, print dots (...) to indicate additional pairs
  if (n_lambda_pairs > 5) {
    cat("...\n")
  }

  # Print the number of non-zero coefficients
  cat('\nNumber of non-zero coefficients:\n')
  cat(x$df)
  cat('\n')
}

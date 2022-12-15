#' Wrapper for \code{genlassoRcpp}
#' 
#' See for details \code{\link{genlassoRcpp}}
#' 
#' @seealso \code{\link{genlassoRcpp}}
#' @export
genlasso_wrapper <- function(y, W, m, c, eta1, eta2, a, rho, max_iter, eps, truncate) { 
  genlassoRcpp(y, W, m, c, eta1, eta2, a, rho, max_iter, eps, truncate) 
}

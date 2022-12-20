#include <Rcpp.h>
using namespace Rcpp;

//' Solving Generalized LASSO with fixed \eqn{\lambda = 1}
//' 
//' Solves efficiently the generalized LASSO problem of the form 
//' \deqn{
//'   \hat{\beta} = \text{argmin } \frac{1}{2} || y - \beta ||_2^2 + ||D\beta||_1 
//' }
//' where \eqn{\beta} and \eqn{y} are \eqn{m}-dimensional vectors and 
//' \eqn{D} is a \eqn{(c \times m)}-matrix where \eqn{c \geq m}. 
//' We solve this optimization problem using an adaption of the ADMM
//' algorithm presented in Zhu (2017). 
//' 
//' @param y The \eqn{y} vector of length \eqn{m}
//' @param W The weight matrix \eqn{W} of dimensions \eqn{m x m}
//' @param m The number of graphs 
//' @param c Number of rows of matrix \eqn{D}, which is equal to 
//'          \eqn{c = m + (m(m-1))/2}   
//' @param eta1 Equals \eqn{\lambda_1 / rho} 
//' @param eta2 Equals \eqn{\lambda_2 / rho} 
//' @param a Value added to the diagonal of \eqn{-D'D} so that
//'          the matrix is positive definite, see 
//'          \code{\link{matrix_A_inner_ADMM}}
//' @param rho The ADMM's parameter
//' @param max_iter Maximum number of iterations
//' @param eps Stopping criterion. If differences 
//'            are smaller than \eqn{\epsilon}, algorithm
//'            is halted
//' @param truncate Values below \code{truncate} are 
//'                 set to \code{0}
//'
//' @return The estimated vector \eqn{\hat{\beta}}
//'
//' @references 
//' Zhu, Y. (2017). An Augmented ADMM Algorithm With Application to the 
//' Generalized Lasso Problem. Journal of Computational and Graphical Statistics, 
//' 26(1), 195–204. https://doi.org/10.1080/10618600.2015.1114491
//' 
//' @seealso \code{\link{genlasso_wrapper}}
// [[Rcpp::export]]
NumericVector genlassoRcpp(const NumericVector y, 
                                const NumericMatrix& W, 
                                const int m, 
                                const int c, 
                                const double eta1, 
                                const double eta2, 
                                double a, 
                                const double rho, 
                                const int max_iter,
                                const double eps, 
                                const double truncate) { 
  /* some frequently used constants */
  a = rho*a ; 
  const double C = 1 / (1 + a) ; 
  
  /* initialize vectors for beta-update step in the ADMM */
  double beta_new[m] ; 
  double beta_old[m] ; 
  double delta[m] ; 
  
  for (int i = 0; i < m; i ++) { 
    *(beta_new + i) = 0;
    *(beta_old + i) = 0;
    *(delta + i) = 0;
  }
  //memset(beta_new, 0, m);
 // memset(beta_old, 0, m);
  //memset(delta, 0, m);
  
  // Eigen::VectorXd beta_new(m);
  // Eigen::VectorXd beta_old(m);
  // Eigen::VectorXd delta(m);
  
  /* initialize vectors for alpha-update step in the ADMM */
  double alpha_new[c] ; 
  double alpha_old1[c] ; 
  double alpha_old2[c] ; 
  double alpha[c] ; 
  
  for (int i = 0; i < c; i ++) { 
    *(alpha_new + i) = 0 ; 
    *(alpha_old1 + i) = 0 ; 
    *(alpha_old2 + i) = 0 ; 
    *(alpha + i) = 0 ; 
  }
  // 
  // memset(alpha_new, 0, c);
  // memset(alpha_old1, 0, c);
  // memset(alpha_old2, 0, c);
  // memset(alpha, 0, c);
  
  // Eigen::VectorXd alpha_new(c);
  // Eigen::VectorXd alpha_old1(c);
  // Eigen::VectorXd alpha_old2(c);
  // Eigen::VectorXd alpha(c);
  
  /* Indices used for the alpha-update step */
  int steps[m-1] ; 
  
  //std::vector<size_t> steps(m - 1) ; 
  //int steps[m-1] ; 
  *steps = 0 ; 
  
  for (int i = 1; i < (m-1); i ++) { 
    *(steps + i) = m - i + *(steps + i - 1) ; 
  }
  
  //return(NumericVector(beta_new ,beta_new + sizeof(beta_new) / sizeof(*beta_new)));
  
  int iter = 0 ;    // number of iterations 
  double diff = 0 ; // absolute difference between beta^(k+1) and beta^k
  
  /* loop until either the max. no. of iterations are reached 
   * or the difference (diff) is smaller then eps
   */
  while (iter < max_iter) { 
    
    /* ------- beta-update step ---------*/
    //alpha = 2*alpha_old1 - alpha_old2 ; 
    for (int i = 0; i < c; i++) {
      //alpha[i] = 2*alpha_old1[i] - alpha_old2[i] ; 
      *(alpha + i) = 2*(*(alpha_old1 + i)) - *(alpha_old2 + i) ;
      //Rcout << "i: " << i << "  alpha[i] = " << alpha[i] << "  alpha_old1[i] = " << alpha_old1[i] << "  alpha_old2[i] = " << alpha_old2[i] << std::endl ; 
    }
    
    // go over all possible pairs (i,j), same as D^T %*% (2 alpha^(k) - alpha^(k-1)) 
    //delta = eta1 * alpha ; 
    //double[] delta (m) ;
    
    for (int i = 0; i < m; i++) { 
      *(delta + i) = eta1 * (*(alpha + i)) ;  
      
      for (int j = i+1; j < m; j++) { 
        *(delta + i) += eta2* W(i,j)*(*(alpha + m + steps[i] + (j - i) - 1)) ; 
      }
      
      for (int j = 0; j < i; j++) { 
        *(delta + i) -= eta2*W(i,j)*(*(alpha + m + steps[j] - (j - i) - 1)) ; 
      }
    }
    
    // update beta with the computed delta and determine difference
    //double[] beta_new (m) ;
    // double[] beta_new (m) ;
    //beta_new = C*(a*beta_old + y - delta) ; 
    //diff = sum(abs(beta_new - beta_old)) ; 
    diff = 0 ;
    for (int i = 0; i < m; i ++) {
      *(beta_new + i) = C*(a*(*(beta_old + i)) + y[i] - *(delta + i)) ; 
      //*(beta_new + i) = C*(a*beta_old[i] + y[i] - delta[i]) ;
      diff += abs(beta_new[i] - beta_old[i]) ;
    }
    
    // determine whether converged or not
    if (diff < eps) { 
      /* Turn to zero when really close */
      for (int i = 0; i < m; i ++) { 
        if (fabs(beta_new[i]) < truncate) { 
          *(beta_new + i) = 0 ;  
        } 
      }
      
      //Eigen::VectorXd res(m);
      //std::copy(beta_new.begin(), beta_new.end(), res.begin()) ; 
      for (int i = 0; i < m; i ++) { 
        Rcout << *(beta_new + i) << " " ; 
      }
      Rcout << std::endl; 
      
      Rcout << "CONVERGED" << std::endl ; 
      
      return(NumericVector(beta_new, beta_new + sizeof(beta_new) / sizeof(*beta_new)));
    }
    
    /* --------- alpha update step ----------- */
    //alpha_new = Rcpp::clone(alpha_old1) + rho * eta1 * Rcpp::clone(beta_new) ; 
    // double[] alpha_new (c) ; // TODO change back
    //alpha_new = Rcpp::clone(alpha_old1) ; 
    for (int i = 0; i < m; i++) {
      *(alpha_new + i) = *(alpha_old1 + i) + rho * eta1 * *(beta_new + i) ;
    }
    
    int k = m; 
    // go over all unique pairs (i,j)
    for (int i = 0; i < m-1; i++) { 
      for (int j = i+1; j < m; j++) { 
        *(alpha_new + k) = *(alpha_old1 + k) + rho * eta2 * W(i,j) * (*(beta_new + i) - *(beta_new + j)) ; 
        //Rprintf("k = %d\t(i,j) = (%d,%d)\tW[%d,%d] = %g\t%g --> %g\n", k, i, j, i+1, j+1, W[i,j], alpha_old1[k], alpha_new[k]) ; 
        k ++; 
      }
    }
    
    /* Threshold alpha. Must lie in [-1, 1] range */ 
    for (int i = 0; i < c; i ++) { 
      if (alpha_new[i] > 1) { 
        alpha_new[i] = 1 ;  
      } 
      if (alpha_new[i] < -1) { 
        alpha_new[i] = -1 ;  
      } 
    }
    
    /* update beta and alpha for the next iteration step */
    for (int i = 0; i < m; i ++) { 
      *(beta_old + i) = *(beta_new + i) ; 
    }
    
    for (int i = 0; i < c; i ++) { 
      *(alpha_old2 + i) = *(alpha_old1 + i) ; 
      *(alpha_old1 + i) = *(alpha_new + i) ; 
    }
    //memcpy(beta_old, beta_new, sizeof(beta_new));
    //memcpy(alpha_old2, alpha_old1, sizeof(alpha_old1));
    //memcpy(alpha_old1, alpha_new, sizeof(alpha_new));
    
    
    // for (int i = 0; i < m; i ++) { 
    //   *(beta_old + i) = *(beta_new + i) ; 
    // }
    // 
    // for (int i = 0; i < c; i ++) { 
    //   *(alpha_old2 + i) = *(alpha_old1 + i) ; 
    //   *(alpha_old1 + i) = *(alpha_new + i) ; 
    // }
    // 
    // beta_old = Map<Eigen::VectorXd>(beta_new) ; 
    // alpha_old2 = Map<Eigen::VectorXd>(alpha_old1) ; 
    // alpha_old1 = Map<Eigen::VectorXd>(alpha_new) ; 
    // 
    // std::copy(beta_new.begin(), beta_new.end(), beta_old.begin()) ; 
    // std::copy(alpha_old1.begin(), alpha_old1.end(), alpha_old2.begin()) ; 
    // std::copy(alpha_new.begin(), alpha_new.end(), alpha_old1.begin()) ; 
    // 
    // 
    // double beta_new[m] ; // beta^(k + 1)
    // double delta [m] ; // aux. vector for beta-update step
    // 
    // /* initialize vectors for alpha-update step in the ADMM */
    // double alpha_new[c] ; 
    // double alpha[c] ; 
    //beta_old   = Rcpp::clone(beta_new) ; 
    //alpha_old2 = Rcpp::clone(alpha_old1) ; 
    //alpha_old1 = Rcpp::clone(alpha_new) ; 
    // 
    // double[] beta_new (m) ;
    // double[] alpha_new (c) ; 
    // 
    // double[] delta (m) ; // aux. vector for beta-update step
    // double[] alpha (c) ; // used to store (2*alpha^k - alpha^(k-1))
    // 
    
    iter ++; 
  }
  
  /* Turn to zero when really close */
  for (int i = 0; i < m; i ++) { 
    if (fabs(beta_new[i]) < truncate) { 
      *(beta_new + i) = 0 ;  
    } 
  }
  
  Rcout << "MAX ITER REACHED" << std::endl ; 
  //Eigen::VectorXd res(m);
  //std::copy(beta_new.begin(), beta_new.end(), res.begin()) ; 
  return(NumericVector(beta_new, beta_new + sizeof(beta_new) / sizeof(*beta_new))); 
}
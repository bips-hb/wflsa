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
 NumericVector genlassoRcpp(const NumericVector Y,
                            const NumericMatrix W,
                            const int m,
                            const double eta1,
                            const double eta2,
                            double a,
                            const double rho,
                            const int max_iter,
                            const double eps,
                            const double truncate) {


   // some frequently used constants
   a = rho*a ;
   const double C = 1 / (1 + a) ;

   // number of rows in matrix D
   const int c = (m*m + m) / 2 ;

   /* Compute y to a double array. I did this to make sure that the C compilation
    * environment does not have an effect */
   double* y = new double[m + 1] ;
   for (int i = 0; i < m; i ++) {
     *(y + i) = (double)(Y[i]) ;
   }

   /* initialize vectors for beta-update step in the ADMM */
   double *beta_new = new double[m+1];
   double *beta_old = new double[m+1] ;
   double *delta    = new double[m+1] ;

   for (int i = 0; i < m; i ++) {
     *(beta_new + i) = 0;
     *(beta_old + i) = 0;
     *(delta + i)    = 0;
   }

   /* initialize vectors for alpha-update step in the ADMM */
   double *alpha_new  = new double[c+1] ;
   double *alpha_old1 = new double[c+1] ;
   double *alpha_old2 = new double[c+1] ;
   double *alpha      = new double[c+1];

   for (int i = 0; i < c; i ++) {
     *(alpha_new + i)  = 0 ;
     *(alpha_old1 + i) = 0 ;
     *(alpha_old2 + i) = 0 ;
     *(alpha + i)      = 0 ;
   }

   /* Indices used for the alpha-update step */
   int steps[m-1] ;

   *steps = 0 ;

   for (int i = 1; i < (m-1); i ++) {
     *(steps + i) = m - i + *(steps + i - 1) ;
   }

   int iter = 0 ;    // number of iterations
   double diff = 0 ; // absolute difference between beta^(k+1) and beta^k

   /* loop until either the max. no. of iterations are reached
    * or the difference (diff) is smaller then eps
    */
   while (iter < max_iter) {

     /* ------- beta-update step ---------*/
     // alpha = 2*alpha_old1 - alpha_old2 ;
     for (int i = 0; i < c; i++) {
       *(alpha + i) = 2*(*(alpha_old1 + i)) - *(alpha_old2 + i) ;
     }

     // go over all possible pairs (i,j), same as D^T %*% (2 alpha^(k) - alpha^(k-1))
     for (int i = 0; i < m; i++) {
       *(delta + i) = eta1 * (*(alpha + i)) ;

       for (int j = i+1; j < m; j++) {
         *(delta + i) += eta2* W(i,j)*(*(alpha + m + steps[i] + (j - i) - 1)) ;
       }

       for (int j = 0; j < i; j++) {
         *(delta + i) -= eta2*W(i,j)*(*(alpha + m + steps[j] - (j - i) - 1)) ;
       }
     }

     diff = 0 ;
     for (int i = 0; i < m; i ++) {
       *(beta_new + i) = C*(a*(*(beta_old + i)) + y[i] - *(delta + i)) ;
       diff += fabs(*(beta_new + i) - *(beta_old + i)) ; //abs(beta_new[i] - beta_old[i]) ;
     }

     // determine whether converged or not
     if (diff < eps) {
       /* Turn to zero when really close */
       for (int i = 0; i < m; i ++) {
         if (fabs(*(beta_new + i)) < truncate) {
           *(beta_new + i) = 0 ;
         }
       }

       // convert results to NumericVector for R
       NumericVector res = NumericVector(beta_new, beta_new + m) ;

       // Clean-up ------------
       delete[] beta_old;
       delete[] alpha_old1;
       delete[] alpha_old2;
       delete[] alpha_new;
       delete[] alpha;
       delete[] y;
       delete[] beta_new ;

       return(res);
     }

     /* --------- alpha update step ----------- */
     for (int i = 0; i < m; i++) {
       *(alpha_new + i) = *(alpha_old1 + i) + rho * eta1 * *(beta_new + i) ;
     }

     int k = m;
     // go over all unique pairs (i,j)
     for (int i = 0; i < m-1; i++) {
       for (int j = i+1; j < m; j++) {
         *(alpha_new + k) = *(alpha_old1 + k) + rho * eta2 * W(i,j) * (*(beta_new + i) - *(beta_new + j)) ;
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

     iter ++;
   }

   /* Turn to zero when really close */
   for (int i = 0; i < m; i ++) {
     if (fabs(*(beta_new + i)) < truncate) {
       *(beta_new + i) = 0 ;
     }
   }

   // convert results to NumericVector for R
   NumericVector res = NumericVector(beta_new, beta_new + m) ;

   // Clean-up ------------
   delete[] beta_old;
   delete[] alpha_old1;
   delete[] alpha_old2;
   delete[] alpha_new;
   delete[] alpha;
   delete[] y;
   delete[] beta_new ;

   return(res);
 }

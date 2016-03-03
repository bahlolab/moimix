#include <Rcpp.h>
using namespace Rcpp;
// moimix - Model 1 Binomial Mixture Model (no error)
// Date: 03/03/2016
// Author: Stuart Lee
// Description: Goal is to remove dependency on flexmix and create fast
// modular EM algorithm that can work on either GDS object or 
// matrix of integer read counts.
// Functions required:
//      eStep - compute expectation step
//      mStep - update parameters based on eStep
//      Q - compute expected log-likelihood
//      moimix - do the work, output moimix object
//      clusters - hard assignments of SNPs to MOI categories (take in moimix object)
//      responsibility - soft assignments of SNPs to MOI categories (take in moimix object)

//' Compute eStep for binomial mixture model
//' @export
// [[Rcpp::export]]
// NumericMatrix eStep(IntegerVector x) {
//     return x;
// }

// double Q(NumericVector mu, NumericVector pi) {
//     return 0.0;
// }

NumericMatrix mStep(IntegerVector x) {
    return x;
}



// moimix moimix(IntegerVector x, IntegerVector N, int k, int niter, double eps) {
//     // starting parameters
//     NumericVector mu;
//     NumericVector pi;
//     double old_loglik;
//     double new_loglik;
//     bool converged;
//     
//     
//     for(int i=1; i <= N; i++) {
//         
//         if(abs(new_loglik - old_loglik) < eps) {
//             break
//         }
//     }
//     
//     
//     
// }




// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
# check eStep

# check mStep

# check moimix
*/

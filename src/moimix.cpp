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
NumericMatrix eStep(NumericVector x, NumericVector N, NumericVector mu, NumericVector pi) {
    int k = mu.size();
    int n = x.size();
    NumericMatrix tau(n, k);
    for(int i = 0; i < n; ++i) {
        for(int j = 0; j < k;  ++j) {
            tau(i, j) = pi(j) * R::dbinom(x(i), N(i), mu(j), false);
        }
    }
    // scaling factor
    for(int i = 0; i < n; ++i) {
        tau(i, _) = tau(i,_)/ sum(tau(i, _));
    }
    return tau;
}

// double Q(NumericVector mu, NumericVector pi) {
//     return 0.0;
// }

// // [[Rcpp::export]]
// List mStep(NumericVector x, NumericVector N, NumericMatrix tau) {
//     
//     int k = tau.size(1);
//     List out(k);
//     
//     
//      
// }



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
# check mStep, takes in suffucient statistics and outputs parameter estimates
mu.true <- c(0.3, 0.7)
pi.true <- c(0.8, 0.2)
# generate a simple mixture of binomials to test
mix2 <- sampleMM(1000, 30, 2, mu.true, pi.true)
mu.init <- c(0.5, 0.5)
pi.init <- c(0.5, 0.5)
priors <- eStep(mix2$obs, rep(30, 1000), mu.init, pi.init)

# check eStep


# check moimix
*/

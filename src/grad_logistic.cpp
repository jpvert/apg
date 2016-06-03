#//#include <Rcpp.h>
#include <RcppArmadillo.h>
#include <math.h>
#include <R.h>
#include <Rdefines.h>
#include <iostream>
#include <stdio.h>
using namespace Rcpp;
using namespace std;

//' Gradient of the logistic regression risk.
//'
//' Computes the gradient of the logistic regression error: \deqn{f(x) = 1/n
//' sum_{i=1}^n ( log(1 + exp(- b_i * x' * A_i) ) ,} for a given design matrix
//' \code{A} and a response binary vector \code{b} with values \code{+1} and
//' \code{-1}.
//'
//' @param x A p-dimensional vector where the gradient is computed.
//' @param opts List of parameters, which must include: \itemize{ \item \code{A},
//'   a n*p design matrix \item \code{b}, a n-dimensional response vector taking
//'   values \eqn{1} or \eqn{-1}. }
//'
//' @return The gradient of the function \eqn{f(x) = 1/n sum_{i=1}^n ( log(1 +
//'   exp(- b_i * x' * A_i) )}, which is: \deqn{f'(x) = 1/n sum_{i=1}^n -b_i A_i
//'   / (1 + exp( b_i * x' * A_i ))}
//'
//' @export
//' @examples grad.logistic(c(1,3,-2), list(A=diag(3), b=c(1,-1,1))
//' # The following give the same (elastic net penalized logistic regression without offset):
//' glmnet(A,b,family="binomial",lambda=0.1,standardize=FALSE,intercept=FALSE, alpha=0.5)$b
//' apg(grad.logistic, prox.elasticnet, ncol(A), list(A=A, b=b, lambda=0.1, alpha=0.5) )$x
//' 
// [[Rcpp::export]]
arma::mat grad_logistic(arma::vec& x, List opts) {
    // TYPEOF allows to get the type of a SEXP object.
    // cout << TYPEOF(opts['A']) << endl;
    // opts['A'] is a SEXP object of SEXPTYPE LANGSXP (language objects)
    // i.e this is a call.
    
    // A list is represented by a SEXP object of SEXPTYPE VECSXP.
    // cout << TYPEOF(VECTOR_ELT(opts, 1)) << endl;
    // The two elements of our list are SEXP objects of SEXPTYPE REALSXP (numeric vector)
    
    double* pA;
    SEXP sexpA = VECTOR_ELT(opts, 0);
    pA = REAL(sexpA);
    // cout << pA << endl;
    SEXP dimMatrix = Rf_getAttrib(sexpA, R_DimSymbol);
    int nrow = INTEGER(dimMatrix)[0];
    int ncol = INTEGER(dimMatrix)[1];
    arma::mat A(pA, nrow, ncol, false, true);
    arma::vec b = as<arma::vec>(opts["b"]);
    arma::vec s = 1/(1+arma::exp((A*x)%b));
    arma::mat grad = A.t()*(s%b);

    return -grad/nrow;
}








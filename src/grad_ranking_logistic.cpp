#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

//' @useDynLib apg
//' @importFrom Rcpp sourceCpp

NumericVector matrix_vector_mult(NumericMatrix A, NumericVector x) {
    // Simple matrix-vector multiplication
    int nrow = A.nrow(), ncol = A.ncol();
    NumericVector b(nrow);

    for (int i = 0; i < nrow; i++){
        b[i] = 0;
        for (int j = 0; j < ncol; j++) {
            b[i] += A(i, j) * x[j];
        }
    }
    return b;
}

//' Gradient of the ranking logistic regression risk.
//'
//' Computes the gradient of the ranking logistic regression error: \deqn{f(x) = 1/|comp|
//' sum_{(i,j)\in comp} ( log(1 + exp(- x' * (A_j - A_i) )} for a given design matrix
//' \code{A} and a set of pairs \eqn{(i,j)\in P} where the j-th sample is larger than the i-th sample.
//'
//' @param x A p-dimensional vector where the gradient is computed.
//' @param opts List of parameters, which must include: \itemize{ \item \code{A},
//'   a n*p design matrix \item \code{comp}, a list of length n, where the i-th entry is a vector containing the indices of the samples j which are larger than sample i.}
//'
//' @return The gradient of the function \eqn{f(x) = 1/|comp|
//' sum_{(i,j)\in comp} ( log(1 + exp(- x' * (A_j - A_i) )}, which is: \deqn{f'(x) = 1/|comp| sum_{(i,j)\in comp} (A_j - A_i) / (1 + exp( x' * ( A_i - A_j) )}
//'
//' @export
//' @examples grad.rankinglogistic(c(1,3,-2), list(A=diag(3), comp=list(as.integer(c(2,3)),as.integer(3),integer(length=0))))
// [[Rcpp::export(name="grad.rankinglogistic")]]
NumericVector grad_ranking_logistic(NumericVector x, List opts) {
    // Gradient of the ranking logistic regression

    // The design matrix (each row is a sample)
    NumericMatrix A = as<NumericMatrix>(opts["A"]);

    // The pairwise comparisons between samples (a list of vectors, the i-th vector in the list contains the indices of the samples which are larger than the i-th sample)
    Rcpp::List jlist = opts["comp"];

    int n=jlist.size();
    Rcpp::NumericVector weight(n);
    int p = A.ncol();
    int N = 0;
    double tmp;
    NumericVector grad(p);
    NumericMatrix grad_mat(p, 1);

    // Compute exp(A*x)
    NumericVector s = matrix_vector_mult(A, x);
    for (int i=0; i<n; i++)
        s[i] = exp(s[i]);

    // Compute the weight of each sample in the gradient
    for (int i = 0; i < n; i++){
        SEXP ll=jlist[i];
        Rcpp::NumericVector y(ll);
        N += y.size(); // total number of pairs
        for (int j=0; j<y.size(); j++){
//            tmp = 1/(1+exp(s[y[j]-1] - s[i]));
            tmp = 1/(1+ s[y[j]-1] / s[i]);
            weight[i] +=tmp;
            weight[y[j]-1] -=tmp;
        }
    }

    // Compute the gradient as a linear combination of the samples
    for (int j = 0; j < p; j++) {
        for (int i=0; i<n; i++) {
            grad[j] += weight[i]*A(i,j);
        }
        grad[j] /= N;
    }

    grad_mat(_, 0) = grad;

    return grad_mat;
}

#' Gradient of the logistic regression risk.
#'
#' Computes the gradient of the logistic regression error: \deqn{f(x) = 1/n
#' sum_{i=1}^n ( log(1 + exp(- b_i * x' * A_i) ) ,} for a given design matrix
#' \code{A} and a response binary vector \code{b} with values \code{+1} and
#' \code{-1}.
#'
#' @param x A p-dimensional vector where the gradient is computed.
#' @param opts List of parameters, which must include: \itemize{ \item \code{A},
#'   a n*p design matrix \item \code{b}, a n-dimensional response vector taking
#'   values \eqn{1} or \eqn{-1}. }
#'
#' @return The gradient of the function \eqn{f(x) = 1/n sum_{i=1}^n ( log(1 +
#'   exp(- b_i * x' * A_i) )}, which is: \deqn{f'(x) = 1/n sum_{i=1}^n -b_i A_i
#'   / (1 + exp( b_i * x' * A_i ))}
#'
#' @export
#' @examples grad.logistic(c(1,3,-2), list(A=diag(3), b=c(1,-1,1))
#' # The following give the same (elastic net penalized logistic regression without offset):
#' glmnet(A,b,family="binomial",lambda=0.1,standardize=FALSE,intercept=FALSE, alpha=0.5)$b
#' apg(grad.logistic, prox.elasticnet, ncol(A), list(A=A, b=b, lambda=0.1, alpha=0.5) )$x


grad.logistic <- function(x, opts) {

    sA <- opts[["A"]]*opts[["b"]] # The sampled multiplied by the response

    return( crossprod( -sA , 1 / (1 + exp(sA %*% x))) / nrow(opts[["A"]]) )
}

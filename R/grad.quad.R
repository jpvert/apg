#' Gradient of the squared error.
#'
#' Computes the gradient of the sum of squared error: \deqn{f(x) = 1/(2n)
#' \sum_{i=1}^n ( x'*A_i - b_i)^2}.
#'
#' @param x A p-dimensional vector where the gradient is computed.
#' @param opts List of parameters, which must include: \itemize{ \item \code{A},
#'   a n*p design matrix, where each row is a sample and each column is a
#'   variable \item \code{b}, a n-dimensional response vector. }
#'
#' @return The gradient of the function \eqn{f(x) = 1/(2n) || Ax - b ||^2},
#'   which is \eqn{A'*(Ax - b)/n}.
#'
#' @export
#' @examples grad.quad(c(1,3,-2), list(A=diag(3), b=rep(1,3)))


grad.quad <- function(x, opts) {

    # The gradient is simply t(A) %*% ( A %*% x - b ) / n

    crossprod(opts[["A"]] , opts[["A"]] %*% x - opts[["b"]]) / nrow(opts[["A"]])

}

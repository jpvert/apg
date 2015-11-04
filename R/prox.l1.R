#' Proximal operator of the scaled L1 norm.
#'
#' Computes the proximal operator of the L1 norm: \deqn{h(x) = \lambda ||x||_1
#' ,} where \eqn{\lambda} is a scaling factor.
#'
#' @param x The input vector
#' @param t The step size
#' @param opts List of parameters, which must include:
#'   \itemize{\code{opts$lambda}, the scaling factor of the L1 norm.}
#'
#' @return The proximal of \eqn{h} at {x} with step size \eqn{t}, given by
#'   \deqn{prox_h(x,t) = argmin_u [ t h(u) + 1/2 || x - u ||^2 ]}.
#'
#' @export
#' @examples prox.l1(c(1,3,-2), 1.5, list(lambda=1))


prox.l1 <- function(x, t, opts) {
    ## Compute the soft-thresholded operator
    thres <- t * opts$lambda
    idx.1 <- which(x < -thres)
    idx.2 <- which(x > thres)
    res <- rep(0,length(x))
    if ( length(idx.1)>0 ) res[idx.1] <- x[idx.1] + thres
    if ( length(idx.2)>0 ) res[idx.2]<- x[idx.2] - thres
    return(res)
}

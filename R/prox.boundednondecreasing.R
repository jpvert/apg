#' Proximal operator for the set of bounded non-decreasing vectors
#' 
#' Computes the proximal operator of the bounded non-decreasing vectors, i.e., 
#' \deqn{prox(x) = argmin_u || u - x ||^2 such that u_1 \le ... \le u_n and 
#' ||u|| \le M} This is simply the projection of $x$ onto the set of
#' non-decreasing vectors, obtained by isotonic regression, followed by
#' projection onto the Euclidean ball of radius $M$.
#' 
#' @param x The input vector
#' @param opts Optional list of parameters, which can include: \itemize{ \item
#'   \code{maxnorm} : the bound on the Euclidean norm of the non-decreasing
#'   vector (default \code{1}).}
#'   
#' @return The projection of \code{x} onto the set of nondecreasing vectors with
#'   norm bounded by \code{maxnorm}.
#'   
#' @export
#' 
#' @examples 
#' a <- rnorm(20)+1
#' plot(a)
#' lines(prox.boundednondecreasing(a), col=2)
#' lines(prox.boundednondecreasing(a, opts=list(maxnorm=10)), col=3)


prox.boundednondecreasing <- function(x, t=0, opts=list()) {
    
    if (exists("maxnorm",where=opts)) {
        M <- opts[["maxnorm"]]
    } else {
        M <- 1
    }
    # We first project onto the set of non-decreasing vectors using isotonic
    # regression, then onto the unit Euclidean ball.
    u <- pava(x)
    unorm <- sqrt(sum(u^2))
    if (unorm > M) {
        return(M*u/unorm)
    } else {
        return(u)
    }
}
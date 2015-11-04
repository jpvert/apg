#' Proximal operator of the isotonic constraint
#'
#' Computes the proximal operator of the isotonic constraint, i.e., the
#' projection onto the set of nondecreasing vectors.
#'
#' @param x The input vector
#'
#' @return The projection of \code{x} onto the set of nondecreasing vectors,
#'   obtained by solving an isotonic regression.
#'
#' @export
#' @examples prox.isotonic(c(1,3,-2,4,5))


prox.isotonic <- function(x, ...) {
    return(pava(x))
}

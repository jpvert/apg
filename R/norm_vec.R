#' Euclidean norm of a vector
#' 
#' Computes the Euclidean norm of a vector
#' 
#' @param x The input vector
#' 
#' @return The Euclidean norm of the input vector
#' 
#' @examples norm_vec(c(1,1))
norm_vec <- function(x) {
    sqrt(sum(x^2))
}

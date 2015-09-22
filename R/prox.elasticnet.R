#' Proximal operator of the scaled elastic net penalty.
#' 
#' Computes the proximal operator of the scaled elastic net penalty: \deqn{h(x) 
#' = \lambda [ (1 - \alpha)/2 ||x||_2^2 + \alpha ||x||_1 ] ,} where 
#' \eqn{\lambda} is a scaling factor and \eqn{\alpha \in [0,1]} balances between
#' the L1 and L2 norms.
#' 
#' @param x The input vector
#' @param t The step size (default is \code{1})
#' @param opts List of parameters, which can include: \itemize{ \item 
#'   \code{lambda} : the scaling factor of the L1 norm (default is
#'   \code{lambda=1})\item \code{alpha} : the balance between L1 and L2 norms.
#'   \code{alpha=0} is the squared L2 (ridge) penalty, \code{alpha=1} is the L1
#'   (lasso) penalty. Default is \code{alpha=1} (lasso).}
#'   
#' @return The proximal of \eqn{h} at {x} with step size \eqn{t}, given by 
#'   \deqn{prox_h(x,t) = argmin_u [ t h(u) + 1/2 || x - u ||^2 ]}.
#'   
#' @export
#' 
#' @examples prox.elasticnet(c(1,3,-2), 1.5, list(lambda=1,alpha=0.5))


prox.elasticnet <- function(x, t=1, opts=list(lambda=1, alpha=1)) {
    
    if (!exists("lambda",where=opts)) {
        lambda <- 1
    } else {
        lambda <- opts[["lambda"]]
    }

    if (!exists("alpha",where=opts)) {
        alpha <- 1
    } else {
        alpha <- opts[["alpha"]]
    }
    
    ## Compute the soft-thresholded operator
    thres <- t * lambda * alpha
    idx.1 <- which(x < -thres)
    idx.2 <- which(x > thres)
    res <- rep(0,length(x))
    if ( length(idx.1)>0 ) res[idx.1] <- x[idx.1] + thres
    if ( length(idx.2)>0 ) res[idx.2]<- x[idx.2] - thres
    return(res / (1 + t * lambda * (1 - alpha)))
}
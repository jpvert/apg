#' Proximal operator of the group lasso penalty
#'
#' Computes the proximal operator of the group lasso penalty: \deqn{h(x) =
#' \sum_{g group} w_g ||x[g]||_2 .} Note that the groups should not
#' overlap.
#'
#' @param x The input vector
#' @param t The step size
#' @param opts List of parameters, which can include: \itemize{ \item
#'   \code{groups} : a list of groups, each group is just a sequence of indices
#'   of the components that form the group (default: all singletons). \item \code{groupweigths} : a vector of weights for the groups. If a single number, all groups have the same weight (default \code{1})
#'   }

#' @return The proximal operator of the group lasso, which is a soft-thresholing
#'   operator applied to the restriction of the \code{x} to each group.
#'
#' @export
#' @examples
#' x <- rnorm(5)
#' # When groups are all the singletons we recover the L1 (lasso) penalty
#' prox.grouplasso(f,1,list(groups=as.list(seq(length(f)))))
#' prox.elasticnet(f,1,list(lambda=1,alpha=1))


prox.grouplasso <- function(x, t, opts=list(groups=as.list(seq(length(xc))))) {

    if (!exists("groups",where=opts))
        stop("No list of groups provided for the group lasso.")
    ngroups <- length(opts$groups)
    if (!exists("groupweights",where=opts)) {
        w <- rep(t, ngroups)
    } else {
        if (length(opts[["groupweights"]]) == ngroups) {
            w <- t*opts[["groupweights"]]
        } else {
            w <- t*rep(opts[["groupweights"]][1], ngroups)
        }
    }

    u <- x
    for (i in seq(ngroups)) {
        g <- opts[["groups"]][[i]]
        u[g] <- max(0, 1 - w[i] / norm_vec(x[g]) ) * x[g]
    }
    return(u)
}

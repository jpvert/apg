#' Fit a generalized linear model with various regularizations.
#'
#' Fit a generalized linear model via penalized maximum likelihood. Fits linear
#' and logistic regression models, with elastic net or isotonic regularization.
#'
#' @param x The input matrix, each row is a sample, each column a feature.
#' @param y The response variable. Quantitative for \code{family="gaussian"}, binary with values \code{+1} and \code{-1} for \code{family="binomial"}, survival data with two columes (follow-up time and event) for \code{family="survival"}
#' @param family The response type. For \code{family="gaussian"}, ...
#' @param penalty The penalty type.
#' @param lambda The scaling of the penalty (default \code{1})
#' @param intercept Should intercept(s) be fitted (default=\code{TRUE}) or set to zero (\code{FALSE})
#' @param opts List of parameters, which must include: \itemize{ \item }
#'
#' @return toto
#'
#' @export
#' @examples
#' n <- 100
#' p <- 5
#' x <- matrix(rnorm(n*p),n,p)
#' y <- rbinom(n,1,0.5)*2-1
#' lambda <- 0.2*max(abs(crossprod(y,x)))/n
#'
#' # Lasso regression with intercept:
#' m <- glm.apg(x, y, lambda=lambda)
#' # same as
#' library(glmnet)
#' m2 <- glmnet(x, y, standardize=FALSE, lambda=lambda)
#'
#' # Ridge regression with intercept:
#' m <- glm.apg(x, y, lambda=lambda, opts=list(alpha=0))
#' # Does the same as
#' m2 <- glmnet(x, y, standardize=FALSE, lambda=lambda, alpha=0)
#'
#' # Elastic net regression with intercept:
#' m <- glm.apg(x, y, lambda=lambda, opts=list(alpha=0.5))
#' # Does the same as
#' m2 <- glmnet(x, y, standardize=FALSE, lambda=lambda, alpha=0.5)
#'
#' # Elastic net regression without intercept:
#' m <- glm.apg(x, y, lambda=lambda, intercept=FALSE, opts=list(alpha=0.5))
#' # Does the same as
#' m2 <- glmnet(x, y, standardize=FALSE, lambda=lambda, alpha=0.5, intercept=FALSE)
#'
#' # Lasso penalized logistic regression with intercept:
#' m <- glm.apg(x, y, family="binomial", lambda=lambda)
#' # Does the same as
#' m2 <- glmnet(x, y, family="binomial", lambda=lambda, standardize=FALSE)
#'
#' # Elastic net penalized logistic regression with intercept:
#' m <- glm.apg(x, y, family="binomial", lambda=lambda, opts=list(alpha=0.5))
#' # Does the same as
#' m2 <- glmnet(x, y, family="binomial", lambda=lambda, standardize=FALSE, alpha=0.5)
#'
#' # Isotonic regression with offset
#' m <- glm.apg(x, y, penalty="isotonic", lambda=lambda)
#'
#' # Isotonic logistic regression with offset
#' m <- glm.apg(x, y, family="binomial", penalty="isotonic", lambda=lambda)
#'
#' # Isotonic logistic regression with offset, with non-decreasing model of bounded norm
#' m <- glm.apg(x, y, family="binomial", penalty="boundednondecreasing", lambda=lambda, opts=list(maxnorm=2))

glm.apg <- function(x, y, family=c("gaussian", "binomial", "survival"), penalty=c("elasticnet", "isotonic", "boundednondecreasing"), lambda=1, intercept=TRUE, opts=list()) {
    family <- match.arg(family)
    if (family=="survival") {
        intercept=FALSE
    }
    penalty <- match.arg(penalty)
    y <- drop(y)
    np <- dim(x)
    if (is.null(np) | (np[2] <= 1))
        stop("x should be a matrix with 2 or more columns")
    n = as.integer(np[1])
    p = as.integer(np[2])
    dimy = dim(y)
    nrowy = ifelse(is.null(dimy), length(y), dimy[1])
    if (nrowy != n)
        stop(paste("number of observations in y (", nrowy, ") not equal to the number of rows of x (",
                   n, ")", sep = ""))
    vnames = colnames(x)
    if (is.null(vnames))
        vnames = paste("V", seq(n), sep = "")

    # Gradient of the smooth part
    gradG <- switch(family, gaussian = grad.quad, binomial = grad.logistic, survival = grad.rankinglogistic)

    # Prox of the nonsmooth part
    proxH <- switch(penalty, elasticnet = prox.elasticnet, isotonic = prox.isotonic, boundednondecreasing = prox.boundednondecreasing)

    o <- opts

    # If a non-penalized intercept is added, just add a constant column to x and modify the prox operator of the penalty to not touch the coefficient corresponding to the constant column. We also work on the centered matrix x to speed up convergence.
    if (intercept) {
        centered.x <- scale(x, scale = FALSE)
        o$A <- cbind(centered.x, rep(1,n))
        myproxH <- function(u, ...) {
            return(c(proxH(u[-length(u)], ...), u[length(u)]))
        }
    } else {
        o$A <- x
        myproxH <- proxH
    }

    if (family=="survival") {
        o$comp <- surv_to_pairs(y[,1],y[,2])
    }
    o$b <- y
    o$lambda <- lambda

    # Maximize the penalized log-likelihood
    res <- apg(gradG, myproxH, ncol(o$A), o)
    w <- res[["x"]]

    # Return the model (vector of weight in b, intercept in a0)
    if (intercept) {
        return(list(b=w[-length(w)], a0=w[length(w)] - sum(w[-length(w)] * attr(centered.x, "scaled:center"))))
    } else {
        return(list(b=w, a0=0))
    }
 }

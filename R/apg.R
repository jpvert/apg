#' Accelerated proximal gradient optimization
#'
#' This function implements an accelerated proximal gradient method (Nesterov
#' 2007, Beck and Teboulle 2009). It solves: \deqn{min_x (f(x) + h(x)), x \in
#' R^dim_x} where \eqn{f} is smooth, convex and \eqn{h} is non-smooth, convex
#' but simple in that we can easily evaluate the proximal operator of h.
#'
#' @param grad_f A function that computes the gradient of f :
#'   \eqn{grad_f(v,opts) = df(v)/dv}
#' @param prox_h A function that computes the proximal operator of h :
#'   \eqn{prox_h(v,t,opts) = argmin_x (t*h(x) + 1/2 * norm(x-v)^2)}
#' @param dim_x The dimension of the unknown \eqn{x}
#' @param opts List of parameters, both for the \code{apg} function and for the
#'   \code{grad_f} and \code{prox_h} functions. For \code{apg}, the list can
#'   contain the following fields: \itemize{ \item \code{X_INIT}: initial
#'   starting point (default \code{0}) \item \code{USE_RESTART} : use adaptive
#'   restart scheme (default \code{TRUE}) \item \code{MAX_ITERS} : maximum
#'   iterations before termination (default \code{2000}) \item \code{EPS} :
#'   tolerance for termination (default \code{1e-6}) \item \code{ALPHA} :
#'   step-size growth factor (default \code{1.01}) \item \code{BETA} : step-size
#'   shrinkage factor (default \code{0.5}) \item \code{QUIET} : if false writes
#'   out information every 100 iters (default \code{FALSE}) \item
#'   \code{GEN_PLOTS} : if true generates plots of norm of proximal gradient
#'   (default \code{TRUE}) \item \code{USE_GRA} : if true uses UN-accelerated
#'   proximal gradient descent (typically slower) (default \code{FALSE}) \item
#'   \code{STEP_SIZE} : starting step-size estimate, if not set then apg makes
#'   initial guess (default \code{NULL}) \item \code{FIXED_STEP_SIZE} : don't
#'   change step-size (forward or back tracking), uses initial step-size
#'   throughout, only useful if good STEP_SIZE set (default \code{FALSE}) } In
#'   addition, \code{opts} can contain other fields that will be passed to the
#'   \code{grad_f} and \code{prox_h} functions. This code was borrowed and
#'   adapted from the MATLAB version of Brendan O'Donoghue available at
#'   \code{https://github.com/bodono/apg}
#'
#' @return A list with \code{x}, the solution of the problem: \deqn{min_x (f(x)
#'   + h(x)), x \in R^dim_x ,} and \code{t}, the last step size.
#'
#' @export
#' @examples # Solve a Lasso problem:
#' # min_x 1/2 norm( A%*%x - b )^2 + lambda ||x||_1
#' n <- 50
#' m <- 20
#' lambda <- 1
#' A <- matrix(rnorm(m*n), nrow=n)
#' b <- rnorm(n)
#' r <- apg(grad.quad, prox.l1, m, list(A=A, b=b, lambda=lambda) )
#' # This gives the same result as:
#' # m <- glmnet(A,b,alpha=1, standardize=FALSE,lambda=1/50,intercept=FALSE)


apg <- function(grad_f, prox_h, dim_x, opts) {

    # Set default parameters
    X_INIT <- numeric(dim_x) # initial starting point
    USE_RESTART <- TRUE # use adaptive restart scheme
    MAX_ITERS <- 2000 # maximum iterations before termination
    EPS <- 1e-6 # tolerance for termination
    ALPHA <- 1.01 # step-size growth factor
    BETA <- 0.5 # step-size shrinkage factor
    QUIET <- FALSE # if false writes out information every 100 iters
    GEN_PLOTS <- TRUE # if true generates plots of norm of proximal gradient
    USE_GRA <- FALSE # if true uses UN-accelerated proximal gradient descent (typically slower)
    STEP_SIZE = NULL # starting step-size estimate, if not set then apg makes initial guess
    FIXED_STEP_SIZE <- FALSE # don't change step-size (forward or back tracking), uses initial step-size throughout, only useful if good STEP_SIZE set

    # Replace the default parameters by the ones provided in opts if any
    for (u in c("X_INIT","USE_RESTART", "MAX_ITERS", "EPS", "ALPHA", "BETA", "QUIET", "GEN_PLOTS", "USE_GRA", "STEP_SIZE", "FIXED_STEP_SIZE")) {
        eval(parse(text=paste('if (exists("',u,'", where=opts)) ',u,' <- opts[["',u,'"]]',sep='')))
    }

    # Initialization
    x <- X_INIT
    y <- x
    g <- grad_f(y, opts)
    theta <- 1

    # Initial step size
    if (is.null(STEP_SIZE)) {

        # Barzilai-Borwein step-size initialization:
        t <- 1 / norm_vec(g)
        x_hat <- x - t*g
        g_hat <- grad_f(x_hat, opts)
        t <- abs(sum( (x - x_hat)*(g - g_hat)) / (sum((g - g_hat)^2)))
    } else {
        t <- STEP_SIZE
    }

    # Main loop
    for (k in seq(MAX_ITERS)) {

        if (!QUIET && (k %% 100==0)) {
            message(paste('iter num ',k,', norm(tGk): ',err1,', step-size: ',t,sep=""))
        }

        x_old <- x
        y_old <- y

        # The proximal gradient step (update x)
        x <- prox_h( y - t*g, t, opts)

        # The error for the stopping criterion
        err1 <- norm_vec(y-x) / max(1,norm_vec(x))
        if (err1 < EPS) break

        # Update theta for acceleration
        if(!USE_GRA)
            theta <- 2/(1 + sqrt(1+4/(theta^2)))
        else
            theta <- 1
        end

        # Update y
        if (USE_RESTART && sum((y-x)*(x-x_old))>0) {
            x <- x_old
            y <- x
            theta <- 1
        } else {
            y <- x + (1-theta)*(x-x_old)
        }

        # New gradient
        g_old <- g
        g <- grad_f(y,opts)

        # Update stepsize by TFOCS-style backtracking
        if (!FIXED_STEP_SIZE) {
            t_hat <- 0.5*sum((y-y_old)^2)/abs(sum((y - y_old)*(g_old - g)))
            t <- min( ALPHA*t, max( BETA*t, t_hat ))
        }
    }
    if (!QUIET) {
        message(paste('iter num ',k,', norm(tGk): ',err1,', step-size: ',t,sep=""))
        if (k==MAX_ITERS) message(paste('Warning: maximum number of iterations reached'))
        message('Terminated')
    }

    # Return solution and step size
    return(list(x=x,t=t))
}

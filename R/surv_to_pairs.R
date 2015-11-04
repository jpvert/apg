#' Given survival information (time and status) for a set of subjects, compute all pairs (i,j) for which patient i has the event before patient j.
#'
#' @param time A vector of follow up time for the subjects
#' @param event A vector of status indicator for the subjects, normally \code{0}=alive, \code{1}=dead. Other choices are \code{TRUE/FALSE} (\code{TRUE}=death). By default, all subjects are assumed to have an event.
#'
#' @return  A list of vectors, where the i-th vector is the set of patients that have the events after patient i.
#' @export
#' @examples
#' time <- c(5,6,7,7,8,9)
#' event <- c(1,0,0,1,1,0)
#' res <- surv_to_pairs(time,event)

surv_to_pairs <- function(time, event=rep(1, length(time))) {

    if (missing(time))
        stop("Must have a time argument")
    if (is.logical(event))
        event <- as.numeric(event)

    l <- lapply (seq(length(time)), function(i){
        which((event[i]==1) & ((time>time[i]) | ((time==time[i]) & !event) ))
    })
    return(l)
}

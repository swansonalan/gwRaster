

#' Test for a vailid doParallel cluster
#'
#' This function tests to see if a dummy function can be executed over a doParallel cluster using
#' 'foreach() %dopar% ...' syntax.  Returns T or F based on whether the foreach executed without
#' error.  If verbose=T, the function prints the number of cores detected.
#'
#' @param verbose whether or not to print cluster size
#' @return logical; whether or not the foreach executes successfully over a cluster
#' @export
test_cluster <- function(verbose=T) {
    # x <- !inherits(nowarnings(try(foreach::foreach(i=1:3,.combine='c') %dopar% i)),"try-error")
    x <- !inherits(try(foreach::foreach(i=1:3,.combine='c') %dopar% i),"try-error")
    y <- foreach::getDoParRegistered()
    n <- foreach::getDoParWorkers()
    z <- x & y & n>1
    if(z) {if(verbose) cat("found",n,"core cluster\n")} else {if(verbose) cat("no valid cluster found\n")}
    z
}

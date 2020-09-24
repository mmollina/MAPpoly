#' Access a remote server to get Counts for recombinant classes
#'
#' @param void internal function to be documented
#' @keywords internal
#' @import RCurl
#' @export get_cache_two_pts_from_web
get_cache_two_pts_from_web <- function(m, url.address = NULL, joint.prob = TRUE, verbose = FALSE) {
    if (is.null(url.address)) {
        if (m == 2)
            pl <- "di" else if (m == 4)
                pl <- "autotetra" else if (m == 6)
                    pl <- "autohexa" else if (m == 8)
                        pl <- "autoocta" else if (m == 10)
                            pl <- "autodeca" else if (m == 12)
                                pl <- "autododeca" else if (m == 14)
                                    pl <- "autotetradeca" else stop("Cache file not found. The ploidy level should be one of the following:\n 2 4 6 8 10 12 14.")
                                url.address.cond = paste("http://152.1.45.19/prev.cache.", pl, "ploid.counts.RData", sep = "")
                                url.address.joint = paste("http://152.1.45.19/prev.joint.cache.", pl, "ploid.counts.RData", sep = "")
    }
    prev.joint.cache.from.web <- prev.cache.from.web <- NULL
    ## checking internet connectivity
    if (try(url.exists(url.address.cond)) & try(url.exists(url.address.joint))) {
        if (verbose) cat("Internet connectivity ok.\nLoading genotype counts from web\n")
        load(url(url.address.cond, method = "libcurl"))
        if (!inherits(prev.cache.from.web, "cache.info"))
            stop(deparse(substitute(prev.cache.from.web)), " is not an object of class 'cache.info'\nTry to use function 'cache.two.pts'")
        load(url(url.address.joint, method = "libcurl"))
        if (!inherits(prev.cache.from.web, "cache.info"))
            stop(deparse(substitute(prev.cache.from.web)), " is not an object of class 'cache.info'\nTry to use function 'cache.two.pts'")
    } else stop("Cache file not found. Try to use get.from.web = FALSE")
    structure(list(cond = prev.cache.from.web, joint = prev.joint.cache.from.web), class = "cache.info")
}

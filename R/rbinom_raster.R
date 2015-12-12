#' Sample Binomial Values in a Raster*
#'
#' Returns a Raster* object whose values are random, Binomial(N,prob.seed).
#'
#' @param size Number, or Raster* object of sample sizes.
#' @param prob Number, or Raster* object of probabilities.
#' @export
#' @return A Raster* object of the same form as the input.
#' At least one of \code{size} or \code{prob} must be Raster* objects.
rbinom_raster <- function (size, prob) {
    if (!inherits(size,"Raster") && !inherits(prob,"Raster")) {
        stop("Either size or prob must be Raster objects.")
    }
    if (!inherits(size,"Raster")) {
        out <- prob
        prob <- raster::values(prob)
        size <- rep_len(size,length(out))
    }
    if (!inherits(prob,"Raster")) {
        out <- size
        size <- raster::values(size)
        prob <- rep_len(prob,length(out))
    }
    nonzero <- ( !is.na(prob) ) & ( !is.na(size) ) & ( prob>0 ) & ( size>0 )
    raster::values(out)[nonzero] <- rbinom(sum(nonzero),size=size[nonzero],prob=prob[nonzero])
    return(out)
}

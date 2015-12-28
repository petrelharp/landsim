#' Sample Poisson Raster
#'
#' Returns a Raster* which is Poisson-values with mean given by the input Raster*.
#'
#' @param lambda Raster* object of means.
#' @param only.values Return only the vector of non-NA values?
#' @export
#' @return A Raster* object of the same form as the input,
#' or if \code{only.values} is \code{TRUE} the vector of values.
rpois_raster <- function (lambda,only.values=FALSE) {
    nonmissing <- !is.na(raster::values(lambda))
    x <- rpois(sum(nonmissing),lambda=raster::values(lambda)[nonmissing])
    if (only.values) {
        return(x)
    } else {
        out <- lambda
        raster::values(out)[nonmissing] <- x
        return(out)
    }
}

#' Sample from a Poisson, Preserving Dimension
#'
#' Works as \code{rpois}, but returns an object of the same form as the input.
#'
#' @param lambda A numeric object whose values can be assigned to as "lambda[] <- ...".
#' @export
#' @return An object of the same form as the input, with integer values.
rpois_matrix <- function (lambda) {
    lambda[] <- rpois(length(lambda),as.numeric(lambda))
    return(lambda)
}


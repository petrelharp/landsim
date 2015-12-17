#' Sample Poisson Raster
#'
#' Returns a Raster* which is Poisson-values with mean given by the input Raster*.
#'
#' @param layer Raster* object of means.
#' @param only.values Return only the vector of non-NA values?
#' @export
#' @return A Raster* object of the same form as the input,
#' or if \code{only.values} is \code{TRUE} the vector of values.
rpois_raster <- function (lambda,only.values=FALSE) {
    nonmissing <- !is.na(values(lambda))
    x <- rpois(sum(nonmissing),lambda=values(lambda)[nonmissing])
    if (only.values) {
        return(x)
    } else {
        out <- lambda
        values(out)[nonmissing] <- x
        return(out)
    }
}

#' Sample from a Poisson, Preserving Dimension
rpois_matrix <- function (lambda) {
    lambda[] <- rpois(length(lambda),as.numeric(lambda))
    return(lambda)
}


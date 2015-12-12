#' Sample Poisson Raster
#'
#' Returns a Raster* which is Poisson-values with mean given by the input Raster*.
#'
#' @param layer Raster* object of means.
#' @export
#' @return A Raster* object of the same form as the input.
rpois_raster <- function (layer) {
    out <- layer
    nonzero <- ( !is.na(values(layer)) ) & ( values(layer)>0 )
    values(out)[nonzero] <- rpois(sum(nonzero),lambda=values(layer)[nonzero])
    return(out)
}

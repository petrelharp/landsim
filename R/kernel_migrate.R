#' Weighted "Adjacency" Matrix from a Raster Layer.
#'
#' This returns the (sparse) Matrix giving the pseudo-adjacency matrix 
#' for cells in \code{layer}, with weights.
#'
#' @inheritParams layer_adjacency
#' @param layer Raster* object consisting of layers to (independently) smooth.
#' @param fun Weighting kernel applied to distances.
#' @param sigma Distance scaling for kernel.
#' @param radius Maximum distance away to truncate the kernel.
#' @param normalize Normalize the kernel so that the total sum of weights is equal to this.
#' @keywords layers
#' @export
#' @return A Raster* of the same form as the input.
#' If the factor \code{normalize} is NULL then the result is approximately stochastic,
#' but may be pretty far off if the discretization is very coarse.
#' It is exactly stochastic if \code{normalize} is 1;
#' the interpretation of \code{normalize} more generally is the total production
#' per unit of \code{layer}.
#'
kernel_migrate <- function (layer,
                            radius,
                            fun="gaussian",
                            sigma=1,
                            normalize=NULL
                 ) {
    if (is.character(fun)) {
        fun <- switch( fun,
                gaussian=function (x) {
                        exp(-x^2/2) / (2*pi)
                    },
                cauchy=function (x) {
                        1/(2*pi^2*x*(1+x^2))
                    },
                get(fun,mode="function") 
            )
    }
    area <- prod(raster::res(layer))
    cell.radius <- ceiling(radius/raster::res(layer))
    w <- matrix(nrow=2*cell.radius[1]+1,ncol=2*cell.radius[2]+1)
    cc <- cell.radius+1
    w[] <- fun( sqrt( (xres(layer)*(row(w)-cc[1]))^2 + (yres(layer)*(col(w)-cc[2]))^2 )/sigma ) * area/sigma^2
    if (!is.null(normalize)) { w <- (normalize/sum(w))*w }
    out <- focal( layer, w=w, na.rm=TRUE, pad=TRUE, pad.value=0 )
    out[is.na(layer)] <- NA
    names(out) <- names(layer)
    return(out)
}


##
# Extend \code{focal} to work on Raster* objects.

methods::setMethod("focal", signature("RasterStack"), function(x,...) {
          x.names <- names(x)
          for (k in seq_len(raster::nlayers(x))) {
              x[[k]] <- focal( x[[k]], ... )
          }
          names(x) <- x.names
          return(x)
        } )
methods::setMethod("focal", signature("RasterBrick"), function(x,...) {
          x.names <- names(x)
          for (k in seq_len(raster::nlayers(x))) {
              x[[k]] <- focal( x[[k]], ... )
          }
          names(x) <- x.names
          return(x)
        } )

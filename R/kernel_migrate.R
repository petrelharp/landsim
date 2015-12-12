#' Weighted "Adjacency" Matrix from a Raster Layer.
#'
#' This returns the (sparse) Matrix giving the pseudo-adjacency matrix 
#' for cells in \code{layer}, with weights.
#'
#' @inheritParams layer_adjacency
#' @param fun Weighting kernel applied to distances.
#' @param sigma Distance scaling for kernel.
#' @param radius Maximum distance away to truncate the kernel.
#' @keywords layers
#' @export
#' @return A Raster* of the same form as the input.
kernel_migrate <- function (layer,
                            radius,
                            fun="gaussian",
                            sigma=1,
                            normalize=FALSE
                 ) {
    if (is.character(fun)) {
        fun <- switch( fun,
                gaussian=function (x) {
                        exp(-x^2) / (2*pi)
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
    w[] <- fun( sqrt( (xres(layer)*(row(w)-cc[1]))^2 + (yres(layer)*(col(w)-cc[2]))^2 )/sigma )*area/sigma^2
    if (normalize) { w <- w/sum(w) }
    out <- layer
    for (k in 1:raster::nlayers(layer)) {
        out[[k]] <- focal( layer[[k]], w=w, na.rm=TRUE, pad=TRUE, pad.value=0 )
        out[[k]][is.na(layer[[k]])] <- NA
    }
    return(layer)
}

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
#' @return A sparse Matrix \code{M}, where \code{M[i,j]} is
#'    fun( d[i,j]/\sigma )
#' or, if \code{fun="gaussian"},
#'    \exp( - d[i,j]^2/\sigma^2 ) / ( 2 \pi \sigma^2 ) * (area of a cell)
#' where \code{d[i,j]} is the distance from \code{from[i]} to \code{to[j]},
#' if this is greater than \code{min.prob}, and is 0 otherwise.
#'
#' Approximates the area of each cell to be constant.
kernel_adjacency <- function (layer,
                              radius,
                              fun="gaussian",
                              sigma=1,
                              from=which(raster::values(!is.na(layer))),
                              to=from,
                              normalize=FALSE
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
    directions <- matrix( 1, nrow=2*cell.radius[1]+1, ncol=2*cell.radius[2]+1 )
    directions[cell.radius[1]+1,cell.radius[2]+1] <- 0
    # note columns are (to,from)
    # does NOT include diagonal
    ij <- raster::adjacent(layer,cells=from,target=to,directions=directions,pairs=TRUE,sorted=TRUE)
    pos <- raster::xyFromCell(layer,cell=from)
    ii <- match(ij[,2],from)
    jj <- match(ij[,1],to)
    adj <- Matrix::sparseMatrix( i=ii, j=jj,
            x=fun(raster::pointDistance(pos[ii,],pos[jj,],lonlat=FALSE,allpairs=FALSE)/sigma)*area/sigma^2 )
    if (normalize) {
        adj <- sweep(adj,1,rowSums(adj),"/")
    }
    return(adj)
}


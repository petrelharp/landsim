#' Gaussian-weighted "Adjacency" Matrix from a Raster Layer.
#'
#' This returns the (sparse) Matrix giving the pseudo-adjacency matrix 
#' for cells in \code{layer}, with Gaussian weights.
#'
#' @inheritParams layer_adjacency
#' @param sigma Distance parameter in Gaussian kernel.
#' @param radius Maximum distance away to truncate the kernel.
#' @keywords layers
#' @export
#' @return A sparse Matrix \code{M}, where \code{M[i,j]} is
#'    \exp( - d[i,j]^2/sigma^2 ) / ( 2 \pi \sigma^2 ) * (area of a cell)
#' where \code{d[i,j]} is the distance from \code{from[i]} to \code{to[j]},
#' if this is greater than \code{min.prob}, and is 0 otherwise.
#'
#' Approximates the area of each cell to be constant.
gauss_adjacency <- function (layer,
                             sigma,
                             from=which(values(!is.na(layer))),
                             to=from,
                             radius=20,
                             normalize=FALSE
                 ) {
    area <- prod(res(layer))
    fun <- function (dd) {
        area * exp(-dd^2/sigma^2) / (2*pi*sigma^2)
    }
    cell.radius <- ceiling(radius/res(layer))
    directions <- matrix( 1, nrow=2*cell.radius[1]+1, ncol=2*cell.radius[2]+1 )
    directions[cell.radius[1]+1,cell.radius[2]+1] <- 0
    # note columns are (to,from)
    # does NOT include diagonal
    ij <- raster::adjacent(layer,cells=from,target=to,directions=directions,pairs=TRUE,sorted=TRUE)
    pos <- xyFromCell(layer,cell=from)
    ii <- match(ij[,2],from)
    jj <- match(ij[,1],to)
    adj <- Matrix::sparseMatrix( i=ii, j=jj,
            x=fun(raster::pointDistance(pos[ii,],pos[jj,],lonlat=FALSE,allpairs=FALSE)) )
    if (normalize) {
        adj <- sweep(adj,1,rowSums(adj),"/")
    }
    return(adj)
}


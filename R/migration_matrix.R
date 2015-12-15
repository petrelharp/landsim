#' Weighted "Adjacency" Matrix from a Raster Layer.
#'
#' This returns the (sparse) Matrix giving the pseudo-adjacency matrix 
#' for cells in \code{layer}, with weights.
#'
#' @inheritParams layer_adjacency
#' @param kern Weighting kernel applied to distances.
#' @param sigma Distance scaling for kernel.
#' @param radius Maximum distance away to truncate the kernel.
#' @keywords layers
#' @export
#' @return A sparse Matrix \code{M}, where \code{M[i,j]} is
#'    kern( d[i,j]/\sigma )
#' or, if \code{kern="gaussian"},
#'    \exp( - d[i,j]^2/\sigma^2 ) / ( 2 \pi \sigma^2 ) * (area of a cell)
#' where \code{d[i,j]} is the distance from \code{from[i]} to \code{to[j]},
#' if this is greater than \code{min.prob}, and is 0 otherwise.
#'
#' Approximates the area of each cell to be constant.
migration_matrix <- function (layer,
                              migration=list(sigma=1,normalize=1),
                              kern=migration$kern,
                              sigma=migration$sigma,
                              radius=migration$radius,
                              normalize=migration$normalize,
                              from=which(raster::values(!is.na(layer))),
                              to=from
                 ) {
    kern <- get_kernel(kern)
    area <- prod(raster::res(layer))
    cell.radius <- ceiling(radius/raster::res(layer))
    directions <- matrix( 1, nrow=2*cell.radius[1]+1, ncol=2*cell.radius[2]+1 )
    directions[cell.radius[1]+1,cell.radius[2]+1] <- 0
    # note columns are (to,from)
    # does NOT include diagonal
    ij <- raster::adjacent(layer,cells=from,target=to,directions=directions,pairs=TRUE)
    from.pos <- raster::xyFromCell(layer,cell=from)
    to.pos <- raster::xyFromCell(layer,cell=to)
    ii <- match(ij[,1],from)
    jj <- match(ij[,2],to)
    adj <- Matrix::sparseMatrix( i=ii, j=jj,
            x=kern(raster::pointDistance(from.pos[ii,],to.pos[jj,],lonlat=FALSE,allpairs=FALSE)/sigma)*area/sigma^2 )
    if (!is.null(normalize)) {
        # adj <- (normalize/Matrix::rowSums(adj)) * adj
        # this is twice as quick:
        adj@x <- (normalize/Matrix::rowSums(adj)[1L+adj@i]) * adj@x
    }
    return(adj)
}

# helper function to find column indices of dgCMatrix objects
p.to.j <- function (p) { rep( seq.int( length(p)-1 ), diff(p) ) }

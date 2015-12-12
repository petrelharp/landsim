#' Adjacency Matrix from a Raster Layer.
#'
#' This returns the (sparse) Matrix giving the adjacency matrix of cells in \code{layer}.
#'
#' @param layer A Raster* object.
#' @param from The indices of the cells corresponding to rows in the output matrix. [default: non-NA cells]
#' @param to The indices of the cells corresponding to columns in the output matrix. [default: same as \code{from}]
#' @param directions the number of directions in which cells should be
#'        connected: 4 (rook's case), 8 (queen's case), 16 (knight and
#'        one-cell queen moves), or 'bishop' to connect cells with
#'        one-cell diagonal moves. Or a neigborhood matrix (see \code{raster::adjacent}).
#'        [default: 4]
#' @param normalize If \code{TRUE}, normalize rows to sum to 1.
#' @keywords layers
#' @export
#' @return A sparse Matrix \code{M}, where \code{M[i,j]} is equal to 1 if cell \code{from[i]} is adjacent to \code{to[j]}, and is 0 otherwise.
layer_adjacency <- function (layer,
                             from=which(raster::values(!is.na(layer))),
                             to=from,
                             directions=4,
                             normalize=FALSE
                 ) {
    ij <- raster::adjacent(layer,cells=from,target=to,directions=directions,pairs=TRUE,sorted=TRUE)
    # stopifnot( all(ij[,1] != ij[,2]) ) ## NO DIAGONAL
    # note columns are (to,from)
    adj <- Matrix::sparseMatrix( i=match(ij[,2],from), j=match(ij[,1],to), x=1.0 )
    if (normalize) {
        adj <- sweep(adj,1,rowSums(adj),"/")
    }
    return(adj)
}


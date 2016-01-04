#' Estimate the Root Mean Square Dispersal Distance from a Migration Matrix
#'
#' FOr a migration matrix M[i,j], with rowSus(M) all equal to 1,
#' the root mean square  dispersal distance from location i is 
#'     sqrt( sum( M[i,] * dist[i,]^2 ) );
#' this estimates the average dispersal distance across locations
#' without computing the full distance matrix
#' by taking a random sample of i and averaging.
#'
#' If the rows of M sum to something less than 1 (because migrants may be lost),
#' for the purposes of the migration matrix we should work with
#' the matrix normalized so that rows do sum to 1,
#' since we want migrants conditioned on not being lost.
#'
#' @param M A nonnegative migration matrix, with rows summing to < 1.
#' @param pop A population object containing the raster and sets of accessible and habitable cells corresponding to M.
#' @param normalized Is the matrix M already row-normalized?
#' @export
#' @return A vector of mean dispersal distances, for each cell in pop$habitable.
dispersal_distance <- function ( M,
                                 pop,
                                 normalized=FALSE
                         ) {
    # make sure we know what kind of object we're working with
    M <- as(M,"dgCMatrix")
    # normalize
    if (!normalized) { M@x <- (1/Matrix::rowSums(M)[1L+M@i]) * M@x }
    # get distances between pairs we care about
    cell.pos <- raster::xyFromCell(pop$habitat,cell=which(pop$accessible))
    jj <- p.to.j(M@p)
    M@x <- M@x * raster::pointDistance(cell.pos[1L+M@i,],cell.pos[jj,],lonlat=FALSE,allpairs=FALSE)^2
    return( sqrt(Matrix::rowSums( M )) )
}



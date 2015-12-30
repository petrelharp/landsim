#' Find Points Corresponding to Habitable Cells
#'
#' Given a logical vector corresponding to the habitable locations in a population,
#' return a two-column matrix of coordinates of the TRUE entries in that vector,
#' suitable for plotting.
#'
#' @param pop A population object.
#' @param x A logical vector of length equal to \code{sum(pop$habitable)}.
#' @export
#' @return A two-column (x,y) matrix of coordinates, suitable for plotting or turning into a SpatialPoints object.
habitable_points <- function ( pop, x ) {
    nonz <- which(pop$habitable)[x]
    raster::xyFromCell(pop$habitat,nonz)
}

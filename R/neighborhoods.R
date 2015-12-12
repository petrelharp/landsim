#' Find Indices of Cells Nearby to Specified Points
#'
#' Returns a list of indices of all cells in \code{layer} that are within Euclidean distance \code{dist} of the points in \code{locations}.
#'
#' @param layer A Raster* object.
#' @param dist Radius about each of \code{locations} that cell centers must be within.
#' @param locations Either a SpatialPoints object or a 2-column matrix of coordinates of points to get neighborhoods of.
#' @param mask Indices of cells to return the index of the neighborhoods with respect to. [defaults to \code{1:length(layer)}]
#' @keywords layers
#' @export
#' @return A list of vectors of indices in \code{mask} so that if the output is \code{out},
#'         then \code{mask[out[[k]]]} is a vector of the indices of all cells in \code{layer} whose centers are closer than \code{dist}
#'         to \code{locations[k]}.
get.neighborhoods <- function ( layer, dist, locations, mask ) {
    if ( class(locations)=="SpatialPoints" ) { locations <- sp::coordinates(locations) }  # this is the first thing distanceFromPoints does anyhow
    if (is.null(dim(locations))) { locations <- matrix(locations,ncol=2) }
    neighborhoods <- lapply( 1:NROW(locations) , function (k) {
            dvec <- raster::distanceFromPoints( layer, locations[k,] ) 
            out <- Which( dvec <= max(dist,minValue(dvec)), cells=TRUE, na.rm=TRUE )
            if (!missing(mask)) { out <- match( out, mask ) }
            out[!is.na(out)]
        } )
    return(neighborhoods)
}




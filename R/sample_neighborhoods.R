#' Sample Neighborhoods in a Raster*
#'
#' Chooses a set of randomly sampled points centered on non-NA cells of the input RasterLayer,
#' and returns a randomly chosen maximal subset of nonoverlapping circles centered at these points.
#'
#' @param x The Raster* object.
#' @param n The number of points to sample.
#' @param radius The radius of the neighborhoods.
#' @param nsegs The number of segments in the polygon delimiting each neighborhood.
#' @export
#' @return A named list containing:
#'    centers = a SpatialPoints object of centers
#'    neighborhoods = SpatialPolygons object of the circle(s)
#'    center.cells = the indices of the cells of the Ratser* object that contain the centers
sample_neighborhoods <- function (x, n, radius, nsegs=20, ...) {
    centers <- raster::sampleRandom( x, size=n, xy=TRUE )[,c("x","y")]
    circles <- sp::SpatialPolygons( lapply( 1:nrow(centers), function (k) {
                            cent <- centers[k,]
                            args <- seq(2*pi,0,length.out=nsegs+1)
                            sp::Polygons( list( 
                               sp::Polygon( cbind( x = cent[1] + radius * cos(args),
                                   y = cent[2] + radius * sin(args) ), hole=FALSE )
                               ), ID=paste("neighborhood",k,sep="_") )
                } ) )
    goodones <- rep(TRUE,length(circles))
    for (k in seq_along(circles)[-length(circles)]) {
        checkthese <- goodones & (seq_along(circles)>k)
        if (any(goodones[checkthese])) {
            goodones[checkthese] <- ( goodones[checkthese] & ! rgeos::gIntersects( circles[k], circles[checkthese], byid=TRUE ) )
        }
    }
    center.cells <- raster::extract(habitat,nhoods$centers,cellnumbers=TRUE)[,1]
    return( list( 
                 centers = sp::SpatialPoints(centers[goodones,]), 
                 neighborhoods = circles[goodones],
                 center.cells = center.cells ) )
}

#' Construct Census Functions for a Set of Neighborhoods
#'
#' Given a population object and a set of neighborhoods as returned by \code{sample_neighborhoods},
#' construct a list of functions, the k-th of which, applied to a population object,
#' will return the vector of total numbers of the k-th genotype occurring in each neighborhood.
#'
#' @param pop A population object.
#' @param neighborhoods A list whose element "neighborhoods" is a SpatialPolygons object.
#' @export
#' @return A list of functions, one for each genotype.
census_neighborhoods <- function (pop, neighborhoods) {
    ii <- lapply( raster::extract(pop$habitat,neighborhoods$neighborhoods,cellnumbers=TRUE), 
                 function (x) { z <- match( x[,1], which(pop$habitable) ); z[!is.na(z)] } )
    census.matrix <- Matrix::sparseMatrix( i=rep(seq_len(length(ii)),sapply(ii,length)), j=unlist(ii), dims=c(length(ii),nrow(pop$N)) )
    ff <- lapply( seq_along(pop$genotypes), function (k) {
                     function (pop) { as.vector(census.matrix %*% pop$N[,k]) }
                 } )
    names(ff) <- pop$genotypes
    return(ff)
}

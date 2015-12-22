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
#' @return A SpatialPolygons object containing the circle(s).
sample_neighborhoods <- function (x, n, radius, nsegs=20, ...) {
    centers <- sampleRandom( x, size=n, xy=TRUE )
    circles <- SpatialPolygons( lapply( 1:nrow(centers), function (k) {
                            cent <- centers[k,]
                            args <- seq(2*pi,0,length.out=nsegs+1)
                            Polygons( list( 
                               Polygon( cbind( x = cent[1] + radius * cos(args),
                                   y = cent[2] + radius * sin(args) ), hole=FALSE )
                               ), ID=paste("neighborhood",k,sep="_") )
                } ) )
    goodones <- rep(TRUE,length(circles))
    for (k in seq_along(circles)[-length(circles)]) {
        checkthese <- goodones & (seq_along(circles)>k)
        if (any(goodones[checkthese])) {
            goodones[checkthese] <- ( goodones[checkthese] & ! gIntersects( circles[k], circles[checkthese], byid=TRUE ) )
        }
    }
    return( circles[goodones] )
}


#' Sample Neighborhoods in a Raster*
#'
#' Chooses a set of randomly sampled points centered on non-NA cells of the input RasterLayer,
#' and returns a randomly chosen maximal subset of nonoverlapping circles centered at these points.
#'
#' @param x The Raster* object.
#' @param n The number of points to sample.
#' @param radius The radius of the neighborhoods.
#' @param separation The minimum distance separating neighborhood centers (if \code{thin} is TRUE).
#' @param thin Whether to thin neighborhoods.
#' @export
#' @return A named list containing:
#'    centers = a SpatialPoints object of centers
#'    neighborhoods = SpatialPolygons object of the circle(s)
#'    center.cells = the indices of the cells of the Ratser* object that contain the centers
sample_neighborhoods <- function (x, 
                                  n, 
                                  radius, 
                                  separation=radius, 
                                  thin=TRUE, 
                                  ...) {
    centers <- raster::sampleRandom( x, size=n, xy=TRUE, cells=TRUE )
    goodones <- rep(TRUE,nrow(centers))
    if (thin) {
        circles <- make_circles( centers[,c("x","y"),drop=FALSE], separation, proj4string=CRS(proj4string(x)) )
        for (k in seq_along(circles)[-length(circles)]) {
            checkthese <- goodones & (seq_along(circles)>k)
            if (any(goodones[checkthese])) {
                goodones[checkthese] <- ( goodones[checkthese] & ! rgeos::gIntersects( circles[k], circles[checkthese], byid=TRUE ) )
            }
        }
    }
    return( list( 
                 centers = sp::SpatialPoints(centers[goodones,c("x","y"),drop=FALSE],proj4string=CRS(proj4string(x)) ), 
                 neighborhoods = make_circles( centers[goodones,c("x","y"),drop=FALSE], radius, proj4string=CRS(proj4string(x)) ),
                 center.cells = centers[goodones,"cell"]  ) )
}

#' Make Circles Centered at a Collection of Coordinates
#'
#'
#'
#' @param centers The centers of the circles.
#' @param radii The radii of the circles.
#' @param nsegs The number of segments in the polygon delimiting each neighborhood.
#' @export
#' @return A SpatialPolygons object of circles.
make_circles <- function (centers, radii, proj4string, nsegs=20 ) {
    if (inherits(centers,"SpatialPoints")) { centers <- coordinates(centers) }
    radii <- rep_len(radii,length(centers))
    args <- seq(2*pi,0,length.out=nsegs+1)
    sp::SpatialPolygons( lapply( 1:nrow(centers), function (k) {
                            cent <- centers[k,]
                            sp::Polygons( list( 
                               sp::Polygon( cbind( 
                                          x = cent[1] + radii[k] * cos(args),
                                          y = cent[2] + radii[k] * sin(args) ), hole=FALSE )
                               ), ID=paste("neighborhood",k,sep="_") )
                } ), proj4string=proj4string )
}

#' Construct Census Functions for a Set of Neighborhoods
#'
#' Given a population object and a set of neighborhoods as returned by \code{sample_neighborhoods},
#' construct a list of functions, the k-th of which, applied to a matrix of genotype counts form a population object,
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
                     function (N) { as.vector(census.matrix %*% N[,k]) }
                 } )
    names(ff) <- pop$genotypes
    return(ff)
}

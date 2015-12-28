#' Overlay a Raster with a Sierpinski Mask
#'
#' This constructs a "middle-ninths" Sierpinski gasket of depth n,
#' starting from a given Raster*.
#'
#' @param layer A Raster* object.
#' @param box The bounding box to Cantor-ize.
#' @param n The number of iterations.
#' @param random Choose edges randomly or at (1/3,2/3)?
#' @export
#' @return A Raster* of the same form as the input.
sierpinski_overlay <- function (
                                layer,
                                box=sp::bbox(layer),
                                n=floor(log(max(nrow(layer),ncol(layer))/log(3))),
                                random=FALSE
                        ) {
    if (random) {
        mid.s1 <- sort( c( box[1,], runif(8,min=box[1,1],max=box[1,2]) ) )[c(1,4,7,10)]
        mid.s2 <- sort( c( box[2,], runif(8,min=box[2,1],max=box[2,2]) ) )[c(1,4,7,10)]
    } else {
        mid.s1 <- sort( c( box[1,], box[1,1]+c(1,2)*diff(box[1,])/3 ) )
        mid.s2 <- sort( c( box[2,], box[2,1]+c(1,2)*diff(box[2,])/3 ) )
    }
    # remove the middle ninth
    cutbox.coords <- cbind( s1=mid.s1[1+c(1,2,2,1,1)],
                            s2=mid.s2[1+c(1,1,2,2,1)] )
    cutpoly <- sp::Polygon( coords=cutbox.coords )
    cutpoly.sp <- sp::SpatialPolygons( list( sp::Polygons( list( cutpoly ), ID="cutbox" ) ) )
    cutmask <- raster::rasterize( cutpoly.sp, layer, background=NA )
    layer <- raster::mask( layer, cutmask, inverse=TRUE )
    # iterate 
    if (n>1) {
        for (ii in 1:3) {
            for (jj in 1:3) {
                if ( (ii!=2) || (jj!=2) ) {
                    thisbox <- rbind( mid.s1[c(ii,ii+1)],
                                      mid.s2[c(jj,jj+1)] )
                    layer <- sierpinski_overlay( layer=layer, box=thisbox, n=n-1, random=random )
                }
            }
        }
    }
    return(layer)
}



#' Simulate a Random Habitat
#'
#' Simulates a random habitat, with holes and varying carrying capacity.
#'
#' @param diam Diameter of the desired habitat (in real units, not number of cells).
#' @param res Spatial resolution of the Raster.
#' @param randfun Function used to sample values for the habitat (negative values will be set to NA).
#' @param kern Smoothing kernel.
#' @param sigma Smoothing kernel scale.
#' @param radius Smoothing kernel maximum radius.
#' @param width Width of the region (default: diam).
#' @param height Height of the region (default: diam).
#' @export
#' @return A RasterLayer with nonnegative and missing values.
random_habitat <- function ( diam=2e4, 
                             res=100, 
                             randfun=function(n)pmin(20,(2+rcauchy(n))),
                             kern="gaussian", 
                             sigma=300, 
                             radius=1500,
                             width=diam,
                             height=diam
                            ) {

    habitat <- raster::raster(
          xmn=-width/2, xmx=width/2, ymn=-height/2, ymx=height/2, 
          resolution=res,
          crs="+proj=utm +zone=11 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
    raster::values(habitat) <- randfun(raster::ncell(habitat))
    habitat <- 20*migrate_raster( habitat, kern=kern, sigma=sigma, radius=radius )
    raster::values(habitat)[raster::values(habitat)<0] <- NA
    return(habitat)
}

#' Create a Flat Habitat
#'
#' Simulates a random environment, with holes and varying carrying capacity.
#'
#' @param diam Diameter of the desired habitat (in real units, not number of cells).
#' @param res Spatial resolution of the Raster.
#' @param value The value(s) to populate the habitat with (will be recycled).
#' @param width Width of the region (default: diam).
#' @param height Height of the region (default: diam).
#' @export
#' @return A RasterLayer with nonnegative and missing values.
flat_habitat <- function ( diam=2e4, 
                           res=100, 
                           value=1,
                           width=diam,
                           height=diam
                         ) {
    habitat <- raster::raster(
          xmn=-width/2, xmx=width/2, ymn=-height/2, ymx=height/2, 
          resolution=res,
          crs="+proj=utm +zone=11 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
    raster::values(habitat) <- value
    return(habitat)
}

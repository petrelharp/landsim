#' Simulate a Random Environment
#'
#' Simulates a random environment, with holes and varying carrying capacity.
#'
#' @param diam Diameter of the desired habitat (in real units, not number of cells).
#' @param res Spatial resolution of the Raster.
#' @param randfun Function used to sample values for the habitat (negative values will be set to NA).
#' @param kern Smoothing kernel.
#' @param sigma Smoothing kernel scale.
#' @param radius Smoothing kernel maximum radius.
#' @export
#' @return A RasterLayer with nonnegative and missing values.
random_habitat <- function ( diam=1e4, 
                             res=100, 
                             randfun=function(n)pmin(20,(2+rcauchy(n))),
                             kern="gaussian", 
                             sigma=300, 
                             radius=1500
                            ) {

    habitat <- raster::raster(xmn=-diam, xmx=diam, ymn=-diam, ymx=diam, 
          resolution=res,
          crs="+proj=utm +zone=11 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
    raster::values(habitat) <- randfun(raster::ncell(habitat))
    habitat <- 20*migrate_raster( habitat, kern=kern, sigma=sigma, radius=radius )
    raster::values(habitat)[raster::values(habitat)<0] <- NA
    return(habitat)
}

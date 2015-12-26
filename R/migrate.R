#' Migrate the Quantities in a Population.
#'
#' This applies a \code{migration} operation to a given population object.
#'
#' @param migration The \code{migration} object.
#' @param population The population object to apply the migration operation to.
#' @param x Values to be migrated, if not those given in \code{population}.
#' @export
#' @return A matrix of values  of the same form as \code{x}.  See \code{migrate} for details.
migrate <- function ( migration, 
                      population,
                      x=population$N
                 ) {
    out <- as( migration$M %*% x, class(x) )  # return a dense matrix if we got one in
    return(out)
}


#' Migrate the Quantities in a Raster.
#'
#' This applies a \code{migration} operation to a given Raster* object.
#'
#' @param migration The \code{migration} object containing the relevant parameters, which may alternatively be passed in individually.
#' @param population The population object to apply the migration operation to.
#' @param kern Weighting kernel applied to distances.
#' @param sigma Distance scaling for kernel.
#' @param radius Maximum distance away to truncate the kernel.
#' @param normalize Normalize the kernel so that the total sum of weights is equal to this; pass NULL to do no normalization.
#' @export
#' @return A Raster* of the same form as the input.
#' If the factor \code{normalize} is NULL then the result is approximately stochastic,
#' but may be pretty far off if the discretization is very coarse.
#' It is exactly stochastic if \code{normalize} is 1;
#' the interpretation of \code{normalize} more generally is the total production
#' per unit of \code{population}.
#'
#' However, note that even if \code{normalize} is 1, the migration will still not be conservative
#' at any raster cells nearby to boundary or NA cells.
migrate_raster <- function (population,
                            migration=list(sigma=1,normalize=1),
                            kern=migration$kern,
                            sigma=migration$sigma,
                            radius=migration$radius,
                            normalize=migration$normalize
                 ) {
    if (inherits(population,"population")) { population <- population$habitat }
    if (!inherits(population,"Raster")) {
        stop("migrate_raster: population must be a population or a Raster* object.")
    } 
    kern <- get_kernel(kern)
    area <- prod(raster::res(population))
    cell.radius <- ceiling(radius/raster::res(population))
        w <- matrix(nrow=2*cell.radius[1]+1,ncol=2*cell.radius[2]+1)
        cc <- cell.radius+1
        w[] <- kern( sqrt( (xres(population)*(row(w)-cc[1]))^2 + (yres(population)*(col(w)-cc[2]))^2 )/sigma ) * area/sigma^2
        if (!is.null(normalize)) { w <- (normalize/sum(w))*w }
        out <- focal( population, w=w, na.rm=TRUE, pad=TRUE, padValue=0 )
        out[is.na(population)] <- NA
    names(out) <- names(population)
    return(out)
}

# helper function used elsewhere as well
get_kernel <- function (kern) {
    if (is.character(kern)) {
        kern <- switch( kern,
                gaussian=function (x) {
                        exp(-x^2/2) / (2*pi)
                    },
                cauchy=function (x) {
                        1/(2*pi^2*x*(1+x^2))
                    },
                get(kern,mode="function") 
            )
    } else { kern }
}

##
# Extend \code{focal} to work on Raster* objects.

require(raster)

methods::setMethod("focal", signature("RasterStack"), function(x,...) {
          x.names <- names(x)
          for (k in seq_len(raster::nlayers(x))) {
              x[[k]] <- focal( x[[k]], ... )
          }
          names(x) <- x.names
          return(x)
        } )
methods::setMethod("focal", signature("RasterBrick"), function(x,...) {
          x.names <- names(x)
          for (k in seq_len(raster::nlayers(x))) {
              x[[k]] <- focal( x[[k]], ... )
          }
          names(x) <- x.names
          return(x)
        } )

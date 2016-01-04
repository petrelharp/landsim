#' Migrate the Quantities in a Population.
#'
#' This applies a \code{migration} operation to a given population object.
#'
#' @param x A numeric matrix, or a population object, in which case migration will be applied to \code{x$N}.
#' @param migration The \code{migration} object.
#' @export
#' @return A matrix of values  of the same form as \code{x}.  See \code{migrate} for details.
migrate <- function ( x,
                      migration
                 ) {
    if (is.null(migration$M)) { stop("Migration matrix does not exist: need to use setup_demography, or migrate_raster()?") }
    if (inherits(x,"population")) { x <- x$N }
    orig.class <- class(x)
    # set up for weighted sum
    Mnx <- matrix( 0.0, nrow=nrow(migration$M), ncol=NCOL(x) )
    Mnx[migration$habitable.inds,] <- x
    zero.weight <- 1-sum(migration$n.weights)
    if (zero.weight<0) { warn("n.weights sum to more than one.") }
    out <- if (zero.weight>0) {
            zero.weight*Mnx[migration$habitable.inds,] 
        } else { 
            numeric(length(migration$habitable.inds)) 
        }
    for (n in seq_along(migration$n.weights)) {
        Mnx <- migration$M %*% Mnx
        if (migration$n.weights[n]>0) {
            out <- out + migration$n.weights[n]*Mnx[migration$habitable.inds,]
        }
    }
    # return a dense matrix if we got one in
    return( as(out,orig.class) )
}


#' Migrate the Quantities in a Raster.
#'
#' This applies a \code{migration} operation to a given Raster* object.
#'
#' @param migration The \code{migration} object containing the relevant parameters, which may alternatively be passed in individually.
#' @param x The Raster* to apply the migration operation to.
#' @param kern Weighting kernel applied to distances.
#' @param sigma Distance scaling for kernel.
#' @param radius Maximum distance away to truncate the kernel.
#' @param normalize Normalize the kernel so that the total sum of weights is equal to this; pass NULL to do no normalization.
#' @param n.weights The resulting operator is x -> (1-sum(n.weights)) * x + n.weights[1] * M x + n.weights[2] * M^2 x + ...
#' @export
#' @return A Raster* of the same form as the input.
#' If the factor \code{normalize} is NULL then the result is approximately stochastic,
#' but may be pretty far off if the discretization is very coarse.
#' It is exactly stochastic if \code{normalize} is 1;
#' the interpretation of \code{normalize} more generally is the total production
#' per unit of \code{x}.
#'
#' However, note that even if \code{normalize} is 1, the migration will still not be conservative
#' at any raster cells nearby to boundary or NA cells.
#'
#' The weights in \code{n.weights} should sum to 1.
migrate_raster <- function (x,
                            migration=list(sigma=1,normalize=1,n.weights=1),
                            kern=migration$kern,
                            sigma=migration$sigma,
                            radius=migration$radius,
                            normalize=migration$normalize,
                            n.weights=migration$n.weights
                 ) {
    if (inherits(x,"population")) { x <- x$habitat }
    if (!inherits(x,"Raster")) {
        stop("migrate_raster: x must be a population or a Raster* object.")
    } 
    kern <- get_kernel(kern)
    area <- prod(raster::res(x))
    cell.radius <- ceiling(radius/raster::res(x))
    w <- matrix(nrow=2*cell.radius[1]+1,ncol=2*cell.radius[2]+1)
    cc <- cell.radius+1
    w[] <- kern( sqrt( (raster::xres(x)*(row(w)-cc[1]))^2 + (raster::yres(x)*(col(w)-cc[2]))^2 )/sigma ) * area/sigma^2
    if (!is.null(normalize)) { w <- (normalize/sum(w))*w }
    x.na <- is.na(x)
    x.names <- names(x)
    # weighted sum
    out <- (1-sum(n.weights))*x
    for (n in seq_along(n.weights)) {
        x <- raster::focal( x, w=w, na.rm=TRUE, pad=TRUE, padValue=0 )
        if (n.weights[n]>0) {
            out <- out + n.weights[n]*x
        }
    }
    x[x.na] <- NA
    names(x) <- x.names
    return(x)
}

# helper function used elsewhere as well
get_kernel <- function (kern) {
    if (is.character(kern)) {
        kern <- switch( kern,
                gaussian=function (x) {
                        exp(-x^2/2) / (2*pi)
                    },
                cauchy=function (x) {
                        1/(pi^2*(1+x^2))
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
              x[[k]] <- raster::focal( x[[k]], ... )
          }
          names(x) <- x.names
          return(x)
        } )
methods::setMethod("focal", signature("RasterBrick"), function(x,...) {
          x.names <- names(x)
          for (k in seq_len(raster::nlayers(x))) {
              x[[k]] <- raster::focal( x[[k]], ... )
          }
          names(x) <- x.names
          return(x)
        } )

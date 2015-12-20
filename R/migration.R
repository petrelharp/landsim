#' Construct a \code{migration} Object.
#'
#' This puts together the parameters necessary to specify a \code{migration} object,
#' specifying the form of the dispersal and the total number of dispersers.
#'
#' @param kern A function of scaled distance giving relative weights. (or, "gaussian", "cauchy")
#' @param sigma A number that distances are scaled by before being passed to \code{kern}.
#' @param radius The maximum distance to migrate over.
#' @param normalize Total number of migrants per unit of input.
#' @param do.M Precompute an explicit migration matrix?
#' @param population A Raster* or a \code{population}; used in computing \code{M}.
#' @param from Used in computing \code{M}.
#' @param to Used in computing \code{M}.
#' @param ... Other parameters that are included verbatim in the output object.
#' @export
#' @return A \code{migration} S3 object (just a named list).
#'
#' @section{Details}
#'
#' If \code{kern} is a migration object, can be used to add a migration matrix to it.
#' See \code{migration_matrix}.
migration <- function (
                       kern, 
                       radius, 
                       sigma=1, 
                       normalize=1, 
                       do.M=FALSE,
                       population, 
                       from=which(raster::values(!is.na(population))),
                       to=from,
                       ...)  {
    if (inherits(kern,"migration")) {
        out <- kern
    } else {
        out <- c( list(
                    kern=kern,
                    radius=radius,
                    sigma=sigma,
                    normalize=normalize ), 
                 list(...) )
    }
    class(out) <- "migration"
    if (do.M) {
        if (missing(population)) { stop("To compute M you need a population or a Raster.") }
        out$M <- migration_matrix( population=population, migration=out, from=from, to=to )
        class(out) <- c(class(out),"migration.matrix")
    }
    return(out)
}


#' Construct a \code{migration} Object.
#'
#' This puts together the parameters necessary to specify a \code{migration} object,
#' specifying the form of the dispersal and the total number of dispersers.
#'
#' @param kern A function of scaled distance giving relative weights. (or, "gaussian", "cauchy")
#' @param sigma A number that distances are scaled by before being passed to \code{kern}.
#' @param radius The maximum distance to migrate over.
#' @param normalize Total number of migrants per unit of input.
#' @param n The number of times to apply the smoother per migration step.
#' @param do.M Precompute an explicit migration matrix?
#' @param population A Raster* or a \code{population}; used in computing \code{M}.
#' @param from Used in computing \code{M}.
#' @param to Used in computing \code{M}.
#' @param ... Other parameters that are included verbatim in the output object.
#' @export
#' @return A \code{migration} S3 object (just a named list).
#'
#' If \code{kern} is a migration object, can be used to add a migration matrix to it.
#' See \code{migration_matrix}.
#'
#' Since the memory usage of the migration matrix increases quadratically with the size of the radius,
#' it can be useful to keep it small and instead apply the migration operator several times;
#' this is achieved by setting \code{n} to something larger than 1.
migration <- function (
                       kern, 
                       radius, 
                       sigma=1, 
                       normalize=1, 
                       n=1,
                       do.M=FALSE,
                       population, 
                       from,
                       to,
                       ...)  {
    if (inherits(kern,"migration")) {
        out <- kern
    } else {
        out <- c( list(
                    kern=kern,
                    radius=radius,
                    sigma=sigma,
                    normalize=normalize,
                    n=n ), 
                 list(...) )
        class(out) <- "migration"
    }
    if (do.M) {
        if (missing(population)) { stop("To compute M you need a population or a Raster.") }
        mm.args <- list(population=population, migration=out)
        if (!missing(from)) { mm.args <- c( mm.args, list( from=from ) ) }
        if (!missing(to)) { mm.args <- c( mm.args, list( to=to ) ) }
        out$M <- do.call( migration_matrix, mm.args )
        class(out) <- c(class(out),"migration.matrix")
    }
    return(out)
}


#' Set up Migration Object for a Concrete Population
#'
#' Migration objects are agnostic about the map they work on,
#' but for application we need to add a migration matrix.
#'
#' @param m The migration object.
#' @param population The population object.
#' @param ... Other parameters passed to \code{migration()}.
#' @export
#' @return A \code{migration} object, as before, but with but with a migration matrix (see \code{migration()}).
#'
setup_migration <- function (m,population,...) {
    m <- migration(m,population=population,do.M=TRUE,...)
}

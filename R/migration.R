#' Construct a \code{migration} Object.
#'
#' This puts together the parameters necessary to specify a \code{migration} object,
#' specifying the form of the dispersal and the total number of dispersers.
#'
#' @param kern A function of scaled distance giving relative weights. (or, "gaussian", "cauchy")
#' @param sigma A number that distances are scaled by before being passed to \code{kern}.
#' @param radius The maximum distance to migrate over.
#' @param normalize Total number of migrants per unit of input.
#' @param ... Other parameters that are included verbatim in the output object.
#' @export
#' @return A \code{migration} S3 object (just a named list).
#'
migration <- function (
                       kern, 
                       radius, 
                       sigma=1, 
                       normalize=1, 
                       ...)  {
    out <- c( list(
                kern=kern,
                radius=radius,
                sigma=sigma,
                normalize=normalize ), 
             list(...) )
    class(out) <- "migration"
    return(out)
}


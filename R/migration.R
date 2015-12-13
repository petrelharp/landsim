#' Construct a \code{migration} Object.
#'
#' This puts together the parameters necessary to specify a \code{migration} object,
#' specifying the form of the dispersal and the total number of dispersers.
#'
#' @param kernel A function of scaled distance giving relative weights. (or, "gaussian", "cauchy")
#' @param sigma A number that distances are scaled by before being passed to \code{kernel}.
#' @param radius The maximum distance to migrate over.
#' @param normalize Total number of migrants per unit of input.
#' @param ... Other parameters that are included verbatim in the output object.
#' @export
#' @return A \code{migration} S3 object (just a named list).
#'
migration <- function (
                       kernel, 
                       radius, 
                       sigma=1, 
                       normalize=1, 
                       ...)  {
    if (is.character(kernel)) {
        kernel <- switch( kernel,
                gaussian=function (x) {
                        exp(-x^2/2) / (2*pi)
                    },
                cauchy=function (x) {
                        1/(2*pi^2*x*(1+x^2))
                    },
                get(kernel,mode="function") 
            )
    }
    out <- c( list(
                kernel=kernel,
                radius=radius,
                sigma=sigma,
                normalize=normalize ), 
             list(...) )
    class(out) <- "migration"
    return(out)
}


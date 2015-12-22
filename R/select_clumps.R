
#' Divide a Raster into Regions Separated by a Minimum Distance
#'
#' Classifies the non-NA values in a Raster* object into the unqiue minimal sets
#' such that if two cells are within a given threshold distance of each other,
#' they are in the same set.
#'
#' @param x The Raster* object.
#' @param threshold The minimum distance defining the groups.
#' @param ... Additional parameters passed to \code{raster::clump}.
#' @export
#' @return A Raster* of the same form as \code{x}.
select_clumps <- function (x, threshold, ...) {
    cats <-  raster::clump( raster::distance(x) <= threshold, ... )
    catord <- order( tabulate(values(cats)), decreasing=TRUE )
    values(x)[!is.na(values(x))] <- catord[values(cats)[!is.na(values(x))]]
    return(x)
}


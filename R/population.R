#' Construct a \code{population} Object.
#'
#' This puts together the parameters necessary to specify a \code{population} object,
#' specifying the reference map (as a Raster), accessible and habitable locations,
#' the genotypes, and the census at each habitable location.
#'
#' @param habitat RasterLayer with values that can be used in other computations.
#' @param accessible Vector of indices of cells in `habitat` that migrants will attempt to go to.
#' @param habitable Vector of indices of cells in `habitat` that may have positive population.
#' @param genotypes Character vector of the genotypes.
#' @param N Matrix indexed by (habitable cells) x (genotypes) giving the number of each genotype in each habitable cell.
#' @param ... Other parameters that are included verbatim in the output object.
#' @export
#' @return A \code{population} S3 object (just a named list).
#'
population <- function (
                        habitat,
                        accessible=1:ncells(habitat),
                        habitable=1:ncells(habitat),
                        genotypes=colNames(N),
                        N=matrix(1,nrow=length(habitable),ncol=length(genotypes)),
                       ...)  {
    out <- c( list(
                habitat=habitat,
                accessible=accessible,
                habitable=habitable,
                genotypes=genotypes,
                N=N,
            ), list(...) )
    class(out) <- "population"
    return(out)
}

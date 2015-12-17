#' Construct a \code{population} Object.
#'
#' This puts together the parameters necessary to specify a \code{population} object,
#' specifying the reference map (as a Raster), accessible and habitable locations,
#' the genotypes, and the census at each habitable location.
#'
#' @param habitat RasterLayer with values that can be used in other computations.
#' @param accessible Logical vector indicating which cells in `habitat` that migrants will attempt to go to.
#' @param habitable Locgical vector indicating which cells in `habitat` that may have positive population.
#' @param genotypes Character vector of the genotypes.
#' @param N Matrix indexed by (habitable cells) x (genotypes) giving the number of each genotype in each habitable cell.
#' @param ... Other parameters that are included verbatim in the output object.
#' @export
#' @return A \code{population} S3 object (just a named list).
population <- function (
                        habitat,
                        accessible=!is.na(values(habitat)),
                        habitable=!is.na(values(habitat)),
                        genotypes=colNames(N),
                        N=matrix(1,nrow=sum(habitable),ncol=length(genotypes)),
                       ...)  {
    if ( (ncol(N) != length(genotypes)) || (nrow(genotypes)!=sum(habitable)) ) {
        stop("N must be (number of habitable cells) x (number of genotypes)")
    }
    colnames(N) <- genotypes
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


#' Convert a population to a RasterBrick.
#'
#' Puts a set of values into the map given by the habitat of the population.
#'
#' @param population The population object.
#' @param x The values to be put into the Raster.
pop_to_raster <- function (population,x=population$N) {
    out <- stack( pop['habitat'][rep(1,NCOL(x))] )
    names(out) <- pop$genotypes
    values(out)[pop$habitable] <- x
    return(out)
}

require(raster)
setAs(from="population",to="Raster",def=pop_to_raster)

plot.population <- function (pop,...) { plot(as(pop,"Raster"),...) }

#' Number of Habitable Cells in a Population.
#'
#' Returns the number of habitable cells in a population.
#'
#' @param population The population object.
#' @export
nhabitable <- function (population) {
    return( nrow(population$N) )
}

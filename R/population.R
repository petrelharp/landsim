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
#' @param description A description of this object.
#' @param ... Other parameters that are included verbatim in the output object.
#' @export
#' @return A \code{population} S3 object (just a named list).
population <- function (
                        habitat,
                        accessible=!is.na(values(habitat)),
                        habitable=!is.na(values(habitat)),
                        genotypes=colNames(N),
                        N=matrix(1,nrow=sum(habitable),ncol=length(genotypes)),
                        description='',
                       ...)  {
    if ( (NCOL(N) != length(genotypes)) || (NROW(N)!=sum(habitable)) ) {
        stop("N must be (number of habitable cells) x (number of genotypes)")
    }
    if ( !is.logical(accessible) || (length(accessible) != ncell(habitat)) 
        || !is.logical(accessible) || (length(accessible) != ncell(habitat)) ) {
        stop("accessible and habitable must be logical vectors of the same length as the number of cells in habitat.")
    }
    colnames(N) <- genotypes
    out <- c( list(
                habitat=habitat,
                accessible=accessible,
                habitable=habitable,
                genotypes=genotypes,
                N=N,
                description=description
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
    values(out)[pop$habitable] <- as.numeric(x)
    return(out)
}

# this only works for S4 methods
# require(raster)
# setAs(from="population",to="Raster",def=function(from)pop_to_raster(population=from))

plot.population <- function (pop,zlim=range(pop$N,finite=TRUE),...) { plot(pop_to_raster(pop),zlim=zlim,...) }
values.population <- function (pop) { pop$N }

#' Number of Habitable Cells in a Population.
#'
#' Returns the number of habitable cells in a population.
#'
#' @param population The population object.
#' @export
nhabitable <- function (population) {
    return( nrow(population$N) )
}


#' Access Values of Habitable Cells in a Population
#'
#' Provides a method to read values to the matrix of population sizes
#' in a population object
#' indexed by cells in the underlying Raster* rather than rows the matrix itself,
#' which records only values at habitable locations.
#'
#' @param x The population object.
#' @param i The rows to read values of, indexed by cells in x$habitat.
#' @param j The columns to read values of, or names of genotypes.
#' @export
get_N <- function (x,i,j) {
    if (is.character(j)) { j <- match(j, x$genotypes) }
    x$N[ cbind(match(i,which(x$habitable)), j) ]
}

#' Assign Values to Habitable Cells in a Population
#'
#' Provides a method to assign values to the matrix of population sizes
#' in a population object
#' indexed by cells in the underlying Raster* rather than rows the matrix itself,
#' which records only values at habitable locations.
#'
#' @param x The population object.
#' @param value The new values.
#' @param i The rows to replace values of, indexed by cells in x$habitat.
#' @param j The columns to replace values of, or names of genotypes.
#' @param ... Arguments passed to "[.matrix"
#' @export
set_N <- function (x,i,j,...,value) {
    if (is.character(j)) { j <- match(j, x$genotypes) }
    x$N[ cbind(match(i,which(x$habitable)), j), ... ] <- value
    return(x)
}



#' Set Up a Population from a Configuration
#'
#' Provides a method to read and store population configuration.
#'
#' @param habitat A RasterLayer, or the path to the file where the raster is stored.
#' @param inaccessible.value The values in the raster that should be marked as inaccessible.
#' @param uninhabitable.value The values in the raster that should be marked as not habitable.
#' @param genotypes A character vector of genotypes.
#' @param N Numerical value(s) to initialize the matrix of genotype counts with.
#' @param extent A raster::extent object (or a vector of xmin,xmax,ymin,ymax) to crop the habitat with.
#' @export
make_population <- function (
                             habitat,
                             inaccessible.value,
                             uninhabitable.value,
                             genotypes,
                             N=0,
                             extent
                         ) {
    if (!inherits(habitat,"Raster")) { habitat <- raster::raster(habitat) }
    if (!missing(extent)) { habitat <- raster::crop(habitat,extent) }
    accessible <- if (is.na(inaccessible.value)) { 
            is.na(raster::values(habitat)) 
        } else { 
            !is.na(values(habitat)) & (raster::values(habitat) != inaccessible.value) 
        }
    habitable <- if (is.na(uninhabitable.value)) { 
            is.na(raster::values(habitat)) 
        } else { 
            !is.na(values(habitat)) & (raster::values(habitat) != uninhabitable.value) 
        }
    population(
              habitat = habitat,
              accessible = accessible, 
              habitable = habitable,
              genotypes = genotypes,
              N = matrix(N,nrow=sum(habitable),ncol=length(genotypes))
          )
}


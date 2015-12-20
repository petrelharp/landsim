#' Simulate Many Generations
#'
#' This steps a population forwards by many generations,
#' according to the parameters in a \code{demography} object,
#' returning the result at specified times in an array.
#'
#' @inheritParams generation
#' @param population A \code{population} object, the initial state (t=0) of the population.
#' @param demography A \code{demography} object, containing the below parameters.
#' @param times A sorted integer vector of times to record the full state of the population at.
#' @param summary.times A sorted integer vector of times to record summaries of the population at.
#' @param summaries A list of functions to apply to the population at generations \code{summary.times}.
#' @param ... Additional parameters that will be passed to \code{generation()}.
#' @export
#' @return An array with dimensions as \code{population$N} x \code{length(times)}.
simulate <- function (
                        population,
                        demography,
                        times,
                        summary.times=1:max(times),
                        summaries=NULL,
                        ...
                ) {
    out <- array(
               dim=c(dim(population$N),length(times)),
               dimnames=list( rownames( population$N ), colnames( population$N ), times )
           )
    summode <- lapply( summaries, function (f) { f(population) } )
    sums <- lapply( summode, function (x) { 
                   s.out <- vector( mode=mode(x), length=length(summary.times)*length(as.vector(x)) )
                   dim(s.out) <- c( length(summary.times), length(as.vector(x)) )
                   colnames(s.out) <- names(x)
                   s.out
           } )
    t <- 0   # the current time
    k <- 1   # index of the last time to record at 
    ks <- 1  # index of the next time to record summaries at 
    stopifnot(is.finite(max(times,summary.times)))
    while ( (k <= length(times)) || ( (length(summaries)>0) && (ks <= length(summary.times)) ) ) {
        if ( (k <= length(times)) && (t>=times[k]) ) {
            out[,,k] <- as.numeric(population$N)
            k <- k+1
        }
        if ( (length(summaries)>0) && (ks <= length(summary.times)) && (t>summary.times[ks]) ) {
            for (j in seq_along(summaries)) {
                sums[[j]][ks,] <- as.vector( summaries[[j]](population) )
                ks <- ks+1
            }
        }
        population$N <- generation(population,demography,t=t,...)
        t <- t+1
    }
    return(list(N=out,summaries=sums,times=times,summary.times=summary.times,final.t=t))
}

#' Convert Simulation Array to RasterBrick(s).
#'
#' From a population object and an array of values 
#' whose dimensions correspond to habitable cells, genotypes, and times (in that order),
#' create a list of RasterBricks, one for each genotype.
#'
#' @param sim An array, as above.
#' @param pop A \code{population} object.
#' @export
#' @return A named list of RasterBrick objects.
sim_to_brick <- function (sim,pop) {
    out <- lapply( pop$genotypes, function (geno) {
                  rb <- do.call( stack, list( pop$habitat )[rep(1,dim(sim$N)[3])] )
                  values(rb)[pop$habitable] <- sim$N[,match(geno,pop$genotypes),]
                  names(rb) <- sim$times
                  rb
           } )
    names(out) <- pop$genotypes
    return(out)
}



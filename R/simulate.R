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
#' @param N Overrides the matrix of genotype counts in \code{population}, if present.
#' @param tinit The current time at the beginning of the simulation (defaults to t=0).
#' @param summary.times A sorted integer vector of times to record summaries of the population at.
#' @param summaries A list of functions to apply to the matrix of genotype counts, population$N, at generations \code{summary.times}.
#' @param ... Additional parameters that will be passed to \code{generation()}.
#' @export
#' @return A named list, with elements
#'       N = an array with dimensions as \code{population$N} x \code{length(times)}.
#'       summaries = A list of matrices, with dimensions (length of summary.times) x (length of output of summary function).
#'       times = As in the input.
#'       summary.times = As in the input.
#'       t = The final generation.
simulate <- function (
                        population,
                        demography,
                        times,
                        N=population$N,
                        tinit=0,
                        summary.times=1:max(times),
                        summaries=NULL,
                        ...
                ) {
    # Note we update N below, not population.
    out <- array(
               dim=c(dim(N),length(times)),
               dimnames=list( rownames( N ), colnames( N ), times )
           )
    summode <- lapply( summaries, function (f) { f(N) } )
    sums <- lapply( summode, function (x) { 
                   s.out <- vector( mode=mode(x), length=length(summary.times)*length(as.vector(x)) )
                   dim(s.out) <- c( length(summary.times), length(as.vector(x)) )
                   colnames(s.out) <- names(x)
                   s.out
           } )
    t <- tinit  # current time
    k <- 1   # index of the last time to record at 
    ks <- 1  # index of the next time to record summaries at 
    stopifnot(is.finite(max(times,summary.times)))
    while ( (k <= length(times)) || ( (length(summaries)>0) && (ks <= length(summary.times)) ) ) {
        if ( (k <= length(times)) && (t>=times[k]) ) {
            # record at first time point greater than or equal to entries in times
            out[,,k] <- as.numeric(N)
            k <- k+1
        }
        if ( (length(summaries)>0) && (ks <= length(summary.times)) && (t>=summary.times[ks]) ) {
            for (j in seq_along(summaries)) {
                sums[[j]][ks,] <- as.vector( summaries[[j]](N) )
            }
            ks <- ks+1
        }
        # don't do the last, unnecessary step...
        if (t<max(times,summary.times)) {
            # generation() really only needs N, not the whole population; but pass it in anyhow
            N <- generation(population,demography,N=N,t=t,...)
            t <- t+1
        }
    }
    outlist <- list(N=out,summaries=sums,times=times,summary.times=summary.times,t=t)
    class(outlist) <- "popsim"
    return(outlist)
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



#' Continue a Previously Initiated Simulation
#'
#' Run a previously initiated popluation for more generations.
#'
#' @inheritParams simulate
#' @param sim The previous \code{simulation} object.
#' @param append If TRUE, append results to results in previous simulation object.
#' @export
#' @return A simulation object, as \code{simulation}.
extend_simulation <- function (
                        sim,
                        population,
                        demography,
                        times,
                        append=TRUE,
                        N=sim$N[,,dim(sim$N)[3]],
                        tinit=sim$t+1,
                        summary.times=(sim$t+1):max(times),
                        summaries=NULL,
                        ...
                ) {
    new.sim <- simulate( population=population, demography=demography, times=times, N=N, tinit=tinit, summary.times=summary.times, summaries=summaries, ... )
    if (append) {
        new.dims <- dim(new.sim$N)
        new.sim$N <- c( sim$N, new.sim$N )
        dim(new.sim$N) <- c( 0, 0, new.dims[3] ) + dim( sim$N )
        for (k in seq_along(new.sim$summaries)) {
            if (ncol(new.sim$summaries[[k]])!=ncol(sim$summaries[[k]])) { stop("summaries must be the same to extend a simulation object") }
            new.sim$summaries[[k]] <- rbind( sim$summaries[[k]], new.sim$summaries[[k]] )
        }
        new.sim$times <- c(sim$times,new.sim$times)
        new.sim$summary.times <- c(sim$summary.times,new.sim$summary.times)
    }
    return(new.sim)
}



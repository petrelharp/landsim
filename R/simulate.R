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
#' @param stop.fun A function applied to N that, if it returns TRUE, will stop the simulation. (remainder of output array will be entirely zeros)
#' @param ... Additional parameters that will be passed to \code{generation()}.
#' @export
#' @return A named list, with elements
#'       N = an array with dimensions as \code{population$N} x \code{length(times)}.
#'       summaries = A list of matrices, with dimensions (length of summary.times) x (length of output of summary function).
#'       times = As in the input.
#'       summary.times = As in the input.
#'       t = The final generation.
simulate_pop <- function (
                        population,
                        demography,
                        times,
                        N=population$N,
                        tinit=0,
                        summary.times=1:max(times),
                        summaries=NULL,
                        stop.fun=NULL,
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
        if (!is.null(stop.fun) && stop.fun(N)) {
            cat("stopping early, at generation ", t, "\n")
            break
        }
        # don't do the last, unnecessary step...
        if (t<max(times,summary.times)) {
            # generation() really only needs N, not the whole population; but pass it in anyhow
            N <- generation(population,demography,N=N,t=t,...)
            t <- t+1
        }
    }
    outlist <- simulation(N=out,summaries=sums,times=times,summary.times=summary.times,t=t)
    return(outlist)
}

#' @describeIn simulate_pop Constructor for a simulation object.
simulation <- function (N, summaries, times, summary.times, t, ...) {
    out <- list(N=N, summaries=summaries, times=times, summary.times=summary.times, t=t, ...)
    class(out) <- "simulation"
    return(out)
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
                  rb <- do.call( raster::stack, list( pop$habitat )[rep(1,dim(sim$N)[3])] )
                  raster::values(rb)[pop$habitable] <- sim$N[,match(geno,pop$genotypes),]
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
#' @inheritParams simulate_pop
#' @param sim The previous \code{simulation} object.
#' @param append If TRUE, append results to results in previous simulation object.
#' @export
#' @return A simulation object, as \code{simulation}.
#' If summary functions are not the same, the initial ones are removed.
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
    new.sim <- simulate_pop( population=population, demography=demography, times=times, N=N, tinit=tinit, summary.times=summary.times, summaries=summaries, ... )
    if (append) {
        new.dims <- dim(new.sim$N)
        new.sim$N <- c( sim$N, new.sim$N )
        dim(new.sim$N) <- c( 0, 0, new.dims[3] ) + dim( sim$N )
        for (k in seq_along(new.sim$summaries)) {
            if (ncol(new.sim$summaries[[k]])!=ncol(sim$summaries[[k]])) { 
                warning("summaries must be the same to extend a simulation object") 
            } else {
                new.sim$summaries[[k]] <- rbind( sim$summaries[[k]], new.sim$summaries[[k]] )
            }
        }
        new.sim$times <- c(sim$times,new.sim$times)
        new.sim$summary.times <- c(sim$summary.times,new.sim$summary.times)
    }
    return(new.sim)
}


#' Crop a Simulation
#'
#' Crop a simulation object to a smaller region.
#'
#' @param sim A simulation object.
#' @param pop The corresponding population object.
#' @param extent A raster::extent object, or anything that can be used to crop a Raster*.
#' @export
#' @return A simulation object.  In the new object, `summaries` are NULL.
crop_simulation <- function (sim, pop, extent) {
    newpop <- crop_population(pop,extent,ind.return=TRUE)
    simulation(
                   N = sim$N[newpop$crop.inds,,,drop=FALSE],
                   summaries=NULL,
                   times=sim$times,
                   summary.times=numeric(0),
                   t=sim$t
               )
}


#' Animate a Simulation
#'
#' Sequentially plot the time steps in a simulation, optionally cropped.
#'
#' @param sim A simulation object.
#' @param pop The corresponding population object.
#' @param max.frames Maximum number of frames to plot.
#' @param zlim Range of values to represent by colors, consistently.
#' @param legend.width,legend.mar Parameters controlling how the legend is displayed.
#' @param animate If not FALSE, instead of plotting the result, package frames into an mp4 animation with this file name.
#' @param duration Duration of the animation in seconds.
#' @param ... Additional parameters passed to plot( ).
#' @export
plot.simulation <- function (sim, 
                             pop, 
                             max.frames=300,
                             zlim=range(sim$N,finite=TRUE),
                             pause=interactive(),
                             legend.width=2,
                             legend.mar=12,
                             animate=NULL,
                             duration=5,
                             ... ) {
    hab <- pop$habitat
    # even out the sampled times
    plot.t <- sort(unique(sim$times))
    plot.cols <- match(plot.t,sim$times)
    plot.inds <- as.numeric( cut( seq(
                                      min(plot.t),
                                      max(plot.t),
                                      length.out=min( max.frames,
                                                     (1+diff(range(plot.t)))/ceiling(1/min(diff(plot.t)[diff(plot.t)>0])) )
                                      ), 
                     breaks=c(plot.t[1]-1,plot.t) ) )
    if (!is.null(animate)) {
        png.base <- file.path( tempdir(), paste("frame",floor(1e6*runif(1)),"_",sep='') )
        png.files <- paste( png.base, "%03d.png", sep='')
        png( png.files, height=2*300, width=length(pop$genotypes)*2*300, res=300, pointsize=10 )
        on.exit(dev.off(),add=TRUE)
    }
    layout(t(1:ncol(sim$N)))
    for (k in plot.inds) {
        mains <- paste( c(sprintf("t=%d",floor(plot.t[k])),rep("",length(pop$genotypes)-1)), pop$genotypes )
        for (j in 1:ncol(sim$N)) {
            values(hab)[pop$habitable] <- sim$N[,j,plot.cols[k]]
            plot(hab,zlim=zlim,main=mains[j],...)
        }
        if (pause && is.null(animate) && length(locator(1))==0) { break }
    }
    if (!is.null(animate)) {
        system2( "ffmpeg", c("-y", "-v 8", "-i", png.files, "-r", length(plot.inds)/duration, animate) )
        unlink( list.files( dirname(png.base), paste(basename(png.base), "[0-9]*[.]png$", sep=''), full.names=TRUE ) )
    }
}


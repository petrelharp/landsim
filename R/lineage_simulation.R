#' Simulate Lineages
#'
#' Given a collection of lineages, all of the same genotype,
#' sample from their locations across previous generations,
#' given the entire set of reproductive numbers in those generation.
#'
#' @param lineages A two-column matrix with columns labeled 'location' and 'genotype'.
#' @param gens A list as returned by \code{simulate_pop(...,return=everything=TRUE)}.
#' @param num.alleles A vector giving the number of alleles of the type followed in each genotype.
#' @param ... Other parameters passed to \code{lineage_generation()}.
#' @export
#' @return An array of dimensions \code{c(dim(lineages),length(gens$times))}, in the SAME order as gens$times (so, [,,1] is the MOST DISTANT time.
simulate_lineages <- function (lineages, gens, num.alleles, ...) {

    if (length(gens$gens)!=length(gens$times)) {
        stop("Don't have 'gens' information: need to simulate with return.everything=TRUE.")
    }
    N <- gens$N[,,dim(gens$N)[3]]
    out <- array( dim=c(dim(lineages),length(gens$times)) )
    out[,,length(gens$times)] <- last.lins <- lineages
    for (k in rev(seq_along(gens$gens)[-1])) {
        N[] <- N + gens$gens[[k]]$death - gens$gens[[k]]$germination
        last.lins[] <- out[,,k-1] <- lineage_generation( last.lins, N=N,
                                          gen=gens$gens[[k]], num.alleles=num.alleles, ... )
    }
    return(out)
}

#' Plot Lineages on a Map
#'
#' Plot the path traced out by lineages on the habitat of the corresponding population.
#'
#' @param lineages An array indexed by (lineage,location-or-genotype,time).
#' @param pop A population object.
#' @param cols A vector of colors.
#' @param ... Other parameters to pass to \code{arrows()}.
#' @export
plot_lineages <- function (lineages,pop,cols=rainbow(dim(lineages)[1])) {
    lin.locs <- xyFromCell(pop$habitat,which(pop$habitable)[lineages[,1,]])
    dim(lin.locs) <- c( dim(lineages)[c(1,3)], 2 )

    plot(pop$habitat)
    for (k in 1:nrow(lin.locs)) {
        now.locs <- lin.locs[k,-1,]
        parent.locs <- lin.locs[k,-dim(lin.locs)[2],]
        do.arrows <- ( rowSums( abs(now.locs-parent.locs) )>0 )
        suppressWarnings( arrows( x0=now.locs[do.arrows,1], x1=parent.locs[do.arrows,1], 
                   y0=now.locs[do.arrows,2], y1=parent.locs[do.arrows,2], 
                   length=0.05, col=cols[k] ) )
    }
}


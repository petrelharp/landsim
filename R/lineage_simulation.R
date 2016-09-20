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



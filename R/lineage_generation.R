#' Simulate a Single Backwards Generation
#'
#' Given a collection of lineages, all of the same genotype,
#' sample from their location in the previous generation,
#' given the entire set of reproductive numbers in that generation.
#'
#' @param lineages A two-column matrix with columns labeled 'location' and 'genotype'.
#' @param N The population state at the start of the previous generation.
#' @param gen A list as returned by \code{generation(...,return.everything=TRUE)}.
#' @param num.alleles A vector giving the number of alleles of the type followed in each genotype.
#' @param debug Whether to run sanity checks.
#' @export
#' @return A matrix of the same form as \code{lineages}.
lineage_generation <- function (lineages, N, gen, num.alleles, debug=FALSE) {

    # current generation
    now.N <- N - gen$death + gen$germination

    # the youth will move
    is.youth <- ( 0 < rbinom( nrow(lineages), 
                         size=1, prob=gen$germination[lineages]/now.N[lineages] ) )
    youth <- lineages[is.youth,,drop=FALSE]

    if (any(is.youth)) {

        # do backwards seed migration for youth
        for (k in unique(youth[,'genotype'])) {
            do.these <- ( youth[,'genotype']==k )
            youth[do.these,'location'] <- backmigrate( youth[do.these,'location'], 
                                                  demog$seed.migration, weights=gen$seed.production[,k] )
        }

        # choose seed or pollen parent
        num.seed.parents <- gen$seeders[youth[,'location'],,drop=FALSE]
        num.pollen.parents <- gen$pollen[youth[,'location'],,drop=FALSE]
        # check this was possible
        if (debug) {
            stopifnot(all(rowSums(num.seed.parents)>0))
            stopifnot(all(rowSums(num.pollen.parents)>0))
        }

        parent.types <- backmeiosis( youth[,'genotype'], mating=demog$mating, weights=num.alleles, 
                                    seed=num.seed.parents, pollen=num.pollen.parents )
        colnames(parent.types) <- c("genotype","which")

        # copy back over lineages that chose 'seed' parents
        youth[,'genotype'] <- parent.types[,'genotype']
        is.pollen <- ( parent.types[,2]==2 )

        # backmigrate pollen parents
        for (k in unique(youth[is.pollen,'genotype'])) {
            do.these <- ( is.pollen & ( youth[,'genotype']==k ) )
            youth[do.these,'location'] <- backmigrate( youth[do.these,'location'], demog$pollen.migration, weights=N[,k] )
        }
    }

    # finally, find the previous generation's locations
    parents <- lineages
    parents[is.youth,] <- youth

    return(parents)
}


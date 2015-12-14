#' Apply a \code{demography} Object to a Population.
#'
#' This steps a population forwards by one generation,
#' according to the parameters in a \code{demography} object,
#' (which can also be passed in separately; see \code{demography}).
#'
#' @inheritParams demography
#' @param population
#' @param demography A \code{demography} object, containing the below parameters.
#' @param ... Additional parameters that will be passed to demographic functions.
#' @export
#' @return A new population of the same form as the old.
#'
generation <- function (
                       population,
                       demography,
                       prob.seed = demography$prob.seed,
                       fecundity = demography$fecundity,
                       prob.germination = demography$prob.germination,
                       prob.survival = demography$prob.survival,
                       pollen.migration = demography$pollen.migration,
                       seed.migration = demography$seed.migration,
                       genotypes = demography$genotypes,
                       mating = demography$mating,
                       ...
        ) {
    # various of these can be simple numbers or more complicated functions
    fun_or_number <- function (f) {
        if (mode(f)=="function") { f(population,...) } else { f }
    }
    # sample number of seed-producing individuals:
    M <- rbinom_raster( size=population, prob=fun_or_number(prob.seed) )
    # find mean pollen flux
    P <- migrate(population,pollen.migration)
    # mean seed production
    S <- seed_production(seeders=M,pollen=P,mating=mating,fecundity=fun_or_number(fecundity))
    # seed dispersal
    SD <- migrate(S,seed.migration)
    # new individuals
    G <- rpois_raster( SD * fun_or_number(prob.germination) )
    # deaths
    V <- rbinom_raster( size=population, prob=fun_or_number(prob.survival) )
    return(V+G)
}

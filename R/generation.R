#' Simulate a Single Generation.
#'
#' This steps a population forwards by one generation,
#' according to the parameters in a \code{demography} object,
#' (which can also be passed in separately; see \code{demography}).
#'
#' @inheritParams demography
#' @param population A \code{population} object, the initial state of the population.
#' @param demography A \code{demography} object, containing the below parameters.
#' @param N Optionally, a matrix of the same form as \code{population$N} (overrides \code{population$N} if present).
#' @param expected If TRUE, do no sampling, returning only expected values.
#' @param return.everything If TRUE, returns the number of seeders, amount of pollen, density of seeds, number of germinating seeds, and number of dying individuals.
#' @param ... Additional parameters that will be passed to demographic functions.
#' @export
#' @return A matrix of the same dimensions as \code{population$N},
#' unless \code{return.everything} is TRUE, in which case it is a list with entries:
#'
#'   seeders : number of seeding individuals
#'   pollen : quantity of pollen produced, by genotype of parent
#'   seed.production : quantity of fertilized seeds produced
#'   seeds.dispersed : quantity of seeds, post-dispersal
#'   germination : number of new seeds that germinate (survive)
#'   death : number of previously alive individuals that die
generation <- function (
                       population,
                       demography,
                       N=population$N,
                       prob.seed = demography$prob.seed,
                       fecundity = demography$fecundity,
                       prob.germination = demography$prob.germination,
                       prob.survival = demography$prob.survival,
                       pollen.migration = demography$pollen.migration,
                       seed.migration = demography$seed.migration,
                       genotypes = demography$genotypes,
                       mating = demography$mating,
                       expected = FALSE,
                       return.everything = FALSE,
                       ...
        ) {
    # various of these can be simple numbers or more complicated functions
    fun_or_number <- function (f) {
        if (mode(f)=="function") { f } else { function(N,t,...){f} }
    }
    # sample number of seed-producing individuals:
    if (!expected) {
        seeders <- rbinom_matrix( size=N, prob=fun_or_number(prob.seed)(N,t=t,...) )
    } else {
        seeders <- ( N * fun_or_number(prob.seed)(N,t=t,...) )
    }
    # find mean pollen flux
    pollen <- migrate(N,pollen.migration)
    # mean seed production
    seed.production <- seed_production(seeders=seeders,
                                       pollen=pollen,
                                       mating=mating,
                                       fecundity=fun_or_number(fecundity)(N,t=t,...) )
    # seed dispersal
    seeds.dispersed <- migrate(seed.production,seed.migration)
    # deaths
    if (!expected) {
        survivors <- rbinom_matrix( size=N, prob=fun_or_number(prob.survival)(N,t=t,...) )
    } else {
        survivors <- ( N * fun_or_number(prob.survival)(N,t=t,...) )
    }
    # have death occur before recruitment to allow just-vacated spots to be filled
    # (wouldn't be necessary if we had a seed bank)
    # new individuals
    if (!expected) {
        germination <- rpois_matrix( seeds.dispersed * fun_or_number(prob.germination)(N=survivors,t=t,...) )
    } else {
        germination <- ( seeds.dispersed * fun_or_number(prob.germination)(N=survivors,t=t,...) )
    }
    if (return.everything) {
        return( list(seeders=seeders, 
                     pollen=pollen, 
                     seed.production=seed.production, 
                     seeds.dispersed=seeds.dispersed,
                     germination=germination,
                     death=N-survivors
                 ) )
    } else {
        return(survivors+germination)
    }
}

#' Apply a \code{demography} Object to a Raster
#'
#' This steps a population, stored as a Raster*, forwards by one generation,
#' according to the parameters in a \code{demography} object,
#' (which can also be passed in separately; see \code{demography}).
#'
#' @inheritParams demography
#' @param population A Raster* object.
#' @param demography A \code{demography} object, containing the below parameters.
#' @param ... Additional parameters that will be passed to demographic functions.
#' @export
#' @return A new population of the same form as the old.
#'
generation_raster <- function (
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
    P <- migrate_raster(population,pollen.migration)
    # mean seed production
    S <- seed_production_raster(seeders=M,pollen=P,mating=mating,fecundity=fun_or_number(fecundity))
    # seed dispersal
    SD <- migrate_raster(S,seed.migration)
    # new individuals
    G <- rpois_raster( SD * fun_or_number(prob.germination) )
    # deaths
    V <- rbinom_raster( size=population, prob=fun_or_number(prob.survival) )
    return(V+G)
}

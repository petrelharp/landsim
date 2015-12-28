#' Find the Intrinsic Growth Rate for Each Genotype
#'
#' This computes the per capita growth rate for each genotype, separately, in a population of constant density.
#'
#' @param population A \code{population} object, the initial state (t=0) of the population.
#' @param demography A \code{demography} object, containing the below parameters.
#' @param density The density to compute growth rates at.
#' @param ... Additional parameters that will be passed to \code{generation()}.
#' @export
#' @return A matrix of the same form as \code{population$N}.
#' For more details, see \code{generation()}, with the option \code{expected=TRUE}.
intrinsic_growth <- function (population, demography, density=1, ...) {
    growth <- sapply( population$genotypes, function (geno) {
                population$N[] <- 0
                population$N[,geno] <- density
                rowSums( generation(population,demography, expected=TRUE, ... ) )/density
         } )
    colnames(growth) <- population$genotypes
    return(growth)
}


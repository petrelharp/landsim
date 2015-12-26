#' Construct a \code{demography} Object.
#'
#' This puts together the parameters necessary to specify a \code{demography} object,
#' specifying pollen and seed production and dispersal, and mortality.
#' The parameters \code{prob.seed}, \code{fecundity}, \code{prob.germination}, 
#' and \code{prob.survival} can be either Raster* objects, numbers, 
#' or functions of the current population state (a Raster* object).
#'
#' @param prob.seed Probability of producing seed, per individual.
#' @param fecundity Mean number of seeds produced, per individual, if seeds are produced.
#' @param prob.germination Probability of germination (establishment) per seed.
#' @param prob.survival Probability of survival per already-extant individual.
#' @param pollen.migration A \code{migration} object for pollen dispersal.
#' @param seed.migration A \code{migration} object for seed dispersal.
#' @param genotypes A character vector of genotypes.
#' @param mating An array giving probabilities for offspring genotypes given parental genotypes.
#' @param ... Other parameters that are included verbatim in the output object.
#' @export
#' @return A \code{demography} S3 object (just a named list).
#'
demography <- function (
                       prob.seed,
                       fecundity,
                       prob.germination,
                       prob.survival,
                       pollen.migration,
                       seed.migration,
                       genotypes,
                       mating=mating_tensor(genotypes),
                       ...)  {
    out <- c( list(
                prob.seed=prob.seed,
                fecundity=fecundity,
                prob.germination=prob.germination,
                prob.survival=prob.survival,
                pollen.migration=pollen.migration,
                seed.migration=seed.migration,
                genotypes=genotypes,
                mating=mating
            ), list(...) )
    class(out) <- "demography"
    return(out)
}

#' Set up a Demography Object for a Concrete Population
#'
#' Demography objects are agnostic about the map they work on,
#' but for application we need to link them to a particular one.
#' This function does that, by setting up any elements that are \code{migration} or \code{vital} objects.
#'
#' @param demog The demography object.
#' @param pop The population object.
#' @param ... Other parameters passed to \code{migration()}.
#' @export
#' @return A \code{demography} object, as before, but with but with any migration objects having a migration matrix.
#'
setup_demography <- function (demog,pop,...) {
    for (k in seq_along(demog)) {
        if (inherits(demog[[k]],"migration")) { demog[[k]] <- setup_migration(demog[[k]],pop,...) }
        if (inherits(demog[[k]],"vital")) { demog[[k]] <- setup_vital(demog[[k]],pop,...) }
    }
    return(demog)
}

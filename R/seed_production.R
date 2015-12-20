#' Mean Seed Production
#'
#' Find the mean seed production, by genotype,
#' given the local numbers of seeding individuals and pollen density,
#' by genotype.
#'
#' @param seeders A numeric matrix of numbers of seeding individuals, with number of columns equal to the number of genotypes.
#' @param pollen A numeric matrix of pollen density, with number of columns equal to the number of genotypes.
#' @param mating An array with probabilities of producing each genotype from each parental genotype.
#' @export
#' @return A matrix of the same form as \code{seeders}.
seed_production <- function (
                         seeders,
                         pollen,
                         mating,
                         fecundity=1
                     ) {
    stopifnot( all( dim(mating) == c( NCOL(seeders), NCOL(pollen), NCOL(seeders) ) ) )
    if (is.null(dim(seeders))) { dim(seeders) <- c(NROW(seeders),NCOL(seeders)) }
    if (is.null(dim(pollen))) { dim(pollen) <- c(NROW(pollen),NCOL(pollen)) }
    total.pollen <- rowSums(pollen)
    out <- 0.0 * seeders
    for (i in 1:dim(mating)[1]) {
        for (j in 1:dim(mating)[2]) {
            for (k in 1:dim(mating)[3]) {
                if ( mating[i,j,k]>0 ) {
                    lprod <- seeders[,i]*pollen[,j]/total.pollen
                    out[,k] <- out[,k]+mating[i,j,k]*lprod
                }
            }
        }
    }
    return(out*fecundity)
}

#' Mean Seed Production
#'
#' Find the Raster* of mean seed production, by genotype,
#' given the local numbers of seeding individuals and pollen density,
#' by genotype.
#'
#' @param seeders A Raster* of numbers of seeding individuals.
#' @param pollen A Raster* of pollen density.
#' @param mating An array with probabilities of producing each genotype from each parental genotype.
#' @export
#' @return A Raster* of the same form as the input.
seed_production_raster <- function (
                         seeders,
                         pollen,
                         mating,
                         fecundity=1
                     ) {
    total.pollen <- sum(pollen)  # this is the Raster method that acts like rowSums( )
    out <- 0.0 * do.call( brick, list(total.pollen)[rep(1,dim(mating)[3])] )
    # values(out)[!is.na(values(out))] <- 0  # this is done by the multiplication by 0
    for (i in 1:dim(mating)[1]) {
        for (j in 1:dim(mating)[2]) {
            for (k in 1:dim(mating)[3]) {
                if ( mating[i,j,k]>0 ) {
                    lprod <- seeders[[i]]*pollen[[j]]/total.pollen
                    out[[k]] <- out[[k]]+mating[i,j,k]*lprod
                }
            }
        }
    }
    names(out) <- dimnames(mating)[[3]]
    return(out*fecundity)
}

#' Comptue the Mating Tensor for a Set of Genotypes.
#'
#' This returns the three-dimensional array whose \code{[i,j,k]}th entry is the probability that maternal genotype \code{i} and paternal genotype \code{j} combine to make offspring genotype \code{k}
#'
#' @param maternal Character vector of genotypes.
#' @param paternal Character vector of genotypes. [default: same as maternal]
#' @param offspring Character vector of genotypes. [default: same as maternal]
#' @param sep Character separating the two alleles. [default: '']
#' @param phased Does order matter for the diploid genotypes?
#' @export
#' @return A three-dimensional array, with dimnames.
#' If \code{phased} is \code{FALSE} then the two paternal alleles are sorted before being combined together into the offspring diploid genotype.
mating_tensor <- function (maternal, paternal=maternal, offspring=maternal, sep='', phased=FALSE) {
    mating <- array( 0, dim=c(length(maternal),length(paternal),length(offspring)) )
    dimnames(mating) <- list( maternal, paternal, offspring )
    mat.zygotes <- strsplit(maternal,split=sep)
    pat.zygotes <- strsplit(paternal,split=sep)
    for (i in seq_along(maternal)) {
        for (j in seq_along(paternal)) {
                mating_fun <- if (phased) {
                        function (x,y) { paste(x,y,sep=sep) }
                    } else {
                        function (x,y) { paste( ifelse( x<y, x, y ), ifelse( x<y, y, x ), sep=sep ) }
                    }
                outputs <- outer( mat.zygotes[[i]],
                                  pat.zygotes[[j]],
                                mating_fun )
                mating[i,j,] <- tabulate( factor( outputs, levels=offspring ), nbins=length(offspring) ) / length(outputs)
                stopifnot(sum(mating[i,j,])==1)
        }
    }
    return(mating)
}

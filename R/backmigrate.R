#' Sample from Parental Locations
#'
#' This takes a set of locations and a nonnegative weighting vector,
#' and samples from the time-reversed migration probabilities, weighted by the weights.
#'
#' @param x A vector of locations (indices).
#' @param migration The \code{migration} object.
#' @param weights A nonnegative vector of weights.
#' @export
#' @return A vector of locations of the same form as \code{x}.
#' If the probability of migrating from i to j is M[i,j], then a sample at j will choose i
#' with probability equal to \code{weights[i] M[i,j]/sum(weights*M[,j])}.
backmigrate <- function ( x,
                      migration,
                      weights
                 ) {
    if (is.null(migration$M)) { stop("Migration matrix does not exist: need to use setup_demography, or migrate_raster()?") }

    # prob.step[i,k] will contain the probability of migrating from i to x[k]
    #  in however many steps we've taken so far
    prob.step <- sparseMatrix( i=x, j=seq_along(x), x=rep(1.0,length(x)), 
                          dims=c(nrow(migration$M),ncol=length(x)) )

    zero.weight <- 1-sum(migration$n.weights)
    if (zero.weight<0) { warning("n.weights sum to more than one.") }
    probs <- zero.weight * prob.step
    for (n in seq_along(migration$n.weights)) {
        # note M is row-stochastic, so %*% does the right thing
        prob.step <- migration$M %*% prob.step 
        if (migration$n.weights[n]!=0) {
            probs <- probs + migration$n.weights[n]*prob.step
        }
    }
    # return samples
    probs <- sweep( probs, 1, weights, "*" )
    out <- apply( probs, 2, function (pp) { sample.int( n=nrow(probs), size=1, prob=pp ) } )
    return( as(out,class(x)) )
}

#' Sample from Parental Genotypes
#'
#' This takes a mating tensor, a vector of genotypes, a vector of weights, 
#' and matrices of numbers of possible seed and pollen parents (with seed[k,j] and pollen[k,j] 
#' being the numbers of possible j-genotype parents of genotype[k]),
#' and samples from parental genotypes.  
#'
#' @param x A vector of genotypes.
#' @param mating A mating tensor.
#' @param weights A numeric vector of length equal to number of genotypes.
#' @param seed A matrix with dimensions \code{length(x)} by number of genotypes.
#' @param pollen A matrix with dimensions \code{length(x)} by number of genotypes.
#' @export
#' @return A \code{length(x)} by 2 matrix whose first column is the genotype index 
#' and the second column is 1 or 2, meaning:
#'   1 : seed (first dimension in mating tensor)
#'   2 : pollen (second dimension in mating tensor)
#'
#' For genotype=k, the probability 
#' of choosing parent genotype pair i,j
#' is proportional to 
#'   \code{X[i,j] = mating[i,j,k]*seed[i]*pollen[j]}
#' and given \code{i,j}, the probability of choosing the seed parent is equal to
#'   \code{weights[i]/(weights[i]+weights[j])} .
backmeiosis <- function ( x,
                         mating,
                         weights,
                         seed,
                         pollen
                     ) {
    n.genotypes <- dim(mating)[1]
    # premultiply by weights
    out <- matrix( nrow=length(x), ncol=2 )
    for (k in seq_along(x)) {
        X <- seed[k,] *  sweep( mating[,,x[k],drop=FALSE], 2, pollen[k,], "*" )
        ij <- arrayInd( sample.int(length(X), size=1, prob=X), .dim=dim(mating)[1:2] )
        ps <- sample.int( 2, size=1, prob=weights[ij] )
        out[k,] <- c( ij[ps], ps )
    }
    return(out)
}

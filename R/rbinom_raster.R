#' Sample Binomial Values in a Raster* or Population
#'
#' Returns a Raster* object whose values are random, Binomial(N,prob.seed),
#' or just the corresponding values.
#'
#' @param size Number, or Raster* object of sample sizes (must be integers).
#' @param prob Number, or Raster* object of probabilities.
#' @param only.values Return only the vector of non-NA values?
#' @export
#' @return A Raster* object of the same form as the input,
#' or if \code{only.values} is \code{TRUE} the vector of values.
#' At least one of \code{size} or \code{prob} must be Raster* objects.
rbinom_raster <- function (size, prob, only.values=FALSE) {
    if (!inherits(size,"Raster")) {
        out <- prob
        prob <- raster::values(prob)
        size <- rep_len(size,length(out))
    }
    if (!inherits(prob,"Raster")) {
        out <- size
        size <- raster::values(size)
        prob <- rep_len(prob,length(out))
    }
    nonmissing <- ( !is.na(prob) ) & ( !is.na(size) )
    if (any(nonmissing)) {
        x <- rbinom(sum(nonmissing),size=size[nonmissing],prob=prob[nonmissing])
        if (only.values) {
            return(x)
        } else {
            raster::values(out)[nonmissing] <- x
        }
    }
    return(out)
}

#' Sample from a Binomial, Preserving Dimension
#'
#' Works as `rbinom`, except that it preserves the dimensions, and dimnames of the input(s).
#'
#' @param size Integer vector, matrix, or array.
#' @param prob Numeric vector, matrix, or array.
#' @export
#' @return An integer object of the same dimensions as the input (recycled to match each other),
#' and the same dimnames if these can be matched.
rbinom_matrix <- function (size, prob) {
    dims <- c( max(NROW(size),NROW(prob)), max(NCOL(size),NCOL(prob)) )
    out <- rbinom( prod(dims), size=as.numeric(size), prob=as.numeric(prob) )
    dim(out) <- dims
    if (!is.null(dim(prob)) && dim(prob)==dims) { dimnames(out) <- dimnames(prob) }
    if (!is.null(dim(size)) && dim(size)==dims) { dimnames(out) <- dimnames(size) }
    return(out)
}

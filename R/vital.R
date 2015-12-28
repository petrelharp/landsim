#' Set Up a Vital Rates Function With Parameters
#'
#' Takes a function and an arbitrary collection of other parameters that are used in that function,
#' and packages the other parameters into the environment of the function.
#'
#' @param fun The function, whose first argument is \code{N}, a matrix of genotype counts.
#' @param ... Other variables used by \code{fun}.
#' @export
vital <- function (
                     fun, 
                     ...
                 ) {
    fun.env <- list2env(list(...),parent=environment(fun))
    environment(fun) <- fun.env
    class(fun) <- c("vital","function")
    return(fun)
}


#' Extract Parameters from Vital Function
#'
#' Unpackages parameters from the result of \code{vital()} (or more generally from the environment of a function).
#'
#' @param fun The function.
#' @param ... Other variables used by \code{fun}.
#' @export
get_vital <- function ( fun, ...) {
    mget(...,envir=environment(fun),inherits=FALSE)
}


as.list.vital <- function (x) { c( as.list(environment(x)), list(fun=x) ) }

"$.vital" <- function (x, name) { get(name,envir=environment(x),inherits=FALSE) }
"$<-.vital" <- "[[<-.vital" <- function (x, i, value) { assign(i,value,envir=environment(x),inherits=FALSE); x }
"[[.vital" <- function (x, i) { get(i,envir=environment(x),inherits=FALSE) }
"[.vital" <- function (x, i) { mget(i,envir=environment(x),inherits=FALSE) }
names.vital <- function (x) { ls(environment(x)) }


#' Set up a Vital Function for a Concrete Population
#'
#' A vital function may include migration objects,
#' which for use we need to link to particular population
#' (precomputing the migration matrix).
#' *In addition*, every element of the population is added
#' to the environment of the function,
#' so that the function can depend on objects only visible in the population.
#'
#' @param fun The vital function.
#' @param pop The population object.
#' @param ... Other parameters passed to \code{migration()}.
#' @export
#' @return A \code{vital} function, as before, but with any migration objects having a migration matrix.
setup_vital <- function (fun,pop,...) {
    for (k in names(fun)) {
        if (inherits(fun[[k]],"migration")) { fun[[k]] <- setup_migration(fun[[k]],pop,...) }
    }
    for (n in names(pop)) {
        fun[[n]] <- pop[[n]]
    }
    return(fun)
}

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


#' Set Up a Vital Rates Function With Parameters
#'
#' Takes a function and an arbitrary collection of other parameters that are used in that function,
#' and packages the other parameters into the environment of the function.
#'
#' @param x The function.
#' @param ... Other variables used by \code{fun}.
#' @export
as.list.vital <- function (x,...) { c( as.list(environment(x)), list(fun=x) ) }

#' Assign or Extract Items from the Environment of a Vital Function.
#'
#' Extraction methods return NULL if the item is not defined, just like for lists.
#'
#' @param x,object The vital function.
#' @param i,name The name of the item to extract.
#' @param value The new value to assign.
#' @export
getElement.vital <- function(object,name) { 
    if (exists(name,environment(object))) {
        get(name,envir=environment(object),inherits=FALSE) 
    } else { NULL }
}

#' @describeIn getElement.vital Extract items from the environment of a vital function.
#' @export
"$.vital" <- function (x, name) { getElement.vital(x,name) }

#' @describeIn getElement.vital Assign to environment of a vital function.
#' @export
"$<-.vital" <- function (x, i, value) { assign(i,value,envir=environment(x),inherits=FALSE); x }

#' @describeIn getElement.vital Assign to environment of a vital function.
#' @export
"[[<-.vital" <- function (x, i, value) { assign(i,value,envir=environment(x),inherits=FALSE); x }

#' @describeIn getElement.vital Extract items from the environment of a vital function.
#' @export
"[[.vital" <- function (x, i) { getElement.vital(x,i) }

#' @describeIn getElement.vital Extract multiple items from the environment of a vital function.
#' @export
"[.vital" <- function (x, i) { mget(i,envir=environment(x),inherits=FALSE,ifnotfound=list(NULL)) }

#' Extract Names of Objects in Environment of a vital Function
#'
#' @param x The vital function.
#' @export
names.vital <- function (x) { ls(environment(x)) }



#' Read in a JSON Configuraiton Object
#'
#' This reads in JSON, using jsonlite::fromJSON,
#' but then does additional parsing of lists with names:
#'    "R"        : evaluate this as an R expression.
#'    "migrate"  : pass this to \code{make_migration()}
#'    "function" : evaluate this as the body of an R function whose arguments are "(N,...)"
#'                 Lines *must* be separated by semicolons.
#' This parsing is done with everything else in the list available as context.
#'
#' @export
#' @param json A file name or character string.
#' @param envir An environment to evaluate expressions in.
#' @return A named list.
#'
json_config <- function (json,
                         envir=environment(),
                         simplifyVector=TRUE
                     ) {
    out <- jsonlite::fromJSON(gsub("\n","",json),simplifyVector=simplifyVector)
    for (k in seq_along(out)) {
        if ( is.list(out[[k]]) && (length(out[[k]])==1) ) {
            out[[k]] <- switch( names(out[[k]]),
                    "R" = eval(parse(text=out[[k]][[1]])),
                    "function" = out[[k]] <- eval(parse(text=paste("function (N,...) {",out[[k]][[1]],"}"))),
                    "migration" = do.call( make_migration, c( out[[k]][[1]], list(do.M=TRUE, population=population) ) ),
                    out[[k]]
                )
        }
    }
    return(out)
}


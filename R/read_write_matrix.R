
#' Write out a Sparse Matrix in Triplet Form
#'
#' Writes out a file (using \code{write.table}) containing columns
#'  row  col   value
#'   i    j     M[i,j]
#' for every nonzero entry of M.
#'
#' @param M The matrix.
#' @param file The file name to write out to.
#' @param row.names Whether to write out row indices.
#' @param ... Additional parameters to be passed to \code{write.table}.
#' @return The result of \code{write.table}.
write_triplets <- function ( M, file, row.names=FALSE,... ) {
    M <- as( M, "dgCMatrix" )
    dothese <- which( M@x > 0 )
    Mtab <- cbind( row=1L+M@i[dothese], col=p.to.j(M@p)[dothese], value=M@x[dothese] )
    write.table( Mtab, file=file, row.names=row.names, ... )
}

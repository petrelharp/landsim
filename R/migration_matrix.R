#' Weighted "Adjacency" Matrix from a Raster Layer.
#'
#' This returns the (sparse) Matrix giving the pseudo-adjacency matrix 
#' from and to for specified cells in the given Raster, with weights.
#'
#' @param population Population object or Raster*.
#' @param accessible Logical vector indicating which cells in `habitat` that migrants will attempt to go to.
#' @param migration Migration object; is overridden by other parameters below.
#' @param kern Weighting kernel applied to distances.
#' @param sigma Distance scaling for kernel.
#' @param radius Maximum distance away to truncate the kernel.
#' @param from The indices of the cells corresponding to rows in the output matrix. [default: non-NA cells]
#' @param to The indices of the cells corresponding to columns in the output matrix. [default: same as \code{from}]
#' @param normalize If not NULL, migration weights to accessible locations are set to sum to this (see details).
#' @export
#' @return A sparse Matrix \code{M}, where for each pair of cells \code{i} in \code{from} and \code{i} in \code{to},
#' the value \code{M[i,j]} is
#'    kern( d[i,j]/sigma )
#' or, if \code{kern="gaussian"},
#'    exp( - d[i,j]^2/sigma^2 ) / ( 2 pi sigma^2 ) * (area of a cell)
#' where \code{d[i,j]} is the distance from \code{from[i]} to \code{to[j]},
#' if this is greater than \code{min.prob}, and is 0 otherwise.
#'
#' Migration is possible to each \code{accessible} cell in the Raster*.
#' If \code{normalize} is not NULL, then for each cell, the set of migration weights to each other accessible cell
#' within distance \code{radius} is normalized to sum to \code{normalize}.
#' The resulting matrix may not have rows summing to \code{normalize},
#' however, if \code{from} is a subset of the \code{accessible} ones.
#'
#' If not normalized, the kernel is multiplied by the area of the cell (taken to be approximately constant).
#'
#' The usual way to use this is to call \code{migration_matrix(pop,mig)},
#' where \code{pop} is a \code{population} object, which contains both the Raster 
#' and the vector of which locations are accessible;
#' and \code{mig} is a \code{migration} object containing the kernel, migration radius, etcetera.
#'
#' Alternatively, \code{population} can be a RasterLayer,
#' in which case by default all non-NA cells are accessible.
#'
#' Inaccessible cells included in \code{to},
#' have zero migration outwards.
#' Inaccessible cells included in \code{from} behave as usual.
migration_matrix <- function (population,
                              accessible,
                              migration=list(sigma=1,normalize=1),
                              kern=migration$kern,
                              sigma=migration$sigma,
                              radius=migration$radius,
                              normalize=migration$normalize,
                              from,
                              to=from
                 ) {
    # Fill in default values.
    if (inherits(population,"population")) {
        if (missing(accessible)) { accessible <- population$accessible }
        if (missing(from)) { from <- which(population$habitable) }
        population <- population$habitat
    } else if (inherits(population,"Raster")) {
        if (missing(accessible)) { accessible <- !is.na(raster::values(population)) }
        if (missing(from)) { from <- which(!is.na(raster::values(population))) }
    }
    if (!is.integer(from)) { stop("migration_matrix: 'from' must be integer-valued (not logical).") }
    if (!is.logical(accessible)) { stop("migration_matrix: 'accessible' must be logical (not a vector of indices).") }
    kern <- get_kernel(kern)
    area <- prod(raster::res(population))
    cell.radius <- ceiling(radius/raster::res(population))
    directions <- matrix( 1, nrow=2*cell.radius[1]+1, ncol=2*cell.radius[2]+1 )
    directions[cell.radius[1]+1,cell.radius[2]+1] <- 0
    # columns are (from,to)
    # does NOT include diagonal
    use.to <- accessible[to]  # only want to go to accessible 'to' locations
    to.acc <- to[use.to]
    ij <- raster::adjacent(population,cells=from,target=to.acc,directions=directions,pairs=TRUE)
    from.pos <- raster::xyFromCell(population,cell=from)
    to.pos <- raster::xyFromCell(population,cell=to.acc)
    # add on the diagonal
    both.ij <- intersect(from,to.acc)
    ii <- match(c(ij[,1],both.ij),from)
    jj <- match(c(ij[,2],both.ij),to.acc)
    M <- Matrix::sparseMatrix( 
            i = ii,
            j = which(use.to)[jj], # map back to index in all given 'to' values
            x = kern(raster::pointDistance(from.pos[ii,],to.pos[jj,],lonlat=FALSE,allpairs=FALSE)/sigma)*area/sigma^2,
            dims=c(length(from),length(to))
        )
    if (!is.null(normalize)) {
        # M <- (normalize/Matrix::rowSums(M)) * M
        # # this is twice as quick:
        M@x <- (normalize/Matrix::rowSums(M)[1L+M@i]) * M@x
    }
    if (any(!is.finite(M@x))) { warning("Some values of migration matrix are not finite.") }
    return(M)
}

# helper function to find column indices of dgCMatrix objects
p.to.j <- function (p) { rep( seq.int( length(p)-1 ), diff(p) ) }

#' Subset a Migration Matrix to Match a Smaller Raster
#'
#' This converts a migration matrix computed for one raster
#' into a migration matrix for another raster that must be a subset of the first.
#' This does not do any renormalization.
#'
#' @param M The migration matrix.
#' @param old The original Raster* the migration matrix was computed for.
#' @param new The Raster* for the new migration matrix.
#' @param from.old The indices of the "from" cells corresponding to rows in the original migration matrix.
#' @param to.old The indices of the "to" cells corresponding to columns in the original migration matrix.
#' @param from.new The indices of the "from" cells corresponding to rows in the resulting migration matrix.
#' @param to.new The indices of the "to" cells corresponding to columns in the resulting migration matrix.
#' @export
#' @return A migration matrix.  See \code{migrate}.
subset_migration <- function (M, old, new, 
                              from.old=which(!is.na(raster::values(old))), to.old=from.old, 
                              from.new=which(!is.na(raster::values(new))), to.new=from.new
                         ) {
    if (any(raster::res(old)!=raster::res(new))) { stop("Resolutions must be the same for the two layers to subset a migration matrix.") }
    if ( raster::xmin(new)<raster::xmin(old) || raster::xmax(new)>raster::xmax(old) || raster::ymin(new)<raster::ymin(old) || raster::ymax(new)>raster::ymax(old) ) {
        stop("New layer must be contained in the old layer to subset a migration matrix.")
    }
    from.locs <- raster::xyFromCell(new,from.new)
    from.inds <- match( raster::cellFromXY(old,from.locs), from.old )
    to.locs <- raster::xyFromCell(new,to.new)
    to.inds <- match( raster::cellFromXY(old,to.locs), to.old )
    return(M[from.inds,to.inds])
}



#' Take a Geometric Power of a Matrix
#'
#' For a matrix M, compute (1-p) ( I + p M + p^2 M^2 + ... ) = (1-p)(I-pM)^{-1} .
#' If M is a transition matrix, the result is the distribution after taking a Geometric(p) number of steps.
#'
#' @param M The migration matrix.
#' @param p The parameter in the geometric distribution (must be between 0 and 1).
#' @param eps The numerical tolerance.
#' @param n The number of terms to use, minus one (defaults so that p^n = eps).
#' @export
#' @return A matrix of the same form as M.
geo_power <- function (M, p, eps=1e-8, n=log(eps)/log(p)) {
    if (p==0) { n <- 0 }
    out <- pM <- p * M
    diag(out) <- diag(out) + 1
    if (n>1) {
        for (k in 2:n) {
            pM <- p * M %*% pM
            out <- out + pM
        }
    }
    return( (1-p)*out )
}

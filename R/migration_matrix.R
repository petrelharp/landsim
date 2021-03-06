#' Create Migration Matrix from a Migration object and a Population object
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
#' @param discretize Whether or not to "discretize" the kernel first using the \code{discretize_kernel} function.
#' @param disc.fact The discretization factor, if used (see \code{discretize_kernel}).
#' @param n.chunks Number of chunks to compute the matrix in (can help with memory usage).
#' @export
#' @return A sparse Matrix \code{M}, where for each pair of cells \code{i} in \code{from} and \code{i} in \code{to},
#' the value \code{M[i,j]} is
#'    kern( d[i,j]/sigma )
#' or, if \code{kern="gaussian"},
#'    exp( - d[i,j]^2/sigma^2 ) / ( 2 pi sigma^2 ) * (area of a cell)
#' where \code{d[i,j]} is the distance from \code{from[i]} to \code{to[j]},
#' if this is greater than \code{min.prob}, and is 0 otherwise.
#'
#' Summary of how to encode different types of boundary: set absorbing boundary elements to accessible, 
#' but do not include them in from (and to). Mark reflecting (non-)boundaries as not accessible. 
#' External boundaries will be reflecting, so if they should be absorbing, you need to extend the raster and set values appropriately.
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
                              from=which(accessible),
                              to=from,
                              discretize=migration$discretize,
                              disc.fact=10,
                              n.chunks=1
                 ) {
    # Fill in default values.
    if (inherits(population,"population")) {
        if (missing(accessible)) { accessible <- population$accessible }
        population <- population$habitat
    } else if (inherits(population,"Raster")) {
        if (missing(accessible)) { accessible <- !is.na(raster::values(population)) }
    }
    if (!is.integer(from)) { stop("migration_matrix: 'from' must be integer-valued (not logical).") }
    if (!is.logical(accessible)) { stop("migration_matrix: 'accessible' must be logical (not a vector of indices).") }
    kern <- get_kernel(kern)
    if (!is.null(discretize) && discretize) { 
        if (!is.null(migration$disc.fact)) { disc.fact <- migration$disc.fact }
        kern <- discretize_kernel( kern, res=raster::res(population), radius=radius, sigma=sigma, fact=disc.fact )
    }
    area <- prod(raster::res(population))
    cell.radius <- ceiling(as.numeric(radius)/raster::res(population))
    directions <- matrix( 1, nrow=2*cell.radius[1]+1, ncol=2*cell.radius[2]+1 )
    directions[cell.radius[1]+1,cell.radius[2]+1] <- 0
    M <- Matrix::sparseMatrix( i = numeric(0), j = numeric(0), x = numeric(0), dims=c(length(from),length(to)))
    dk <- ceiling(length(from)/n.chunks)
    for (k in 1:n.chunks) {
        these.from <- from[seq(from=min((k-1)*dk+1,length(from)), to=min(k*dk,length(from)))]
        # columns are (from,to)
        # does NOT include diagonal
        use.to <- accessible[to]  # only want to go to accessible 'to' locations
        to.acc <- to[use.to]
        ij <- raster::adjacent(population, cells=these.from, target=to.acc, 
                               directions=directions, pairs=TRUE)
        if (is.null(dim(ij)) && length(ij)) {
            # adjacent returns a list in the case of only one adjacent cell =( =(
            dim(ij) <- c(1,2)
        }
        # add on the diagonal
        both.ij <- intersect(these.from, to.acc)
        # it is NOT a mistake that from.pos does not need indexing by ii yet to.pos does by jj!
        from.pos <- raster::xyFromCell(population, cell=these.from)[match(c(ij[,1],both.ij), these.from),]
        to.pos <- raster::xyFromCell(population, cell=to.acc)
        # (i,j) in the *full* matrix
        ii <- match(c(ij[,1],both.ij), from)
        jj <- match(c(ij[,2],both.ij), to.acc)
        # raster::pointDistance just does sqrt(dx^2+dy^2) 
        M[cbind(ii, which(use.to)[jj])] <- area/sigma^2 * kern( sqrt(rowSums((from.pos-to.pos[jj,])^2))/sigma )
    }
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

#' Subset a Migration Matrix to Match a Smaller Raster of the Same Resolution
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

#' Aggregate a Migration Matrix to Match a Coarser Raster
#'
#' This converts a migration matrix computed for one raster
#' into a migration matrix for another raster that is coarser than the first.
#' This function does no re-normalization.
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
#' The resulting matrix has [i,j]th entry equal to the average over all rows in the old matrix that correspond to `i`
#' of the sum of columns in that row corresponding to j.  In other words, if the old matrix M is indexed by [u,v],
#' then the new matrix's [i,j]th entry is:
#'     (1/n(j)) sum_{v: j(v)=j} sum_{u : i(u)=i}  M[u,v] 
#' where n(j) = #{ v : j(v)=j }.
aggregate_migration <- function (M, old, new, 
                              from.old=which(!is.na(raster::values(old))), 
                              to.old=from.old, 
                              from.new=which(!is.na(raster::values(new))), 
                              to.new=from.new
                         ) {
    # an easy mistake
    if (is.logical(from.old) || is.logical(to.old) || is.logical(from.new) || is.logical(to.new)) {
        stop("From and to must be indices, not a boolean vector.")
    }
    # which new cells do old cell centers fall in
    ## for 'from'
    from.old.locs <- raster::xyFromCell(old,from.old)
    from.old.in.new <- match( raster::cellFromXY(new,from.old.locs), from.new )
    ## and 'to'
    to.old.locs <- raster::xyFromCell(old,to.old)
    to.old.in.new <- match( raster::cellFromXY(new,to.old.locs), to.new )
    if ( (length(from.old.in.new)==0) || (length(to.old.in.new)==0) ) {
        stop( "Resulting matrix would be empty." )
    }
    if ( any(is.na(from.old.in.new)) || any(is.na(from.old.in.new)) ) {
        stop( "New raster does not cover old raster.")
    }
    # output M
    from.projmat <- Matrix::sparseMatrix( i=from.old.in.new, j=seq_along(from.old), x=rep(1.0,length(from.old)), dims=c(length(from.new),nrow(M)) )
    from.projmat@x <- 1/Matrix::rowSums(from.projmat)[1L+from.projmat@i]  # make this an averaging matrix (over rows)
    to.projmat <- Matrix::sparseMatrix( i=to.old.in.new, j=seq_along(to.old), x=rep(1.0,length(to.old)), dims=c(length(to.new),ncol(M)) )
    new.M <- Matrix::tcrossprod( from.projmat %*% M, to.projmat )
    return(new.M)
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

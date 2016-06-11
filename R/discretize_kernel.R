#' Create a "Discretized" Migration Function from a Continuous One
#'
#' Given a function kern(r) that gives per-individual migration rates and a resolution, 
#' returns the (approximate) function giving per-individual migration rates *per unit area*
#' from squares of side length `resolution` centered at distance r.  (This depends on the 
#' direction the cell is in, but the result averages over directions.)
#'
#' @param kern The migration rate function.
#' @param res The resolution of the grid.
#' @param radius The maximum distance to extrapolate out to.
#' @param sigma Distance scaling for kernel.
#' @param fact The number of divisions of each grid cell side when approximating the integral.
#' @param sim Whether to obtain the result by Monte Carlo integration (currently only implemented for the Gaussian kernel).
#' @export
#' @return A function as output by `approxfun()`, determined by numerical integration.
#' If the required number of grid cells in the approximation is too large,
#' falls back on Monte Carlo integration.
#'
#' The result is designed to be used in migration_matrix, which scales input to the kernel by sigma
#' and assumes that the kernel is a density function, so multiplies the value by the area of a cell; 
#' so the resulting function scales as a density function, and takes input in units of sigma.
#'
#' For good results, the 'fine' grid used in numerical integration should be smaller than sigma;
#' since the fine grid has a scale of res/fact, you should set fact >= 3 * res / sigma, say.
discretize_kernel <- function (kern, 
                               res, 
                               radius, 
                               sigma=1, 
                               fact=10,
                               sim=(kern=="gaussian"&&fact>10)
                           ) {
    cell.radius <- max(ceiling(radius/res))
    stencil <- raster::raster( xmn=-(cell.radius+0.5)*res[1], 
                               xmx=(cell.radius+0.5)*res[1], 
                               ymn=-(cell.radius+0.5)*res[2],
                               ymx=(cell.radius+0.5)*res[2],
                               resolution=res,
                               # without CRS string coordinates are lonlat not meters
                               crs=sp::CRS("+proj=utm +datum=NAD83 +units=m") )
    distvec <- values( raster::distanceFromPoints( stencil, c(0,0) ) )
    center.loc <- raster::cellFromXY(stencil,c(0,0))
    if (!sim) {
        stencil.fine <- raster::disaggregate(stencil, fact=fact)
        distvec.fine <- raster::values( raster::distanceFromPoints( stencil.fine, c(0,0) ) )
        # which new cells do old cell centers fall in
        fine.locs <- raster::xyFromCell(stencil.fine,seq_along(stencil.fine))
        fine.in.coarse <- raster::cellFromXY(stencil,fine.locs)
        M.fine <- migration_matrix( stencil.fine, 
                           accessible=rep(TRUE,length(stencil.fine)),
                           kern=kern, 
                           sigma=sigma,
                           radius=radius, 
                           normalize=NULL,
                           from=which(fine.in.coarse==center.loc), 
                           to=seq_along(stencil.fine) )
        M.aggr <- aggregate_migration( M.fine, 
                           old=stencil.fine, 
                           new=stencil,
                           from.old=which(fine.in.coarse==center.loc), 
                           to.old=seq_along(stencil.fine),
                           from.new=center.loc, 
                           to.new=seq_along(stencil) )
        M.probs <- as.vector(M.aggr)
    } else {
        npts <- 1e6
        xy0 <- cbind( 
                 x=(runif(npts)-0.5)*cell.radius*res[1],
                 y=(runif(npts)-0.5)*cell.radius*res[2]
             )
        xy1 <- cbind(
                 x=xy0[,"x"] + rnorm(npts,sd=sigma),
                 y=xy0[,"y"] + rnorm(npts,sd=sigma)
             )
        M.probs <- as.vector(tabulate( raster::cellFromXY( stencil, SpatialPoints( xy1, proj4string=sp::CRS(sp::proj4string(stencil)) ) ), nbins=prod(dim(stencil)) ))/npts
    }
    #   must divide by the target area to get a density function
    # outfun <- approxfun( distvec/sigma, as.vector(M.aggr)/prod(res/sigma), rule=2 )
    out.ord <- order(distvec) # y must be increasing for method='hyman'
    # use.these <- c(TRUE, diff(distvec[out.ord])/diff(range(distvec)) > 1e-3 )
    use.these <- TRUE   # if not using 'hyman' don't need to remove redundants
    outfun <- splinefun( distvec[out.ord][use.these]/sigma, 
                        M.probs[out.ord][use.these]/prod(res/sigma), 
                        method="monoH.FC" )
    attr(outfun,"res") <- res
    attr(outfun,"sigma") <- sigma
    return( outfun )
}

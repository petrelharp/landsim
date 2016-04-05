#' Create a "Discretized" Migration Function from a Continuous One
#'
#' Given a function kern(x,y) that gives per-individual migration rates and a resolution, 
#' returns the (approximate) function giving per-individual migration rates *per unit area*
#' from squares of side length `resolution` centered at (x,y).
#'
#' @param kern The migration rate function.
#' @param res The resolution of the grid.
#' @param radius The maximum distance to extrapolate out to.
#' @param sigma Distance scaling for kernel.
#' @param fact The number of divisions of each grid cell side when approximating the integral.
#' @export
#' @return A function as output by `approxfun()`.
#'
#' The result is designed to be used in migration_matrix, which scales input to the kernel by sigma
#' and assumes that the kernel is a density function, so multiplies the value by the area of a cell; 
#' so the resulting function scales as a density function, and takes input in units of sigma.
discretize_kernel <- function (kern, res, radius, sigma=1, fact=10) {
    cell.radius <- max(ceiling(radius/res))
    stencil <- raster::raster( xmn=-(cell.radius+0.5)*res[1], 
                               xmx=(cell.radius+0.5)*res[1], 
                               ymn=-(cell.radius+0.5)*res[2],
                               ymx=(cell.radius+0.5)*res[2],
                               resolution=res,
                               # without CRS string coordinates are lonlat not meters
                               crs=CRS("+proj=utm+datum=NAD83+units=m") )
    distvec <- values( raster::distanceFromPoints( stencil, c(0,0) ) )
    stencil.fine <- raster::disaggregate(stencil, fact=fact)
    distvec.fine <- raster::values( raster::distanceFromPoints( stencil.fine, c(0,0) ) )
    # which new cells do old cell centers fall in
    fine.locs <- raster::xyFromCell(stencil.fine,seq_along(stencil.fine))
    fine.in.coarse <- raster::cellFromXY(stencil,fine.locs)
    center.loc <- raster::cellFromXY(stencil,c(0,0))
    M.fine <- migration_matrix( stencil.fine, 
                       accessible=rep(TRUE,length(stencil.fine)),
                       kern=kern, 
                       sigma=sigma,
                       radius=cell.radius, 
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
    #   must divide by the target area to get a density function
    # outfun <- approxfun( distvec/sigma, as.vector(M.aggr)/prod(res/sigma), rule=2 )
    outfun <- splinefun( distvec/sigma, as.vector(M.aggr)/prod(res/sigma), method="hyman" )
    attr(outfun,"res") <- res
    attr(outfun,"sigma") <- sigma
    return( outfun )
}

context("testing aggregate migration function")

set.seed(42)
habitat <- raster(xmn=0, xmx=5, ymn=0, ymx=5, 
      resolution=1,
      crs="+proj=utm +zone=11 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
values(habitat)[] <- 1.0

naive.migration <- migration( kern="gaussian",
                             sigma=0.5,
                             radius=5,
                             normalize=NULL )
center.loc <- as.integer(cellFromXY( habitat, c(0,0) ))
disag.fact <- 5
habitat.fine <- disaggregate(habitat, fact=disag.fact)
fine.locs <- xyFromCell(habitat.fine,seq_len(prod(dim(habitat.fine))))
fine.in.coarse <- cellFromXY(habitat,fine.locs)
M.fine <- migration_matrix( habitat.fine, 
                           kern=naive.migration$kern, 
                           sigma=naive.migration$sigma, 
                           radius=naive.migration$radius, 
                           normalize=naive.migration$normalize )
                           # from=which(fine.in.coarse==center.loc), to=seq_len(prod(dim(habitat.fine))) )
M.aggr <- aggregate_migration( M.fine, old=habitat.fine, new=habitat )
                           # from.old=which(fine.in.coarse==center.loc), to.old=seq_len(prod(dim(habitat.fine))),
                           # from.new=center.loc, to.new=seq_len(prod(dim(habitat))) )


# aggregate_migration is supposed to sum over columns and average over rows
fine.mat <- as.matrix( habitat.fine )
fine.row <- row(fine.mat)
fine.col <- col(fine.mat)

coarse.mat <- as.matrix(habitat)
coarse.row <- row(coarse.mat)
coarse.col <- col(coarse.mat)

tf <- function (n) { as.vector(1+((n-1)%/%disag.fact)) } 
coarse.ind <- coarse.mat
coarse.ind[] <- 1:25
fine.to.coarse <- coarse.ind[ cbind( tf(fine.row), tf(fine.col) ) ]
dim(fine.to.coarse) <- dim(fine.mat)

expect_true( all( table(fine.to.coarse) == disag.fact^2 ) )

M.colsums <- t( apply( M.fine, 1, tapply, fine.to.coarse, sum ) )
M.rowmeans <- apply( M.colsums, 2, tapply, fine.to.coarse, mean )

expect_true( all( abs(M.rowmeans-M.aggr) < 10*.Machine$double.eps ) )

testthat::context("Agreement between migration_matrix and migrate.")

na_to_zero <- function (x) {
    x[is.na(x)] <- 0
    return(x)
}

# check we understand focal()

habitat <- raster(xmn=0, xmx=5, ymn=0, ymx=5, 
      resolution=1,
      crs="+proj=utm +zone=11 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
values(habitat) <- sample( c(1,2), length(habitat), replace=TRUE )
habitat.NA <- habitat
values(habitat.NA)[c(1,12)] <- NA

w <- matrix(1/16,nrow=3,ncol=3)
w[2,2] <- 1/2
hf1 <- focal( habitat, w=w, pad=TRUE, na.rm=TRUE, padValue=0 )
hf1.NA <- focal( habitat.NA, w=w, pad=TRUE, na.rm=TRUE, padValue=0 )

hmat <- matrix(values(habitat),nrow=nrow(habitat),ncol=ncol(habitat))
hmat.NA <- matrix(values(habitat.NA),nrow=nrow(habitat),ncol=ncol(habitat))
testthat::expect_equivalent(as.numeric(hmat),values(habitat))
testthat::expect_equivalent(as.numeric(hmat.NA),values(habitat.NA))

hmat.pad <- rbind( 0, cbind( 0, hmat, 0 ), 0 )
hmat.NA.pad <- rbind( 0, cbind( 0, hmat.NA, 0 ), 0 )
hmf.NA <- hmf <- 0*hmat
for (dx in c(-1,0,1)) { 
    for (dy in c(-1,0,1)) {
        hmf <- hmf + hmat.pad[dx+1+(1:nrow(hmf)),dy+1+(1:ncol(hmf))]*w[2+dx,2+dy]
        add.these <- hmat.NA.pad[dx+1+(1:nrow(hmf)),dy+1+(1:ncol(hmf))]
        hmf.NA <- hmf.NA + ifelse(is.na(add.these),0,add.these)*w[2+dx,2+dy]
    }
}

testthat::expect_equivalent(as.numeric(hmf),values(hf1))
testthat::expect_equivalent(as.numeric(hmf.NA),values(hf1.NA))


# check simple smoothing matrix construction

M <- migration_matrix( habitat, kern=function(x) { ifelse(x>0,1/16,1/2) }, sigma=1, radius=1, normalize=NULL )
matrix(M[3,],nrow=5)
testthat::expect_equal( hmf, matrix(M%*%as.numeric(hmat),nrow=nrow(hmf)) )
testthat::expect_equal( hmf.NA, matrix(M%*%na_to_zero(as.numeric(hmat.NA)),nrow=nrow(hmf)) )


########
testthat::context("Absorbing external, not internal, boundaries")

wm.fun <- function (x, na.rm) { weighted.mean(x,w=w,na.rm=na.rm) }
hf2 <- focal( habitat, w=(w>0), pad=TRUE, fun=wm.fun, na.rm=TRUE, padValue=0 )
hf2.NA <- focal( habitat.NA, w=(w>0), pad=TRUE, fun=wm.fun, na.rm=TRUE, padValue=0 )

pad.extent <- extent(habitat)+c(-1,1,-1,1)*res(habitat) # extends by one cell
pad.habitat <- extend(habitat,pad.extent,value=0)
M <- migration_matrix( pad.habitat, kern=function(x) { ifelse(x>0,1/16,1/2) }, sigma=1, radius=1, normalize=1 )
M <- subset_migration( M, old=pad.habitat, new=habitat )
matrix(M[3,],nrow=nrow(habitat))

testthat::expect_equal( values(hf2), as.numeric(M%*%values(habitat)) )

pad.habitat.NA <- extend(habitat.NA,pad.extent,value=0)
M.NA <- migration_matrix( pad.habitat.NA, kern=function(x) { ifelse(x>0,1/16,1/2) }, sigma=1, radius=1, normalize=1 )
M.NA <- subset_migration( M.NA, old=pad.habitat.NA, new=habitat.NA  )

testthat::expect_equal( values(hf2.NA)[!is.na(values(habitat.NA))], as.numeric(M.NA%*%(values(habitat.NA)[!is.na(values(habitat.NA))])) )

########
testthat::context("Non-boundaries")


wm.fun <- function (x, na.rm) { weighted.mean(x,w=w,na.rm=na.rm) }
hf3 <- focal( habitat, w=(w>0), pad=TRUE, fun=wm.fun, na.rm=TRUE, padValue=NA )
hf3.NA <- focal( habitat.NA, w=(w>0), pad=TRUE, fun=wm.fun, na.rm=TRUE, padValue=NA )

M <- migration_matrix( habitat, kern=function(x) { ifelse(x>0,1/16,1/2) }, sigma=1, radius=1, normalize=1 )
matrix(M[3,],nrow=nrow(habitat))
testthat::expect_equal( values(hf3), as.numeric(M%*%values(habitat)) )

M.NA <- migration_matrix( habitat.NA, kern=function(x) { ifelse(x>0,1/16,1/2) }, sigma=1, radius=1, normalize=1 )
testthat::expect_equal( values(hf3.NA)[!is.na(values(habitat.NA))], as.numeric(M.NA%*%(values(habitat.NA)[!is.na(values(habitat.NA))])) )



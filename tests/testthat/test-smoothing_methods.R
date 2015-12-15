testthat::context("Agreement between migration_matrix and migrate.")

# setup

this.migration <- migration(
                    kern = "gaussian",
                    sigma = 100,
                    radius = 400,
                    normalize = NULL
             )

habitat <- raster(xmn=-1000/2, xmx=1000/2, ymn=-1000/2, ymx=1000/2,
      resolution=100, # m
      crs="+proj=utm +zone=11 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
values(habitat) <- sample( 100*c(1,2,5,10), length(habitat), replace=TRUE )

migrate.habitat <- migrate(habitat,this.migration )

M <- migration_matrix( habitat, this.migration )
M.habitat <- habitat
nonzeros <- !is.na(values(M.habitat))
values(M.habitat)[nonzeros] <- ( M %*% values(M.habitat)[nonzeros] )@x

## we DO NOT expect these to be the same at because normalization 
## and treatment of the boundary is different.

# ( M.habitat-migrate.habitat )

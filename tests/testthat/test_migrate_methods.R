library(raster)
library(Matrix)

context("testing normalization of migrate()")

set.seed(23)
habitat <- random_habitat(diam=2000, sigma=200,
                  randfun=function(n) pmin(20, (rcauchy(n))) )

###
context("without accessible, uninhabitable locations")

pop <- population( 
                  habitat = habitat,
                  accessible = (!is.na(values(habitat))),
                  habitable = (!is.na(values(habitat))),
                  genotypes = c("aa")
             )
pop$N[,"aa"] <- rpois(sum(pop$habitable),values(pop$habitat)[pop$habitable])

###
context("with n.weights=1")

migr <- migration(
                kern = "gaussian",
                sigma = 200,
                radius = 300,
                normalize = 1
         )

migr <- setup_migration( migr, pop )

x <- rbinom(nrow(pop$N),1,0.1)
x <- x/sum(x)
Mx <- migrate( x, migr )

test_that( "basic properties of migration matrix with n.weights=1", {
          expect_equal( Matrix::rowSums(migr$M) , rep(1.0,sum(pop$accessible)) )
          expect_equal( sum(x), 1.0 )
          expect_equal( sum(Mx), 1.0 )
     } )

###
context("with other n.weights")

w.migr <- setup_migration( migration(
                kern = "gaussian",
                sigma = 200,
                radius = 300,
                normalize = 1,
                n.weights=0.1 * 0.9^(1:10)
         ), pop )


Mx <- migrate( x, w.migr )

test_that( "basic properties of migration matrix with n.weights not 1", {
          expect_equal( Matrix::rowSums(migr$M) , rep(1.0,sum(pop$accessible)) )
          expect_equal( sum(Mx), 1.0 )
     } )


###
context("migrate with n.weights matches weighted matrix sum")

Mn <- Diagonal(nrow(w.migr$M),x=1.0)
MM <- (1-sum(w.migr$n.weights))*Mn
for (n in seq_along( w.migr$n.weights ) ) {
    Mn <- w.migr$M %*% Mn
    MM <- MM + w.migr$n.weights[n] * Mn
}

test_that( "migrate with n.weights matches weighted sum matrix", {
          expect_equal( Matrix::rowSums(MM) , rep(1.0,sum(pop$accessible)) )
          expect_equal( as.numeric(Matrix::crossprod(MM,x)), Mx )
     } )



###
context("with accessible, uninhabitable locations")

pop <- population( 
                  habitat = habitat,
                  accessible = (!is.na(values(habitat))),
                  habitable = (!is.na(values(habitat)) & values(habitat)>2),
                  genotypes = c("aa")
             )
pop$N[,"aa"] <- rpois(sum(pop$habitable),values(pop$habitat)[pop$habitable])

###
context("with n.weights=1")

migr <- migration(
                kern = "gaussian",
                sigma = 200,
                radius = 300,
                normalize = 1
         )

migr <- setup_migration( migr, pop )

x <- rep(1/nrow(pop$N),nrow(pop$N))
Mx <- migrate( x, migr )

test_that( "basic properties of migration matrix with n.weights=1", {
          expect_equal( Matrix::rowSums(migr$M) , rep(1.0,sum(pop$accessible)) )
          expect_true( sum(Mx) <= 1.0 )
     } )

###
context("with other n.weights")

w.migr <- setup_migration( migration(
                kern = "gaussian",
                sigma = 200,
                radius = 300,
                normalize = 1,
                n.weights=0.1 * 0.9^(1:10)
         ), pop )


x <- rep(1/nrow(pop$N),nrow(pop$N))
Mx <- migrate( x, w.migr )

test_that( "basic properties of migration matrix with n.weights not 1", {
          expect_equal( Matrix::rowSums(migr$M) , rep(1.0,sum(pop$accessible)) )
          expect_true( sum(Mx) <= 1.0 )
     } )


###
context("migrate with n.weights matches weighted matrix sum")

Mn <- Diagonal(nrow(w.migr$M),x=1.0)
MM <- (1-sum(w.migr$n.weights))*Mn
for (n in seq_along( w.migr$n.weights ) ) {
    Mn <- w.migr$M %*% Mn
    MM <- MM + w.migr$n.weights[n] * Mn
}

tmp.x <- numeric(nrow(MM))
tmp.x[ w.migr$habitable.inds ] <- x
MMx <- Matrix::crossprod(MM, tmp.x)

test_that( "migrate with n.weights matches weighted sum matrix", {
          expect_equal( Matrix::rowSums(MM) , rep(1.0,sum(pop$accessible)) )
          expect_equal( as.numeric(MMx)[w.migr$habitable.inds], Mx )
     } )


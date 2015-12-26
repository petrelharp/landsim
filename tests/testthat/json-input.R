context("Input and output of JSON.")
require(jsonlite)
library(Matrix)
library(raster)


#########
context("population config")

pop.list <- list(
                 habitat = system.file("inst/extdata/test_raster.gri",package="landsim"),
                 inaccessible.value = "NA",
                 uninhabitable.value = 0.0,
                 genotypes = c( "aa", "aA", "AA" ),
                 N=1
             )
pop.json <- toJSON(pop.list,pretty=TRUE)

pop.from.json <- make_population(json=pop.json)
pop.from.list <- do.call( make_population, pop.list )

testthat::expect_equal( pop.from.json, pop.from.list )


#########
context("migration config")

migr.list <- list(
                 kern = "gaussian",
                 sigma = 100,
                 radius = 1000,
                 normalize = 1.0
             )
migr.json <- toJSON(migr.list,prett=TRUE)

migr.from.json <- make_migration( json=migr.json )
migr.from.list <- do.call( make_migration, migr.list )

testthat::expect_equal( migr.from.json, migr.from.list )



#########
context("demography config")

population <- pop.from.json

vital.list <- list(
             r0 = 0.01,  # one in ten seeds will germinate at low densities
             s = 1.5,    # multiplicative selective benefit of the A allele
             carrying.capacity = 10*habitat,
             competition = migration(
                                     kern="gaussian",
                                     sigma=100,
                                     radius=300,
                                     normalize=1,
                                     do.M=TRUE,
                                     population=population
                                 ),
             value = function (N, ...) {
                out <- r0 / ( 1 + migrate(competition,x=rowSums(N))/carrying.capacity )
                return( cbind( aa=out, aA=s*out, AA=s^2*out ) )
             }
         )


vital.json <- '
{
    "value" : { "function" : 
              " out <- r0 / ( 1 + migrate(competition,x=rowSums(N))/carrying.capacity ) ;
                return( cbind( aa=out, aA=s*out, AA=s^2*out ) ) " 
        },
    "r0" : 0.01,
    "s" : 1.5,
    "carrying.capacity" : { "R" : "10*habitat" },
    "competition" : { "migration" : { 
            "kern" : "gaussian",
            "sigma" : 100,
            "radius" : 300,
            "normalize" : 1.0
        } }
}
'

## this mostly just checks there aren't errors
## rather than checking the output

dir.create("testcache",showWarnings=FALSE)

require(raster)

##########
context("Set up random habitat")


seedish <- floor(1e6*runif(1))
# reproducible landscape
set.seed(42)

habitat <- random_habitat(diam=2000,res=200)
carrying.capacity <- values(habitat)[(!is.na(values(habitat)) & values(habitat)>0)]

context("Set up population")

pop <- population( 
                  habitat = habitat,
                  accessible = !is.na(values(habitat)),
                  habitable = (!is.na(values(habitat)) & values(habitat)>0),
                  genotypes = c("aa","aA","AA"),
                  carrying.capacity = carrying.capacity,
                  N = cbind( aa=rpois(length(carrying.capacity),carrying.capacity),
                             aA=0, 
                             AA=0 ),
             )

test_that("make pop", {
        expect_equal_to_reference(pop,"testcache/saved_pop_test.Rd") 
     } )

##########
context("Set up demography")

migr <- migration( kern="gaussian", sigma=100, radius=400, normalize=1 )
demog <- demography(
        prob.seed = 0.2,
        fecundity = 100,
        carrying.capacity=10,
        prob.germination = vital( 
                 function (N, ...) {
                     out <- r0 / ( 1 + rowSums(N)/carrying.capacity )
                     return( cbind( aa=out, aA=s*out, AA=s^2*out ) )
                 },
                 r0 = 0.01,  # one in ten seeds will germinate at low densities
                 s = 1.5     # multiplicative selective benefit of the A allele
             ),
        prob.survival = 0.9,
        pollen.migration = migr,
        seed.migration = migr,
        genotypes = c("aa","aA","AA"),
    )




##########
context("set up demography for population")

demog <- setup_demography( demog, pop )


test_that("make demog", {
        expect_equal_to_reference(demog,"testcache/saved_demog_test.Rd") 
     } )


##########
context("check for population stability")

growths <- intrinsic_growth(pop,demog)

test_that("pop growth", {
        expect_equal( as.numeric(cor(growths)), rep(1.0,9) )
        expect_true( all( ( growths[,"aa"] > 1.0 ) & ( growths[,"aa"] < 1.2 ) ) )
        expect_true( all( ( growths[,"AA"] > 1.2 ) & ( growths[,"AA"] < 1.5 ) ) )
     } )

#########
context("initialize population values")

mutloc <- list( cell.number=sample.int(nhabitable(pop),5) )
mutloc$center <- SpatialPoints(xyFromCell(pop$habitat,which(pop$habitable)[mutloc$cell.number]))
pop <- set_N( pop, i=which(pop$habitable), j="aA", 
         value=ifelse(1:nrow(pop$N)%in%mutloc$cell.number,10,0) )

test_that( "using set_N", {
        expect_equal( c(0,10)[1+((1:nrow(pop$N)) %in% mutloc$cell.number)], pop$N[,2] )
        expect_equal( pop$N[mutloc$cell.number,"aA"], rep(10,length(mutloc$cell.number)) )
     } )

##########
context("run a simulation")

sim <- simulate_pop( pop, demog, times=seq(0,100,length.out=11),
                 summaries= list( totals=function(N){colSums(N)} ),
                 stop.fun=function(N){ sum(N)==0 } )

simbr <- sim_to_brick( sim, pop )

sim <- extend_simulation( sim, pop, demog, times=seq(sim$t,200,length.out=11),
                 summaries= list( totals=function(N){colSums(N)} ) )


test_that("run sims", {
        expect_equal( sapply(simbr,dim), 
             structure(c(10L, 10L, 11L, 10L, 10L, 11L, 10L, 10L, 11L), .Dim = c(3L, 3L), .Dimnames = list(NULL, c("aa", "aA", "AA"))) )
        expect_equal_to_reference(sim,"testcache/saved_sim_test.Rd") 
     } )


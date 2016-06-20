
##########
context("Test writing out matrices in triplet form")


set.seed(42)

habitat <- random_habitat(diam=2000,res=200)
carrying.capacity <- values(habitat)[(!is.na(values(habitat)) & values(habitat)>0)]


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


migr <- migration( kern="gaussian", sigma=100, radius=400, normalize=1 )

M <- migration_matrix( population=pop, migration=migr )
ref.M <- read.csv("write_matrix_ref.csv",header=TRUE)

write_triplets(M=M, file="write_matrix_ref2.csv", sep=',')
write.M <- read.csv("write_matrix_ref2.csv",header=TRUE)

expect_equal(ref.M, write.M)


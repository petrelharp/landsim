context("Construction of mating tensor.")

haploid.gt <- c("a","A")
diploid.gt <- c("aa","aA","AA")

diploid.mt <- array( c(
        # aa offspring
        1,   1/2,   0, # aa
        1/2, 1/4,   0, # aA
        0,     0,   0, # AA
        # aA offspring
        0,   1/2,   1, # aa
        1/2, 1/2, 1/2, # aA
        1,   1/2,   0, # AA
        # AA offspring
        0,     0,   0, # aa
        0,   1/4, 1/2, # aA
        0,   1/2,   1  # AA
    ), dim=c(3,3,3) )
dimnames(diploid.mt) <- list( diploid.gt, diploid.gt, diploid.gt )

haplodip.mt <- array( c(
        # aa offspring
        # a    A
        1,     0,  # aa
        1/2,   0,  # aA
        0,     0,  # AA
        # aA offspring
        # a    A
        0,     1,  # aa
        1/2, 1/2,  # aA
        1,     0,  # AA
        # AA offspring
        # a    A
        0,     0,  # aa
        0,   1/2,  # aA
        0,     1   # AA
    ), dim=c(2,3,3) )
dimnames(haplodip.mt) <- list( haploid.gt, diploid.gt, diploid.gt )

expect_identical( diploid.mt, mating_tensor(diploid.gt) )
expect_identical( diploid.mt, mating_tensor(diploid.gt,diploid.gt,diploid.gt) )

expect_identical( haplodip.mt, mating_tensor(haploid.gt, diploid.gt, diploid.gt) )
expect_identical( aperm(haplodip.mt,c(2,1,3)), mating_tensor(diploid.gt, haploid.gt) )

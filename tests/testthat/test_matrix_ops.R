context("Geometric power")

M <- matrix( runif( 25 ), nrow=5 )
diag(M) <- 5-rowSums(M)
M <- M/rowSums(M)

p <- 0.2
eps <- 1e-8
pM <- (1-p) * MASS::ginv(diag(nrow(M))-p*M)

gM <- geo_power(M,p=p,eps=eps)
expect_true( all( abs(gM-pM) < eps ) )

eps <- 1e-12
gM <- geo_power(M,p=p,eps=eps)
expect_true( all( abs(gM-pM) < eps ) )

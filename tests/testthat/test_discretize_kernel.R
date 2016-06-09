
context("testing kernel discretization")

# Discretization should be finer than sigma:
# it seems to be enough so that res/fact = sigma/3, i.e., if fact = 3*res/sigma

### first test at a not-terribly-small value of sigma

sigma <- 0.1
res <- 1
f <- function (fact) { discretize_kernel('gaussian', res=c(res,res), radius=3*sigma, sigma=sigma, fact=fact) }
fact.vec <- c(5,10,15,20,25,30)
names(fact.vec) <- fact.vec
f.list <- lapply( fact.vec, f )

# check at one value
x.ends <- seq(-0.5,0.5,length.out=51)
xx <- x.ends[-1]-diff(x.ends)/2
dx <- mean(diff(xx))

# offset by (1,0)
right.grid <- expand.grid( 
                         x0=xx, y0=xx,
                         x1=1+xx, y1=0+xx
                     )
right.grid$z <- dx^2 * (dnorm(right.grid$x1-right.grid$x0,sd=sigma)*dnorm(right.grid$y1-right.grid$y0,sd=sigma))
right.check <- sum(right.grid$z)/length(xx)^2

# simulation check
sim <- data.frame( x0=runif(1e6)-1/2, y0=runif(1e6)-1/2 )
sim$x1 <- sim$x0 + rnorm(1e6,sd=sigma)
sim$y1 <- sim$y0 + rnorm(1e6,sd=sigma)
sim.counts <- table( cut(sim$x1,((-1:2)-1/2)*res), 
                     cut(sim$y1,((-1:2)-1/2)*res) )

# offset by (1,1)
diag.grid <- expand.grid( 
                         x0=xx, y0=xx,
                         x1=1+xx, y1=1+xx
                     )
diag.grid$z <- dx^2 * (dnorm(diag.grid$x1-diag.grid$x0,sd=sigma)*dnorm(diag.grid$y1-diag.grid$y0,sd=sigma))
diag.check <- sum(diag.grid$z)/length(xx)^2

# compare
results <- rbind(
    'right'=c( sapply( f.list, function (g) g(1/sigma)*res^2/sigma^2 ), check=right.check, sim=mean(sim.counts[c(2,4,6,8)])/1e6 ),
    'diag'=c( sapply( f.list, function (g) g(sqrt(2)/sigma)*res^2/sigma^2 ), check=diag.check, sim=mean(sim.counts[c(1,3,7,9)])/1e6 )
    )

expect_true( 
            all( abs(results[,"30"]-results[,"sim"]) < 5*sqrt(1e6*results[,"sim"]) ) 
            && all( abs(results[,"30"]-results[,"check"]) < 5*sqrt(1e6*results[,"check"]) ) )


### now test at a smaller value of sigma

sigma <- 0.01
res <- 1
f <- function (fact) { discretize_kernel('gaussian', res=c(res,res), radius=3*sigma, sigma=sigma, fact=fact) }
fact.vec <- (1:3)*res/sigma
names(fact.vec) <- fact.vec
f.list <- lapply( fact.vec, f )

# simulation check
sim <- data.frame( x0=runif(1e6)-1/2, y0=runif(1e6)-1/2 )
sim$x1 <- sim$x0 + rnorm(1e6,sd=sigma)
sim$y1 <- sim$y0 + rnorm(1e6,sd=sigma)
sim.counts <- table( cut(sim$x1,((-1:2)-1/2)*res), 
                     cut(sim$y1,((-1:2)-1/2)*res) )

# compare
results <- rbind(
    'right'=c( sapply( f.list, function (g) g(1/sigma)*res^2/sigma^2 ), sim=mean(sim.counts[c(2,4,6,8)])/1e6 ),
    'diag'=c( sapply( f.list, function (g) g(sqrt(2)/sigma)*res^2/sigma^2 ), sim=mean(sim.counts[c(1,3,7,9)])/1e6 )
    )

expect_true( 
            all( abs(results[,names(fact.vec)[length(fact.vec)]]-results[,"sim"]) < 5*sqrt(1e6*results[,"sim"]) ) 
        )


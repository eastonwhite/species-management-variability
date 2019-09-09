# --Dispersal by stochastic diffusion--
# A version of DiffSt.r that allows for overdispersion in the p[y->x]
# probabilities and also for modeling the dispersal via a multinomial
# distribution (equivalent to DiffSt.r). Overdispersion is added to 
# the multinomial by compounding with a Dirichlet distribution. This
# compound distribution is known as a multivariate Polya distribution.
#
# Updated Aug 2015 with left-right boundaries
#
# Requires Diff.r
#
library(gtools) #for rdirichlet

DiffSt_Polya <- function(N,tt,D,s,theta_D=0) {
  tol = .Machine$double.eps ^ 0.5
  if ( sum(N) < tol ) return(N) #If landscape is extinct, return N, which will be all zeros
  np <- length(N) #number of patches
  occupied <- which(N >= 0.5)
  left <- min(occupied) #leftmost extent of population
  right <- max(occupied) #rightmost extent
  Nnew <- rep(0,np)
  for ( i in left:right ) {
    if ( N[i] == 0 ) next   #for sparse populations
    ls <- rep(0,np)
    ls[i] <- 1
    kbinompois <- Diff(ls,tt,D)          #dispersal kernel
    tmp <- rdirichlet(1,s*kbinompois)   #overdispersion
    sim <- rmultinom(1,N[i],tmp)        #multinomial
    Nnew <- Nnew + sim
  }
  return(as.vector(Nnew))
}


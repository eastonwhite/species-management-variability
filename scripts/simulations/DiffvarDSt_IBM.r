#----Dispersal by negative binomial stochastic diffusion. IBM.----------------------
# The stochastic version models the dispersal events out of each patch,
# by finding the time to the next dispersal event (Poisson process).
# Between or within individual heterogeneity in dispersal rate is modelled as a
# gamma distribution (D is the mean, k is the shape parameter).
# The scale of variation is specified as "between" or "within" individual.
# For DI dispersal, we could model each individual separately. However, we retain
# the patch structure to allow for building in DD interactionsn later.
# The model is based on the Poisson IBM, with only a couple of lines different.
#
# Let Im(x) indicate the presence (1) or absence (0) of individual m at location x. 
# P[ Im(x) != Im(x)-1 ] = exp(-Dm*holes*t)
# D is gamma with shape parameter k.
# Im(x+j)=Im(x+j)+1 | Im(x)=Im(x)-1, P[j] = 0.5, j in (-1,1), for 2 hole patches
# i.e. random direction from two hole patches. This is equivalent to modeling
# the event times for holes individually.
#
#

DiffvarDSt_IBM <- function(N,tt,D,k,vscale="between") {
#  if ( sum(N %% trunc(N)) > 0 ) stop("StDiff_IBM needs integer N") #check for fractional individuals
# Above line gives NaN for N=0. Need to fix that.
  if (length(D) > 1) stop("StDiff_IBM accepts only scalar D")
  np <- length(N) #number of patches
  nI <- sum(N) #number of individuals
#check for empty landscape
  if ( nI==0 ) {
    warning("There are no individuals to disperse, returning empty landscape") 
    return(N)
  }

# Initialize individuals
  Ipos <- NA*(1:nI)
  m <- 1
  for (i in 1:np) { #there's probably a faster way than reading one by one
    if ( N[i] > 0 ) {
      for (j in 1:N[i]) {
        Ipos[m] <- i
        m <- m + 1
      }
    }
  }
  holes <- rep(2,nI) - 1 * ((Ipos == 1) | (Ipos == np))
  Dm <- rgamma( nI, shape = k, scale = D / k ) #random D for individuals
  evT <- rexp(nI) / (Dm * holes) #Initial times for dispersal events.

# Initialize clock and scheduler
  next_I <- which.min(evT) #schedule next event (which individual is it?)
  t <- min(evT) #set clock to time of first event

# Disperse individuals
  while (t < tt) {
  # Carry out the dispersal event
    if ( runif(1) < 0.5 ) j <- -1 else j <- 1 #random direction
    if (Ipos[next_I] == 1) j <- 1 #right if first patch
    if (Ipos[next_I] == np) j <- -1 #left if last patch
    Ipos[next_I] <-  Ipos[next_I] + j #move to neighbouring patch
  # Update holes
    holes[next_I] <- 2 - 1 * ((Ipos[next_I] == 1) | (Ipos[next_I] == np))
  #
  # Update D (for within-individual gamma variation in D).
    if (vscale=="within") {
      Dm[next_I] <- rgamma( 1, shape = k, scale = D / k ) #random D for individuals
    }
  # Update event time
    evT[next_I] <- t + rexp(1) / (Dm[next_I] * holes[next_I])
  # Reset clock and scheduler
    next_I <- which.min(evT)
    t <- min(evT)
  }

# Recollate the landscape
  N <- N*0
  for (m in 1:nI) {
    N[Ipos[m]] = N[Ipos[m]] + 1    
  }
  N #return the new landscape
}

# --Dispersal by stochastic diffusion--
# The stochastic version models the dispersal events out of each patch
# by finding the time to the next dispersal event (Poisson process).
# P[ N(x) != N(x)-1 ] = exp(-D*N(x)*holes*t), N(x) >= 0
# N(x+j)=N(x+j)+1 | N(x)=N(x)-1, P[j] = 0.5, j in (-1,1), for 2 hole patches
# i.e. random direction from two hole patches. This is equivalent to modeling
# the event times for holes individually.
#
DiffSt <- function(N,tt,D,theta_D=0) {
  np <- length(N) #number of patches
  if (length(D)==1) D <- D + 0*(1:np) #Expand D to a vector if it isn't already
  #fN <- D * holes * y^theta_D #Density-dependent emigration rate, theta_D = 0 is DI        
  holes <- c(1,rep(2,(np-2)),1)     #Number of holes
  evT <- rexp(np) / (D * N * holes) #Initial times for dispersal events. n.b. INF for N=0.

# Initialize clock and scheduler
  next_x <- which.min(evT) #schedule next event (which patch is it in?)
  t <- min(evT) #set clock to time of first event

# Disperse individuals
  while (t < tt) {
  # Carry out the dispersal event
    if (N[next_x] < 1) {  #Allows the last bit to be dispersed if we are not modeling discrete individuals
      indsz <- N[next_x]
    } else {
      indsz <- 1
    }
    N[next_x] <- N[next_x] - indsz            #emigrate from patch
    if ( runif(1) < 0.5 ) j <- -1 else j <- 1 #random direction
    if (next_x == 1) j <- 1                   #right if first patch
    if (next_x == np) j <- -1                 #left if last patch
    N[next_x+j] <-  N[next_x+j] + indsz       #immigrate to new patch
  # Update event times
    evT[next_x] <- t + rexp(1) / (D[next_x] * N[next_x] * holes[next_x])
    evT[next_x+j] <- t + rexp(1) / (D[next_x+j]*N[next_x+j]*holes[next_x+j])
  # Reset clock and scheduler
    next_x <- which.min(evT)
    t <- min(evT)
  }
  N #return the new landscape
}

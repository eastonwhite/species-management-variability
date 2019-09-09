#Solve the diffusion model for between-individual variation in D.
#The inner integral solves the basic diffusion model, a linear ODE system, 
#using linear algebra. The matrix A specifies the system.
#The outer integral uses Eulers-McLaurin trapezoid geometry on P[D=d] - i.e. 
#equal increments on the probability scale. This effectively gives us an
#adaptive step size - we take a smaller step size on the D scale where P[D=d]
#is high. This version adds an automatic stopping rule to LAint3.
#
#ls:      Initial conditions for the landscape
#tt:      Total time
#Dbar, k: Parameters of the gamma distribution (mean D, shape).
#ptol:    Determines the limits of the integration on the p-scale. This is a
#         crude but workable way to deal with both the singularity at D=0 and
#         the infinite upper limit.
#eps:     Relative precision of the solution.
#nmax:    Determines max increments (20 -> 2^20, about 1 million increments)
#
DiffvarD_LAint <- function(ls,tt,Dbar,k,ptol=1e-6,eps=1e-5,nmax=20){
  np <- length(ls) #number of patches
  A <- matrix(c(-1,1,rep(0,np-2),rep(c(1,-2,1,rep(0,np-2)),np-2),1,-1),np,np)
  eig <- eigen(A);
  const <- solve(eig$vectors,ls); #ls is initial conditions

# Inner integral: integrate over time for one value of D
  diffuse <- function( D ) {
    eig$vectors %*% (const*exp(tt*eig$values*D))
  }
# Outer integral: integrate over D using an extended trapezoid rule
# Following Numerical Recipes in C, trapzd.c and qtrap.c.
  maxi <- which(ls == max(ls))[1] #The most abundant patch. Take the first if multiple.
  oldlsint <- -1.0e30
  Dl <- qgamma(ptol, shape = k, scale = Dbar / k )
  Du <- qgamma(1-ptol, shape = k, scale = Dbar / k )
  Dhalf <- qgamma(0.5, shape = k, scale = Dbar / k )
# One whole increment first and the halfway point
  lssum <- 0.5*(diffuse(Dl)+diffuse(Du))
  lssum <- lssum + diffuse(Dhalf)
# Progressively smaller increments (halving size each time)
  for (n in 2:nmax) {
    totincs <- 2^n                #The total number of increments (old + new)
    newpts <- 2^(n-1)             #The number of new points
    nwdel <- (1-2*ptol)/newpts    #New delta
    p_evals <- seq(ptol+0.5*nwdel,1-ptol-0.5*nwdel,length.out=newpts)
    D_evals <- qgamma(p_evals, shape = k, scale = Dbar / k )
    lsmids <- ls*0                #Initialize landscape to zero
    for (D in D_evals) {
      lsmids <- lsmids + diffuse(D)
    }
    lssum <- lssum + lsmids       #Add the new points to previous ones
    pinc <- (1-2*ptol)/totincs    #Current increment size on probability scale
    lsint <- lssum * pinc         #Current approximation for integrated landscape
  # Check for convergence (based on the most abundant patch)
    if (n > 5){                   #Avoid spurious early convergence
      if ( abs(lsint[maxi]-oldlsint[maxi]) < eps*abs(oldlsint[maxi]) ) {
        return(lsint)
      }
    }
    oldlsint <- lsint
  }
  stop("Too many steps in integration")
}


#Eulers-McLaurin trapezoid geometry on p(D) - i.e. equal increments on p(D) scale.
#This is more accurate than LAint1 or 2 for a given number of function evaluations, and
#therefore faster because fewer evaluations are needed.
DiffvarD_LAint3 <- function(ls,tt,Dbar,k,incs,ptol){
  np <- length(ls) #number of patches
  A <- matrix(c(-1,1,rep(0,np-2),rep(c(1,-2,1,rep(0,np-2)),np-2),1,-1),np,np)
  eig <- eigen(A);
  const <- solve(eig$vectors,ls); #ls is initial conditions

# Integrate over time for one value of D
  diffuse <- function( D ) {
    eig$vectors %*% (const*exp(tt*eig$values*D))
  }
# Outer integration. Eulers-like, but a built in adaptive step (increments are on probability scale)
  pinc <- (1-2*ptol)/incs #nb inc on the p scale is the probability that D is in the increment on the D scale
# Endpoints first
  Dl <- qgamma(ptol, shape = k, scale = Dbar / k )
  Du <- qgamma(1-ptol, shape = k, scale = Dbar / k )
  lsends <- (pinc/2)*(diffuse(Dl)+diffuse(Du))
# Interior points
  lsmids <- ls*0
  p_evals <- seq(ptol+pinc,1-ptol-pinc,length.out=incs-1)
  D_evals <- qgamma(p_evals, shape = k, scale = Dbar / k )
  for (D in D_evals) {
    lsmids <- lsmids + diffuse(D)
  }
  lsmids <- lsmids*pinc
  return(lsends+lsmids)
}

#Eulers-McLaurin trapezoid geometry on D.
#Equal increments on D scale
#Weighting by p(D) fixes the singularity issue.
#This works because the approximation is diffuse(D)*f(D)*deltaD, and 
#f(D)*deltaD is the approximation for f(D)dD, the area under the density 
#function f(D), for which an alternative and more accurate (exact)
#approximation is p(D), the probability that D is in the interval spanned by 
#deltaD. Intuitively, it is the sumD( diffuse(d) * P[D=d] ).
#
DiffvarD_LAint2 <- function(ls,tt,Dbar,k,incs,ptol){
  np <- length(ls) #number of patches
  A <- matrix(c(-1,1,rep(0,np-2),rep(c(1,-2,1,rep(0,np-2)),np-2),1,-1),np,np)
  eig <- eigen(A);
  const <- solve(eig$vectors,ls); #ls is initial conditions

# Define the outer integral as solution to ODE dN/dD = diffusion(D*p[D=D])
  diffuse <- function( D ) { 
    eigvals <- eig$values*D
    eig$vectors %*% (const*exp(tt*eigvals))
  }
# Eulers
  Dl <- qgamma(ptol, shape = k, scale = Dbar / k )
  Du <- qgamma(1-ptol, shape = k, scale = Dbar / k )
  inc <- (Du-Dl)/incs
# First we do the initial increment, which is half sized
  cdf <- pgamma(c(Dl,Dl+0.5*inc), shape = k, scale = Dbar / k )
  p_in_int <- cdf[2]-cdf[1] #probability D is in the interval spanned by the increment
  lsnew <- diffuse(Dl)*p_in_int #This also initializes lsnew
# Then we do the intermediate increments
  for (D in seq(Dl+inc,Du-inc,length.out=incs-1)) {
    cdf <- pgamma(c(D-0.5*inc,D+0.5*inc), shape = k, scale = Dbar / k )
    p_in_int <- cdf[2]-cdf[1] #probability D is in the interval spanned by the increment
    lsnew <- lsnew + diffuse(D)*p_in_int
  }
# Finally we add the last increment, which is half sized
  cdf <- pgamma(c(Du-0.5*inc,Du), shape = k, scale = Dbar / k )
  p_in_int <- cdf[2]-cdf[1] #probability D is in the interval spanned by the increment
  lsnew <- lsnew + diffuse(Du)*p_in_int
  return(lsnew)
}


#Eulers-McLaurin trapezoid geometry on D.
#Equal increments on D scale.
#Proper equation for the expected value, i.e. weighted by f(D), not p(D)
#This blows up because of the singularity in dgamma (dgamma=Inf at D=0, k<1)
#
DiffvarD_LAint1 <- function(ls,tt,Dbar,k,incs,ptol){
  np <- length(ls) #number of patches
  A <- matrix(c(-1,1,rep(0,np-2),rep(c(1,-2,1,rep(0,np-2)),np-2),1,-1),np,np)
  eig <- eigen(A);
  const <- solve(eig$vectors,ls); #ls is initial conditions

# Inner integral: integrate over time for one value of D
  diffuse <- function( D ) { 
    eigvals <- eig$values*D
    fin <- eig$vectors %*% (const*exp(tt*eigvals))
    fin * (fin > 0) #Catch small negative values and set to zero
  }
# Outer integration. Trapezoidal rule.
# Endpoints first
  Dl <- qgamma(ptol, shape = k, scale = Dbar / k )
  Du <- qgamma(1-ptol, shape = k, scale = Dbar / k )
  inc <- (Du-Dl)/incs
  lsends <- exp(log(diffuse(Dl))+dgamma(Dl,shape=k,scale=Dbar/k,log=TRUE))
  lsends <- lsends + exp(log(diffuse(Du))+dgamma(Du,shape=k,scale=Dbar/k,log=TRUE))
  lsends <- lsends * (inc/2)
# Interior points
  lsmids <- ls*0
  D_evals <- seq(Dl+inc,Du-inc,length.out=incs-1)
  for (D in D_evals) {
    lsmids <- lsmids + exp(log(diffuse(D))+dgamma(D,shape=k,scale=Dbar/k,log=TRUE))
  }
  lsmids <- lsmids * inc
  return(lsends+lsmids)
}

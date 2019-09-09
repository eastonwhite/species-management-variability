# Probability mass functions (pmf) for the Ricker model

#----Poisson pmf for Ricker model
#
Ricker_dpois <- function(Ntp1,Nt,R,alpha,log=FALSE){
  lambda <- Nt*R*exp(-alpha*Nt)
  if (log==FALSE) {
    return( dpois(Ntp1,lambda) )
  } #else
  return( dpois(Ntp1,lambda,log=TRUE) )
}

#----Negative binomial (demographic) pmf for Ricker model
#
Ricker_dnbinom_d <- function(Ntp1,Nt,R,alpha,k,log=FALSE){
  mu <- Nt*R*exp(-alpha*Nt)
  if (log==FALSE) {
    return( dnbinom(Ntp1,size=k*Nt+(Nt==0),mu=mu) )
    #Add 1 to size parameter of NB when Nt=0 to avoid NaN's.
  } #else
  return( dnbinom(Ntp1,size=k*Nt+(Nt==0),mu=mu,log=TRUE) )
}

#----Negative binomial (environmental) pmf for Ricker model
#
Ricker_dnbinom_e <- function(Ntp1,Nt,R,alpha,k,log=FALSE){
  mu <- Nt*R*exp(-alpha*Nt)
  if (log==FALSE) {
    return( dnbinom(Ntp1,size=k,mu=mu) )
  } #else
  return( dnbinom(Ntp1,size=k,mu=mu,log=TRUE) )
}

#----Negative binomial-gamma pmf for Ricker model
#
#This version requires a standard number of increments and has no automatic
#resolution setting. Because of this, the whole vector of Ntp1 can be done at
#once and it is faster. Use version 2 (slower) if higher accuracy is needed.
#The integral uses Eulers-McLaurin trapezoid geometry on P[R=r] - i.e. 
#equal increments on the probability scale. This effectively gives us an
#adaptive step size - we take a smaller step size on the R scale where P[R=r]
#is high. This follows Numerical Recipes in C, trapzd.c and qtrap.c.
#The ptol and nmax values are quite course but sufficient for model fitting.
#
Ricker_dnbinomgamma <- function(Ntp1,Nt,R,alpha,kD,kE,log=FALSE,ptol=1e-3,nmax=6) {
# Integrate over R using an extended trapezoid rule
# Set initial increments
  Rl <- qgamma(ptol, shape = kE, scale = R / kE )
  Ru <- qgamma(1-ptol, shape = kE, scale = R / kE )
  Rhalf <- qgamma(0.5, shape = kE, scale = R / kE )
# One whole increment first and the halfway point
  psum <- 0.5 * ( Ricker_dnbinom_d(Ntp1,Nt,Rl,alpha,kD)
                  + Ricker_dnbinom_d(Ntp1,Nt,Ru,alpha,kD) )
  psum <- psum + Ricker_dnbinom_d(Ntp1,Nt,Rhalf,alpha,kD)
# Progressively smaller increments (halving size each time)
  for (n in 2:nmax) {
    totincs <- 2^n                #The total number of increments (old + new)
    newpts <- 2^(n-1)             #The number of new points
    nwdel <- (1-2*ptol)/newpts    #New delta
    p_evals <- seq(ptol+0.5*nwdel,1-ptol-0.5*nwdel,length.out=newpts)
    R_evals <- qgamma(p_evals, shape = kE, scale = R / kE )
    pmids <- 0                    #Initialize to zero
    for (r in R_evals) {
      pmids <- pmids + Ricker_dnbinom_d(Ntp1,Nt,r,alpha,kD)
    }
    psum <- psum + pmids          #Add the new points to previous ones
    pinc <- (1-2*ptol)/totincs    #Current increment size on probability scale
    pint <- psum * pinc           #Current approximation for p
  }
  if (log==FALSE) {
    return( pint )
  } # else
  return( log(pint) )
}

#----Negative binomial-gamma pmf for Ricker model (version 2)
#
#This version has an automatic stopping rule. With the default eps and ptol settings
#this ends up always stopping at n=6. It is slower than version 1 due to looping.
#The integral uses Eulers-McLaurin trapezoid geometry on P[R=r] - i.e. 
#equal increments on the probability scale. This effectively gives us an
#adaptive step size - we take a smaller step size on the R scale where P[R=r]
#is high. There is an automatic stopping rule. This is based on the function
#DiffvarD_LAint.r, which in turn follows Numerical Recipes in C, trapzd.c and qtrap.c.
#It is a clear candidate for rewriting in C because the integral
#is calculated separately for each cell in the vector Ntp1.
#The ptol and eps values are relatively imprecise but sufficient for model fitting.
#
Ricker_dnbinomgamma2 <- function(Ntp1,Nt,R,alpha,kD,kE,log=FALSE,
                                 ptol=1e-3,eps=1e-2,nmax=20) {
# Do each N in the vector Ntp1 separately
  p <- NA*(1:length(Ntp1)) #To hold the probabilities
# Integrate over R using an extended trapezoid rule
# Set initial increments
  Rl <- qgamma(ptol, shape = kE, scale = R / kE )
  Ru <- qgamma(1-ptol, shape = kE, scale = R / kE )
  Rhalf <- qgamma(0.5, shape = kE, scale = R / kE )
  for (i in 1:length(Ntp1)) {
  # One whole increment first and the halfway point
    psum <- 0.5 * ( Ricker_dnbinom_d(Ntp1[i],Nt[i],Rl,alpha,kD)
                    + Ricker_dnbinom_d(Ntp1[i],Nt[i],Ru,alpha,kD) )
    psum <- psum + Ricker_dnbinom_d(Ntp1[i],Nt[i],Rhalf,alpha,kD)
  # Progressively smaller increments (halving size each time)
    oldp <- -1.0e30
    for (n in 2:nmax) {
      totincs <- 2^n                #The total number of increments (old + new)
      newpts <- 2^(n-1)             #The number of new points
      nwdel <- (1-2*ptol)/newpts    #New delta
      p_evals <- seq(ptol+0.5*nwdel,1-ptol-0.5*nwdel,length.out=newpts)
      R_evals <- qgamma(p_evals, shape = kE, scale = R / kE )
      pmids <- 0                    #Initialize to zero
      for (r in R_evals) {
        pmids <- pmids + Ricker_dnbinom_d(Ntp1[i],Nt[i],r,alpha,kD)
      }
      psum <- psum + pmids          #Add the new points to previous ones
      pinc <- (1-2*ptol)/totincs    #Current increment size on probability scale
      pint <- psum * pinc           #Current approximation for p
#print(cbind(n,pint)) #debug
    # Check for convergence
      if (n > 5){                   #Avoid spurious early convergence
      # For the most important values, this is the stopping rule
        if ( abs(pint-oldp) <= eps*abs(oldp) ) {
          p[i] <- pint
          break
        } #else
      # For tiny p, accept the result so far
        if ( pint < ptol ) {
          p[i] <- pint
          break
        } #else
      }
      oldp <- pint
      if (n==nmax) stop("Too many steps in integration")
    }
  }
  if (log==FALSE) {
    return( p )
  } # else
  return( log(p) )
}

#----Poisson-binomial pmf for Ricker model
#Translate to C? Can't avoid the loop, and this algorithm involves
#lots of memory allocation operations in R.
#
Ricker_dpoisbinom <- function(Ntp1,Nt,p,R,alpha,log=FALSE){
  sump <- rep(NA,length(Ntp1)) #Sum of the probabilities
  for (i in 1:length(Ntp1)){
    f <- 0:Nt[i] #f is number of females
    sump[i] <- sum( exp( dbinom(f,Nt[i],p,log=TRUE) + 
                       dpois(Ntp1[i],(1/p)*R*f*exp(-alpha*Nt[i]),log=TRUE) ) )
    #We use the log scale to ensure small probs multiply accurately
  }
  if (log==FALSE) {
    return( sump )
  } #else
  return( log(sump) )
}

#----Negative binomial-binomial (demographic) pmf for Ricker model
#Translate to C? Can't avoid the loop, and this algorithm involves
#Lots of memory allocation operations in R.
#
Ricker_dnbinombinom_d <- function(Ntp1,Nt,p,R,alpha,k,log=FALSE){
  sump <- rep(NA,length(Ntp1)) #Sum of the probabilities
  for (i in 1:length(Ntp1)){
    f <- 0:Nt[i] #f is number of females
    sump[i] <- sum( exp( dbinom(f,Nt[i],p,log=TRUE) + 
        dnbinom(Ntp1[i],size=k*f+(f==0),
                mu=(1/p)*R*f*exp(-alpha*Nt[i]),log=TRUE) ) )
    #We use the log scale to ensure small probs multiply accurately
    #Add 1 to size parameter of NB when females=0 to avoid NaN's.
  }
  if (log==FALSE) {
    return( sump )
  } #else
  return( log(sump) )
}

#----Negative binomial-binomial (environmental) pmf for Ricker model
#Translate to C? Can't avoid the loop, and this algorithm involves
#lots of memory allocation operations in R.
#
Ricker_dnbinombinom_e <- function(Ntp1,Nt,p,R,alpha,k,log=FALSE){
  sump <- rep(NA,length(Ntp1)) #Sum of the probabilities
  for (i in 1:length(Ntp1)){
    f <- 0:Nt[i] #f is number of females
    sump[i] <- sum( exp( dbinom(f,Nt[i],p,log=TRUE) + 
          dnbinom(Ntp1[i],size=k,mu=(1/p)*R*f*exp(-alpha*Nt[i]),log=TRUE) ) )
    #We use the log scale to ensure small probs multiply accurately
  }
  if (log==FALSE) {
    return( sump )
  } #else
  return( log(sump) )
}

#----Negative binomial-binomial-gamma pmf for Ricker model
#
#Requires a standard number of increments and has no automatic resolution
#setting. Because of this, the whole vector of Ntp1 can be done at once and
#it is faster. See notes for Ricker_dnbinomgamma if higher accuracy is needed.
#The integral uses Eulers-McLaurin trapezoid geometry on P[R=r] - i.e. 
#equal increments on the probability scale. This effectively gives us an
#adaptive step size - we take a smaller step size on the R scale where P[R=r]
#is high. This follows Numerical Recipes in C, trapzd.c and qtrap.c.
#The ptol and nmax values are quite course but sufficient for model fitting.
#
Ricker_dnbinombinomgamma <- function(Ntp1,Nt,p,R,alpha,kD,kE,log=FALSE,ptol=1e-3,nmax=6) {
# Integrate over R using an extended trapezoid rule
# Set initial increments
  Rl <- qgamma(ptol, shape = kE, scale = R / kE )
  Ru <- qgamma(1-ptol, shape = kE, scale = R / kE )
  Rhalf <- qgamma(0.5, shape = kE, scale = R / kE )
# One whole increment first and the halfway point
  psum <- 0.5 * ( Ricker_dnbinombinom_d(Ntp1,Nt,p,Rl,alpha,kD)
                  + Ricker_dnbinombinom_d(Ntp1,Nt,p,Ru,alpha,kD) )
  psum <- psum + Ricker_dnbinombinom_d(Ntp1,Nt,p,Rhalf,alpha,kD)
# Progressively smaller increments (halving size each time)
  for (n in 2:nmax) {
    totincs <- 2^n                #The total number of increments (old + new)
    newpts <- 2^(n-1)             #The number of new points
    nwdel <- (1-2*ptol)/newpts    #New delta
    p_evals <- seq(ptol+0.5*nwdel,1-ptol-0.5*nwdel,length.out=newpts)
    R_evals <- qgamma(p_evals, shape = kE, scale = R / kE )
    pmids <- 0                    #Initialize to zero
    for (r in R_evals) {
      pmids <- pmids + Ricker_dnbinombinom_d(Ntp1,Nt,p,r,alpha,kD)
    }
    psum <- psum + pmids          #Add the new points to previous ones
    pinc <- (1-2*ptol)/totincs    #Current increment size on probability scale
    pint <- psum * pinc           #Current approximation for p
  }
  if (log==FALSE) {
    return( pint )
  } # else
  return( log(pint) )
}


# Likelihoods for the Ricker model
# These return negative log likelihoods (nll)
# The .nll functions require the global variables Nt and Ntp1 

#source("Ricker_pmf.r") #For the pmfs (not used for first three .nlls)


#----Likelihood for the Poisson Ricker model
#
Ricker_pois.nll <- function(lnR,lnalpha){
  #Probability of the data
  lambda <- Nt*exp(lnR)*exp(-exp(lnalpha)*Nt)
  -sum(dpois(Ntp1,lambda,log=TRUE))
}

#----Likelihood for the Negative binomial (demographic) Ricker model
#
Ricker_nbinom_d.nll <- function(lnR,lnalpha,lnk){
  #Probability of the data
  mu <- Nt*exp(lnR)*exp(-exp(lnalpha)*Nt)
  -sum(dnbinom(Ntp1,size=exp(lnk)*Nt,mu=mu,log=TRUE))
}

#----Likelihood for the Negative binomial (environmental) Ricker model
#
Ricker_nbinom_e.nll <- function(lnR,lnalpha,lnk){
  #Probability of the data
  mu <- Nt*exp(lnR)*exp(-exp(lnalpha)*Nt)
  -sum(dnbinom(Ntp1,size=exp(lnk),mu=mu,log=TRUE))
}

#----Likelihood for the Negative binomial-gamma Ricker model
#
Ricker_nbinomgamma.nll <- function(lnR,lnalpha,lnkD,lnkE){
  #Probability of the data
  -sum(Ricker_dnbinomgamma(Ntp1,Nt,exp(lnR),exp(lnalpha),exp(lnkD),exp(lnkE),log=TRUE))
}

#----Likelihood for the Poisson-binomial Ricker model
#
Ricker_poisbinom.nll <- function(lnR,lnalpha){
  #Probability of the data
  -sum(Ricker_dpoisbinom(Ntp1,Nt,0.5,exp(lnR),exp(lnalpha),log=TRUE))
}

#----Likelihood for the Negative binomial-binomial (demographic) Ricker model
#
Ricker_nbinombinom_d.nll <- function(lnR,lnalpha,lnk){
  #Probability of the data
  -sum(Ricker_dnbinombinom_d(Ntp1,Nt,0.5,exp(lnR),exp(lnalpha),exp(lnk),log=TRUE))
}

#----Likelihood for the Negative binomial-binomial (environmental) Ricker model
#
Ricker_nbinombinom_e.nll <- function(lnR,lnalpha,lnk){
  #Probability of the data
  -sum(Ricker_dnbinombinom_e(Ntp1,Nt,0.5,exp(lnR),exp(lnalpha),exp(lnk),log=TRUE))
}

#----Likelihood for the Negative binomial-binomial-gamma Ricker model
#
Ricker_nbinombinomgamma.nll <- function(lnR,lnalpha,lnkD,lnkE){
  #Probability of the data
  -sum(Ricker_dnbinombinomgamma(Ntp1,Nt,0.5,exp(lnR),exp(lnalpha),
                                          exp(lnkD),exp(lnkE),log=TRUE))
}

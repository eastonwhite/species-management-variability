# This model allows births in the absence of males (i.e. when all are females)

# Ricker with stochasticity in birth (B), survival (S), and sex ratio.
# One might consider this to be the full model of demographic stochasticity.
# This has a Poisson-binomial distribution. Birth and survival are together
# in the composite parameter R, so this function cannot parse out birth
# and density independent survival. 
# p is the sex ratio (f/m).
#

RickerStSexBS <- function(Nt, R, alpha, p=0.5){
  females <- rbinom(length(Nt),Nt,p)
# Births and density independent survival (R = births*(1-mortality) and
# density dependent survival ( e^-alpha*Nt ); mortality is
# binomial, so compound distribution is Poisson with mean Re^-aNt).
  survivors <- rpois( length(Nt), females * (1/p) * R * exp( -alpha * Nt ) )
  return(survivors)
}


#Separates out density dependent survival
RickerStSexBS_longversion <- function(Nt, R, alpha, p=0.5){
  females <- rbinom(length(Nt),Nt,p)
# Births and density independent survival (R = births*(1-mortality); mortality is
# binomial, so compound distribution is Poisson with mean R).
  births <- rpois( length(Nt), females * (1/p) * R )
# Density dependent survival
  survivors <- rbinom( length(Nt), births, exp( -alpha * Nt ) )
  return(survivors)
}
# Ricker with stochasticity in birth (B), survival (S),
# and between individual (ie demographic) heterogeneity in B
# and spatio-temporal environmental heterogeneity in B.
# Demographic stochasticity (but not sex) plus demographic heterogeneity (Dh)
# and environmental heterogeneity (Eh).
# kD is the shape parameter of the gamma distribution for demographic
# heterogeneity.
# kE is the shape parameter of the gamma distribution for environmental
# heterogeneity.
# This results in a Negative binomial-gamma distribution.
#
RickerStBS_DEhB <- function(Nt, R, alpha, kD, kE) {
# Heterogeneity in birth rate between times or locations
  Rtx <- rgamma(length(Nt),shape=kE,scale=R/kE)
# Heterogeneity in individual birth rate plus density independent survival
# (R = births*(1-mortality); mortality is binomial, so compound distribution
# is negative binomial with mean R*Nt)
  births <- rnbinom( length(Nt), size = kD * Nt + (Nt==0), mu = Rtx * Nt )
  #Add 1 to size when Nt=0 to avoid NaN's.
# Density dependent survival
  survivors <- rbinom( length(Nt), births, exp( -1 * alpha * Nt ) )
  return(survivors)
}

# Full individual based model with gamma variation.
# This is for demonstration only. Use the above equivalent model.
#
RickerStBS_DhB_full <- function(Nt, R, alpha, kD, kE) {
# Heterogeneity in birth rate between times or locations
  Rtx <- rgamma(length(Nt),shape=kE,scale=R/kE)
# Heterogeneity in individual birth rate
  births <- NA*Nt
  for (tx in 1:length(Nt)) {
    births[tx] <- sum(  rpois( Nt[tx],rgamma(Nt[tx],shape=kD,scale=Rtx[tx]/kD) )  )
  }
  survivors <- rbinom( length(Nt), births, exp( -1 * alpha * Nt ) )
  return(survivors)
}
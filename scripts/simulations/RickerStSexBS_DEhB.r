# Ricker with stochasticity in birth (B), survival (S), sex ratio,
# and between individual (ie demographic) heterogeneity in B
# and spatio-temporal environmental heterogeneity in B.
# Demographic stochasticity plus demographic heterogeneity (Dh)
# and environmental heterogeneity (Eh).
# kD is the shape parameter of the gamma distribution for demographic
# heterogeneity.
# kE is the shape parameter of the gamma distribution for environmental
# heterogeneity.
# p is the sex ratio (f/m).
# This results in a Negative binomial-binomial-gamma distribution.
#
RickerStSexBS_DEhB <- function(Nt, R, alpha, kD, kE, p=0.5) {
# Heterogeneity in birth rate/DI mortality between times or locations
  Rtx <- rgamma(length(Nt),shape=kE,scale=R/kE)
  females <- rbinom(length(Nt),Nt,p)
# Heterogeneity in individual birth rate plus density independent survival
# (R = births*(1-mortality)*(1/p); mortality is binomial, so compound distribution
# is negative binomial with mean (1/p)*R*Nt)
  births <- rnbinom( length(Nt), size = kD * females + (females==0),
            mu = (1/p) * females * Rtx ) # Add 1 to size when females=0 to avoid NaN's.
# Density dependent survival
  survivors <- rbinom( length(Nt), births, exp( -1 * alpha * Nt ) )
  return(survivors)
}


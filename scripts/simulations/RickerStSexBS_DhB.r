# Ricker with stochasticity in birth (B), survival (S), sex ratio,
# and between individual (ie demographic) heterogeneity in B.
# In other words, demographic stochasticity plus demographic heterogeneity (Dh).
# k is the shape parameter of the gamma distribution.
# p is the sex ratio (f/m).
# This results in a Negative-binomial-binomial distribution.
#
RickerStSexBS_DhB <- function(Nt, R, alpha, k, p=0.5) {
  females <- rbinom(length(Nt),Nt,p)
# Heterogeneity in individual birth rate.
  births <- rnbinom( length(Nt), size = k * females + (females==0),
            mu = females * (1/p) * R ) #Add 1 to size when females=0 to avoid NaN's.
# Survival
  survivors <- rbinom( length(Nt), births, exp( -1 * alpha * Nt ) )
  return(survivors)
}
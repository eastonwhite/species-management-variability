# Ricker with stochasticity in birth (B), survival (S), sex ratio,
# and between patch or time (ie environmental) heterogeneity in B.
# In other words, demographic stochasticity plus environmental heterogeneity (Eh).
# k is the shape parameter of the gamma distribution.
# p is the sex ratio (f/m).
# This results in a Negative-binomial-binomial distribution.
#
RickerStSexBS_EhB <- function(Nt, R, alpha, k, p=0.5) {
  females <- rbinom(length(Nt),Nt,p)
# Heterogeneity in individual birth rate.
  births <- rnbinom( length(Nt), size = k, mu = females * (1/p) * R )
# Survival
  survivors <- rbinom( length(Nt), births, exp( -1 * alpha * Nt ) )
  return(survivors)
}
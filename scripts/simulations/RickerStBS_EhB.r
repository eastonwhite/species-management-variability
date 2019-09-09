# Ricker with stochasticity in birth (B), survival (S),
# and between patch or time (ie environmental) heterogeneity in B.
# Demographic stochasticity (but not sex) plus environmental heterogeneity (Eh).
# k is the shape parameter of the gamma distribution.
# This results in a Negative binomial distribution with mu=Ricker, k'=k,
# and theoretical variance: Ricker + Ricker^2/(k).
#
RickerStBS_EhB <- function(Nt, R, alpha, k) {
# Heterogeneity in individual birth rate plus density independent survival
# (R = births*(1-mortality); mortality is binomial, so compound distribution
# is negative binomial with mean R*Nt)
  births <- rnbinom( length(Nt), size = k, mu = R * Nt )
# Density dependent survival
  survivors <- rbinom( length(Nt), births, exp( -1 * alpha * Nt ) )
  return(survivors)
}

# Full individual based model with gamma variation.
# This is for demonstration only. Use the above equivalent model.
#
RickerStBS_EhB_full <- function(Nt, R, alpha, k) {
# Heterogeneity in birth rate between times or locations
  Rx <- rgamma(length(Nt),shape=k,scale=R/k)
  births <- NA*Nt   #vector to hold births
  for (x in 1:length(Nt)) {
    births[x] <- sum(  rpois( Nt[x],Rx[x] )  )
  }
  survivors <- rbinom( length(Nt), births, exp( -1 * alpha * Nt ) )
  return(survivors)
}
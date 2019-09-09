# Ricker with stochasticity in birth (B), survival (S),
# and between individual (ie demographic) heterogeneity in B.
# Demographic stochasticity (but not sex) plus demographic heterogeneity (Dh).
# k is the shape parameter of the gamma distribution.
# This results in a Negative binomial distribution with k'=k*Nt,
# and theoretical variance: Ricker + Ricker^2/(k*Nt).
#
RickerStBS_DhB <- function(Nt, R, alpha, k) {
# Heterogeneity in individual birth rate plus density independent survival
# (R = births*(1-mortality); mortality is binomial, so compound distribution
# is negative binomial with mean R*Nt).
  births <- rnbinom( length(Nt), size = k * Nt + (Nt==0), mu = R * Nt )
  #Add 1 to size when Nt=0 to avoid NaN's.
# Density dependent survival
  survivors <- rbinom( length(Nt), births, exp( -1 * alpha * Nt ) )
  return(survivors)
}

# Full individual based model with gamma variation.
# This is for demonstration only. Use the above equivalent model.
#
RickerStBS_DhB_full <- function(Nt, R, alpha, k) {
# Heterogeneity in individual birth rate (x indicates patch or time).
  births <- NA*Nt
  for (x in 1:length(Nt)) {
    births[x] <- sum(  rpois( Nt[x],rgamma(Nt[x],shape=k,scale=R/k) )  )
  }
  survivors <- rbinom( length(Nt), births, exp( -1 * alpha * Nt ) )
  return(survivors)
}
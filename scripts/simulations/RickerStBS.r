# Ricker with stochasticity in birth (B) and survival (S).
# This is the basic Ricker model with demographic stochasticity.
# This has a Poisson distribution with lambda = deterministic Ricker.
#
RickerStBS <- function(Nt, R, alpha){
# Births and density independent survival (R = births*(1-mortality); mortality is
# binomial, so compound distribution is Poisson with mean R).
  births <- rpois( length(Nt), Nt * R )
# Density dependent survival
  survivors <- rbinom( length(Nt), births, exp( -1 * alpha * Nt ) )
  return(survivors)
}
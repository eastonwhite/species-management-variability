#dirichletmultinom.r
#The pmf for the Dirichlet-multinomial compound distribution, otherwise known
#as the multivariate Polya distribution. In Bayesian modeling, the Dirichlet
#is the conjugate prior for the multinomial. Here the full compound distribution
#is implemented.
#
#             sum(x)!        gamma(sum(a))       prod(gamma(x+a))
# p(x|a) = ------------- ---------------------- ------------------
#           prod(x[i]!)   gamma(sum(x)+sum(a))    prod(gamma(a))
#
#See Minka(2003) web only publication (not peer reviewed) but note that
#his pmf is missing the multinomial coefficient (the first term above).              
#For the Dirichlet and multinomial in the Bayesian context see Clark (2007).
#
#ddirichlet and rdirichlet are available in the gtools package.
#
#This function currently only accepts a vector x and therefore only returns a
#a scalar p. If we wanted to look at multiple multinomial trials, we need to
#modify to allow a matrix x, with each trial in a row.
#See, e.g, the implementation of ddirichlet in gtools.
#
dpolya <- function(x, alpha, log=FALSE) {
  suma <- sum(alpha)
  sumx <- sum(x)
  lnp <- lgamma(suma) - lgamma(sumx+suma) + sum(lgamma(x+alpha)) - sum(lgamma(alpha))
                                              #are the latter 2 more accurate summed first, as here?
  lnp <- lnp + lfactorial(sumx) - sum(lfactorial(x)) #multinomial coefficient
  if (log == FALSE){
    return(exp(lnp))
  }else{
    return(lnp)
  }
}

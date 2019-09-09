#This is an extract of the R code used in the Nature 2008 paper.
#Brett Melbourne
#10 Jul 14

mydata <- read.csv("ricker_data.csv")
#Round abundance to nearest integer. We'll use Nt for adults in publications.
#mydata$Nt <- round(mydata$At)
#mydata$Ntp1 <- round(mydata$Atp1)
head(mydata)
attach(mydata)
plot(Nt,Ntp1)
max(Ntp1)
detach(mydata)

xm <- 1000 #maximum density (for x axis limit)
ym <- 370 #maximum density (for y axis limit)
xmt <- 100 #Intervals for x minor ticks
ymt <- 50 #Intervals for y minor ticks

#Try leaving out influential/outliers (case 19, 24, 58)
#mydata <- mydata[-58,]
#mydata <- mydata[-c(19,24,58),]
#Or try leaving out all the points at the max of prod func (viz reviewer 1)
#mydata <- mydata[!(mydata$Nt %in% c(358,363,360)),]


#Fitting a Ricker by maximum likelihood

attach(mydata)
source("Ricker.r")
source("Ricker_nll.r")
library(stats4) #mle

#Likelihood for Poisson
llfit <- mle( Ricker_pois.nll,start=list(lnR=log(1.5),lnalpha=log(0.001)) )
summary(llfit)
phat <- exp(coef(llfit))
names(phat) <- list("R","alpha")
phat
logLik(llfit)
vcov(llfit)
plot(profile(llfit), absVal=FALSE)
confint(llfit)
Poiss_AIC <- AIC(llfit)
Poiss_AIC

#Likelihood for Negative binomial (demographic)
llfit <- mle( Ricker_nbinom_d.nll,start=list(lnR=log(3),lnalpha=log(0.001),lnk=log(1)) )
summary(llfit)
phat <- exp(coef(llfit))
names(phat) <- list("R","alpha","k")
phat
logLik(llfit)
vcov(llfit)
plot(profile(llfit), absVal=FALSE)
confint(llfit)
NB_AIC <- AIC(llfit)
NB_AIC

#Likelihood for Negative binomial (environmental)
llfit <- mle( Ricker_nbinom_e.nll,start=list(lnR=log(3),lnalpha=log(0.003),lnk=log(1)) )
summary(llfit)
phat <- exp(coef(llfit))
names(phat) <- list("R","alpha","k")
phat
logLik(llfit)
vcov(llfit)
plot(profile(llfit), absVal=FALSE)
confint(llfit)
NB_AIC <- AIC(llfit)
NB_AIC

#Likelihood for Negative binomial-gamma
llfit <- mle( Ricker_nbinomgamma.nll,start=list(lnR=log(3),lnalpha=log(0.003),
             lnkD=log(1),lnkE=log(1)) )
summary(llfit)
phat <- exp(coef(llfit))
names(phat) <- list("R","alpha","kD","kE")
phat
logLik(llfit)
vcov(llfit)
plot(profile(llfit), absVal=FALSE)
confint(llfit)
NBg_AIC <- AIC(llfit)
NBg_AIC

#Likelihood for Poisson-binomial
llfit <- mle( Ricker_poisbinom.nll,start=list(lnR=log(3),lnalpha=log(0.003)) )
summary(llfit)
phat <- exp(coef(llfit))
names(phat) <- list("R","alpha")
phat
logLik(llfit)
vcov(llfit)
plot(profile(llfit), absVal=FALSE)
confint(llfit)
PB_AIC <- AIC(llfit)
PB_AIC


#Likelihood for Negative binomial-binomial (demographic)
llfit <- mle( Ricker_nbinombinom_d.nll,start=list(lnR=log(2.64),lnalpha=log(0.001),lnk=log(1)) )
summary(llfit)
phat <- exp(coef(llfit))
names(phat) <- list("R","alpha","k")
phat
logLik(llfit)
vcov(llfit)
plot(profile(llfit), absVal=FALSE)
confint(llfit)
NBBd_AIC <- AIC(llfit)
NBBd_AIC

#Likelihood for Negative binomial-binomial (environmental)
llfit <- mle( Ricker_nbinombinom_e.nll,start=list(lnR=log(2.64),lnalpha=log(0.001),lnk=log(1)) )
summary(llfit)
phat <- exp(coef(llfit))
names(phat) <- list("R","alpha","k")
phat
logLik(llfit)
vcov(llfit)
plot(profile(llfit), absVal=FALSE)
confint(llfit)
NBBe_AIC <- AIC(llfit)
NBBe_AIC

#Likelihood for Negative binomial-binomial-gamma
llfit <- mle( Ricker_nbinombinomgamma.nll,start=list(lnR=log(1),lnalpha=log(0.001),
             lnkD=log(1),lnkE=log(1)) )
summary(llfit)
phat <- exp(coef(llfit))
names(phat) <- list("R","alpha","kD","kE")
phat
logLik(llfit)
vcov(llfit)
plot(profile(llfit), absVal=FALSE)
confint(llfit)
NBBg_AIC <- AIC(llfit)
NBBg_AIC



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


endings=RickerStSexBS_DEhB(initial_pop, 3.57, 0.00199, 2.55, 17.69, p=0.5) 
plot(initial_pop,endings)
points(Nt,Ntp1,col='red')

#This program is to analyse the dispersal data for the Science paper.
#24 trials are used. These are stock culture samples for 3 generations
#from the LT 2 experiment. They are from 3 later generations.
#We use only the "Run 1" assays, which correspond to the experimental
#conditions.
#Models are fitted by maximum likelihood. The likelihoods are multinomial,
#or dirichlet-multinomial (=Polya).

source("Diff.r") #The "Poisson" diffusion model
source("DiffvarD_LAint.r") #The "Poisson-gamma" diffusion model
source("DiffSt.r") #Poisson simulator
source("dirichletmultinom.r") #for dpolya
library(stats4)
library(emdbook) #for betabinomial
library(gtools)  #for rdirichlet


mydata <- read.csv("EVEstock_assays.csv")
mydata[1,]

attach(mydata)
plot(jitter(Patch,0.2),Num_beetles)
boxplot(Num_beetles~Patch)
detach(mydata)

gen <- unique(mydata$Generation) #the 3 generations in which stock runs were done
col <- 1
plot(1,1,type="n",xlim=c(1,7),ylim=c(0,70),xlab="Patch",ylab="Abundance")
for (g in gen) {
  gendata <- subset(mydata,Generation==g)
  for (i in unique(gendata$Landscape)) {
    plotdata <- subset(gendata,Landscape==i)
    attach(plotdata)
    p <- jitter(Patch,0.2)
    points(p,Num_beetles,col=col)
    lines(p,Num_beetles,col=col)
    detach(plotdata)
  }
  col <- col + 1
}
legend(4,70, gen, lty=rep(1,3), pch=rep(1,3),col=1:3)

#Reshape the data into a matrix: each row is 1 landscape, cols are patches
mydatawide <- reshape(mydata[,-c(2:5,7,10)],
                      idvar=c("Generation","Landscape"),
                      timevar="Patch",direction="wide")
datamat <- as.matrix(mydatawide[,-(1:2)])
rownames(datamat) <- colnames(datamat) <- NULL #cleanup


#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# Likelihood for Poisson diffusion kernel with fixed initial
# abundance, which gives a multinomial distribution.
#
diffmultinompois.nll <- function(D){
  kbinompois <- Diff(c(1,rep(0,6)),1,D) 
  lnlik <- apply(datamat,1,dmultinom,prob=kbinompois,log=TRUE)
  return(-sum(lnlik))
}

llfit <- mle( diffmultinompois.nll,start=list(D=0.5) )
summary(llfit)
phat <- coef(llfit)
logLik(llfit)
ci <- confint(llfit)
ci
AIC(llfit)
matplot(1:7,t(datamat),type="n",xlab="Patch",ylab="Abundance")
#Overlay simulation
Ninit <- rowSums(datamat)
for (i in 1:20){
  n <- sample(Ninit,1,replace=TRUE) #Bootstrap n
  points(jitter((1:7),0.2),DiffSt(c(n,rep(0,6)),1,phat,theta_D=0),col="red")
}
matpoints(1:7,t(datamat),pch=16,col=1)
lines(1:7,Diff(c(70,rep(0,6)),1,phat,theta_D=0))

#Likelihood profile
grid <- seq(0.7,0.99,length.out=200)
profile <- grid*NA
for (i in 1:length(grid)) {
  profile[i] <- diffmultinompois.nll(grid[i])
}
plot(grid,profile,xlab="D",ylab="-Log likelihood")
grid()
cfint <- min(profile) + 1.92
abline(h=cfint,col="red") #Asymptotic confidence interval


#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# Likelihood for Poisson diffusion kernel with fixed initial
# abundance, allowing for overdispersion using a Polya
# distribution. The following parameterization is used for the
# Polya:
#   alpha = s * p
# where s is the "dispersion parameter" and p are the mean
# probabilities in each class (here p = kbinompois).
#
diffpolyapois.nll <- function(D,s){
#  s <- exp(lns)  #Sometimes helpful to constrain positive
  kbinompois <- Diff(c(1,rep(0,6)),1,D) 
  lnlik <- apply(datamat,1,dpolya,alpha=s*kbinompois,log=TRUE)
  return(-sum(lnlik))
}

llfit <- mle( diffpolyapois.nll,start=list(D=0.5,s=10) )
summary(llfit)
phat <- coef(llfit)
logLik(llfit)
ci <- confint(llfit) #The profile fails on s. Need to do by hand.
ci
AIC(llfit)
set.seed(4) #4 for science paper
matplot(1:7,t(datamat),type="n",ylim=c(0,55),xlab="Patch",ylab="Abundance")
#Overlay simulation
Ninit <- rowSums(datamat)
kbinompois <- Diff(c(1,rep(0,6)),1,phat[1]) #Kernel can be outside loop
for (i in 1:24){
  n <- sample(Ninit,1,replace=TRUE) #Bootstrap n
  tmp <- rdirichlet(1,phat[2]*kbinompois)
  sim <- rmultinom(1,n,tmp)
  lines(jitter((1:7),0.2),sim,col="red")
}
matlines(1:7,t(datamat),pch=1,col="blue",lty=1)
lines(1:7,Diff(c(70,rep(0,6)),1,phat[1]),col="black",lwd=2)
legend(5,52,c("Data","Kernel","Simulation"),lty=1,
       col=c("blue","black","red"))

#The 2D profile
library(rgl)
lnDseq <- seq(log(0.7), log(0.98), length=100)
lnsseq <- seq(log(10), log(35), length=100)
par_grid <- expand.grid( lnD = lnDseq, lns = lnsseq )
nll_grid <- rep(NA,nrow(par_grid))
for (i in 1:nrow(par_grid)){
  nll_grid[i] <- diffpolyapois.nll(exp(par_grid$lnD[i]),exp(par_grid$lns[i]))
}
res <- cbind(par_grid,nll_grid)
zlen <- max(-res$nll_grid)- min(-res$nll_grid) + 1
colorlut <- terrain.colors(zlen)
col <- colorlut[-res$nll_grid - min(-res$nll_grid) + 1]
plot3d(res$lnD, res$lns, res$nll_grid,size=1,col=col)


#I checked the profile CI for D and it matches the confint call above.

#Profile confidence interval for s
#We don't need a new nll function.

#Profile s
nll.s <- rep(NA,50)
s_range <- seq(12,32,length.out=length(nll.s))
par <- 0.83  #starting value for D
i <- 1
for (s in s_range){
  nll.s[i] <- optim( par, diffpolyapois.nll, s=s)$value
  print(i) #Monitor progress
  i <- i + 1
}

#Plot the profile
plot(s_range,nll.s,xlab="s",ylab="Neg Log Lik",
     type="l",col="blue",main="Confidence interval for s")
cbind(s_range,nll.s)

#Calculate the value of the nll that defines the interval
conf_lev <- 0.95  #0.95 will give a 95% confidence interval
min_nll <- - logLik(llfit)
conf_lim <- min_nll + qchisq(p=conf_lev,df=1)/2

abline(h=min_nll,col="red")
abline(h=conf_lim,col="red")

#Calculate the lower limit
nll_lo <- nll.s[1:which.min(nll.s)]
s_lo <- s_range[1:which.min(nll.s)]
conf_lo <- approx(nll_lo,s_lo,xout=conf_lim)$y
nll_hi <- nll.s[which.min(nll.s):length(nll.s)]
s_hi <- s_range[which.min(nll.s):length(nll.s)]
conf_hi <- approx(nll_hi,s_hi,xout=conf_lim)$y

#Add CI to graph
abline(v=conf_lo,col="red",lty=2)
abline(v=conf_hi,col="red",lty=2)

text(conf_lo,min_nll+2,round(conf_lo,4),pos=4)
text(conf_hi,min_nll+2,round(conf_hi,4),pos=2)

#The confidence interval
cbind(conf_lo,conf_hi)


#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# Likelihood for Poisson-gamma diffusion kernel with fixed
# initial abundance, which gives a multinomial distribution.
#
diffmultinompoisgamma.nll <- function(lnD,lnk){
  kpoisgamma <- DiffvarD_LAint(c(1,rep(0,6)),1,exp(lnD),exp(lnk))
  lnlik <- apply(datamat,1,dmultinom,prob=kpoisgamma,log=TRUE)
  return(-sum(lnlik))
}

llfit <-mle( diffmultinompoisgamma.nll,start=list(lnD=log(0.8),lnk=log(10000000)))
             #,method="Nelder-Mead" ) #if BFGS barfs
summary(llfit)
phat <- exp(coef(llfit))
phat
logLik(llfit)
cil <- confint(llfit) #Inf is the upper limit, so barfs. Need to do by hand.
cil
AIC(llfit)

#The 2D profile
library(rgl)
lnDseq <- seq(log(0.75), log(0.95), length=100)
lnkseq <- seq(log(55), log(100000), length=100)
par_grid <- expand.grid( lnD = lnDseq, lnk = lnkseq )
nll_grid <- rep(NA,nrow(par_grid))
for (i in 1:nrow(par_grid)){
  nll_grid[i] <- diffmultinompoisgamma.nll(par_grid$lnD[i],par_grid$lnk[i])
}
res <- cbind(par_grid,nll_grid)
zlen <- max(-res$nll_grid)- min(-res$nll_grid) + 1
colorlut <- terrain.colors(zlen)
col <- colorlut[-res$nll_grid - min(-res$nll_grid) + 1]
plot3d(res$lnD, res$lnk, -res$nll_grid,size=1,col=col)
#The surface is extremely flat in the lnk dimension - i.e. Inf is a likely value

#Profile confidence interval for k
diffmultinompoisgamma.nll.lnk <- function(lnD,lnk){
  kpoisgamma <- DiffvarD_LAint(c(1,rep(0,6)),1,exp(lnD),exp(lnk))
  lnlik <- apply(datamat,1,dmultinom,prob=kpoisgamma,log=TRUE)
  return(-sum(lnlik))
}

#Profile lnk
nll.lnk <- rep(NA,50)
lnk_range <- seq(3.3,10,length.out=length(nll.lnk))
par <- log(0.84)  #starting value for lnD
i <- 1
for (lnk in lnk_range){
  nll.lnk[i] <- optim( par, diffmultinompoisgamma.nll.lnk, lnk=lnk)$value
  print(i) #Monitor progress
  i <- i + 1
}

#Plot the profile
plot(lnk_range,nll.lnk,xlab="lnk",ylab="Neg Log Lik",
     type="l",col="blue",main="Confidence interval for lnk")
cbind(lnk_range,nll.lnk)

#It looks like the nll flatlines at large values. We can't directly check
#to see if infinity is the upper bound for lnk but we don't get any change
#in the likelihood for k = 1e100.
nll_Inf <- optim( par, diffmultinompoisgamma.nll.lnk, lnk=log(1e100))$value
nll_Inf

#Calculate the value of the nll that defines the interval
conf_lev <- 0.95  #0.95 will give a 95% confidence interval
min_nll <- - logLik(llfit)
conf_lim <- min_nll + qchisq(p=conf_lev,df=1)/2

abline(h=min_nll,col="red")
abline(h=conf_lim,col="red")

#Calculate the lower limit
nll_lo <- nll.lnk[1:which.min(nll.lnk)]
lnk_lo <- lnk_range[1:which.min(nll.lnk)]
conf_lo <- approx(nll_lo,lnk_lo,xout=conf_lim)$y
exp(conf_lo) #30.19 to greater accuracy with finer grid

#Profile confidence interval for D
diffmultinompoisgamma.nll.D <- function(D,lnk){
  kpoisgamma <- DiffvarD_LAint(c(1,rep(0,6)),1,D,exp(lnk))
  lnlik <- apply(datamat,1,dmultinom,prob=kpoisgamma,log=TRUE)
  return(-sum(lnlik))
}

#Profile D
nll.D <- rep(NA,50)
D_range <- seq(0.75,0.93,length.out=length(nll.D))
par <- 10  #starting value for lnk
i <- 1
for (D in D_range){
  nll.D[i] <- optim( par, diffmultinompoisgamma.nll.D, D=D)$value
  print(i) #Monitor progress
  i <- i + 1
}

#Plot the profile
plot(D_range,nll.D,xlab="D",ylab="Neg Log Lik",
     type="l",col="blue",main="Confidence interval for D")
cbind(D_range,nll.D)

#Calculate the value of the nll that defines the interval
conf_lev <- 0.95  #0.95 will give a 95% confidence interval
min_nll <- - logLik(llfit)
conf_lim <- min_nll + qchisq(p=conf_lev,df=1)/2

abline(h=min_nll,col="red")
abline(h=conf_lim,col="red")

#Calculate the limits
nll_lo <- nll.D[1:which.min(nll.D)]
D_lo <- D_range[1:which.min(nll.D)]
conf_lo <- approx(nll_lo,D_lo,xout=conf_lim)$y
nll_hi <- nll.D[which.min(nll.D):length(nll.D)]
D_hi <- D_range[which.min(nll.D):length(nll.D)]
conf_hi <- approx(nll_hi,D_hi,xout=conf_lim)$y

#Add CI to graph
abline(v=conf_lo,col="red",lty=2)
abline(v=conf_hi,col="red",lty=2)

text(conf_lo,min_nll+2,round(conf_lo,4),pos=4)
text(conf_hi,min_nll+2,round(conf_hi,4),pos=2)

#The confidence interval
cbind(conf_lo,conf_hi)


#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# Likelihood for Poisson-gamma diffusion kernel with fixed
# initial abundance, allowing for overdispersion using a Polya
# distribution. The following parameterization is used for the
# Polya:
#   alpha = s * p
# where s is the "dispersion parameter" and p are the mean
# probabilities in each class (here p = kpoisgamma).
#
diffpolyapoisgamma.nll <- function(lnD,lnk,lns){
  kpoisgamma <- DiffvarD_LAint(c(1,rep(0,6)),1,exp(lnD),exp(lnk))
  lnlik <- apply(datamat,1,dpolya,alpha=exp(lns)*kpoisgamma,log=TRUE)
  return(-sum(lnlik))
}

llfit <-mle( diffpolyapoisgamma.nll,start=list(lnD=log(0.8),lnk=log(1000000),
             lns=log(20)))
             #,method="Nelder-Mead" ) #if BFGS barfs
summary(llfit)
phat <- exp(coef(llfit))
names(phat) <- c("D","k","s")
phat
logLik(llfit)
ci <- confint(llfit) #Inf is the upper limit, so barfs. Need to do by hand.
ci
AIC(llfit)


#Profile confidence interval for k
#Modified nll function for call to optim
diffpolyapoisgamma.nll.lnk <- function(par,lnk){
  lnD <- par[1]; lns <- par[2]
  kpoisgamma <- DiffvarD_LAint(c(1,rep(0,6)),1,exp(lnD),exp(lnk))
  lnlik <- apply(datamat,1,dpolya,alpha=exp(lns)*kpoisgamma,log=TRUE)
  return(-sum(lnlik))
}


#Profile lnk
nll.lnk <- rep(NA,50)
lnk_range <- seq(2,10,length.out=length(nll.lnk))
par <- c(log(0.83),log(19.4))  #starting value for lnD and lns
i <- 1
for (lnk in lnk_range){
  nll.lnk[i] <- optim( par, diffpolyapoisgamma.nll.lnk, lnk=lnk)$value
  print(i) #Monitor progress
  i <- i + 1
}

#Plot the profile
plot(lnk_range,nll.lnk,xlab="lnk",ylab="Neg Log Lik",
     type="l",col="blue",main="Confidence interval for lnk")
cbind(lnk_range,nll.lnk)

#It looks like the nll flatlines at large values. We can't directly check
#to see if infinity is the upper bound for lnk but we don't get any change
#in the likelihood for k = 1e100.
nll_Inf <- optim( par, diffpolyapoisgamma.nll.lnk, lnk=log(1e100))$value
nll_Inf

#Calculate the value of the nll that defines the interval
conf_lev <- 0.95  #0.95 will give a 95% confidence interval
min_nll <- - logLik(llfit)
conf_lim <- min_nll + qchisq(p=conf_lev,df=1)/2

abline(h=min_nll,col="red")
abline(h=conf_lim,col="red")

#Calculate the lower limit
nll_lo <- nll.lnk[1:which.min(nll.lnk)]
lnk_lo <- lnk_range[1:which.min(nll.lnk)]
conf_lo <- approx(nll_lo,lnk_lo,xout=conf_lim)$y

#Add CI to graph
abline(v=conf_lo,col="red",lty=2)
text(conf_lo,min_nll+2,round(conf_lo,4),pos=4)

#The lower limit
exp(conf_lo) #30.19 to greater accuracy with finer grid


#Profile confidence interval for D
#Modified nll function for call to optim
diffpolyapoisgamma.nll.lnD <- function(par,lnD){
  lnk <- par[1]; lns <- par[2]
  kpoisgamma <- DiffvarD_LAint(c(1,rep(0,6)),1,exp(lnD),exp(lnk))
  lnlik <- apply(datamat,1,dpolya,alpha=exp(lns)*kpoisgamma,log=TRUE)
  return(-sum(lnlik))
}

#Profile D
nll.lnD <- rep(NA,50)
lnD_range <- seq(log(0.7),log(1),length.out=length(nll.lnD))
par <- c(log(1e6),log(19.36))  #starting value for lnk & lns
#par <- c(log(10),log(19.36))  #same answer with smaller k
i <- 1
for (lnD in lnD_range){
  nll.lnD[i] <- optim( par, diffpolyapoisgamma.nll.lnD, lnD=lnD)$value
  print(i) #Monitor progress
  i <- i + 1
}

#Plot the profile
plot(lnD_range,nll.lnD,xlab="lnD",ylab="Neg Log Lik",
     type="l",col="blue",main="Confidence interval for lnD")
cbind(lnD_range,nll.lnD)

#Calculate the value of the nll that defines the interval
conf_lev <- 0.95  #0.95 will give a 95% confidence interval
min_nll <- - logLik(llfit)
conf_lim <- min_nll + qchisq(p=conf_lev,df=1)/2

abline(h=min_nll,col="red")
abline(h=conf_lim,col="red")

#Calculate the limits
nll_lo <- nll.lnD[1:which.min(nll.lnD)]
lnD_lo <- lnD_range[1:which.min(nll.lnD)]
conf_lo <- approx(nll_lo,lnD_lo,xout=conf_lim)$y
nll_hi <- nll.lnD[which.min(nll.lnD):length(nll.lnD)]
lnD_hi <- lnD_range[which.min(nll.lnD):length(nll.lnD)]
conf_hi <- approx(nll_hi,lnD_hi,xout=conf_lim)$y

#Add CI to graph
abline(v=conf_lo,col="red",lty=2)
abline(v=conf_hi,col="red",lty=2)

text(conf_lo,min_nll+2,round(conf_lo,4),pos=4)
text(conf_hi,min_nll+2,round(conf_hi,4),pos=2)

#The confidence interval
cbind(conf_lo,conf_hi)
exp(cbind(conf_lo,conf_hi))


#Profile confidence interval for s
#Modified nll function for call to optim
diffpolyapoisgamma.nll.lns <- function(par,lns){
  lnD <- par[1]; lnk <- par[2]
  kpoisgamma <- DiffvarD_LAint(c(1,rep(0,6)),1,exp(lnD),exp(lnk))
  lnlik <- apply(datamat,1,dpolya,alpha=exp(lns)*kpoisgamma,log=TRUE)
  return(-sum(lnlik))
}

#Profile s
nll.lns <- rep(NA,50)
lns_range <- seq(log(12),log(31),length.out=length(nll.lns))
par <- c(log(0.83),log(1e6))  #starting value for lnD & lnk
#par <- c(log(0.83),log(10))  #almost same answer with smaller k
i <- 1
for (lns in lns_range){
  nll.lns[i] <- optim( par, diffpolyapoisgamma.nll.lns, lns=lns)$value
  print(i) #Monitor progress
  i <- i + 1
}

#Plot the profile
plot(lns_range,nll.lns,xlab="lns",ylab="Neg Log Lik",
     type="l",col="blue",main="Confidence interval for lns")
cbind(lns_range,nll.lns)

#Calculate the value of the nll that defines the interval
conf_lev <- 0.95  #0.95 will give a 95% confidence interval
min_nll <- - logLik(llfit)
conf_lim <- min_nll + qchisq(p=conf_lev,df=1)/2

abline(h=min_nll,col="red")
abline(h=conf_lim,col="red")

#Calculate the limits
nll_lo <- nll.lns[1:which.min(nll.lns)]
lns_lo <- lns_range[1:which.min(nll.lns)]
conf_lo <- approx(nll_lo,lns_lo,xout=conf_lim)$y
nll_hi <- nll.lns[which.min(nll.lns):length(nll.lns)]
lns_hi <- lns_range[which.min(nll.lns):length(nll.lns)]
conf_hi <- approx(nll_hi,lns_hi,xout=conf_lim)$y

#Add CI to graph
abline(v=conf_lo,col="red",lty=2)
abline(v=conf_hi,col="red",lty=2)

text(conf_lo,min_nll+2,round(conf_lo,4),pos=4)
text(conf_hi,min_nll+2,round(conf_hi,4),pos=2)

#The confidence interval
cbind(conf_lo,conf_hi)
exp(cbind(conf_lo,conf_hi))









#I am pretty sure that diffpolyapois.nll is the code that was the best fit in the 2009 paper
#It is clear there were much lower levels of dispersal in our work. Is this just a species thing?

# Prepare data for analysis
dispersal48 = subset(dispersal,dispersal$dTime==48)
dispersal24 = subset(dispersal,dispersal$dTime==24)
dispersal12 = subset(dispersal,dispersal$dTime==12)
dispersal6 = subset(dispersal,dispersal$dTime==6)
dispersal1 = subset(dispersal,dispersal$dTime==1)

#dispersal6 = subset(dispersal6,dispersal6$landscape != 7)

#plot(jitter(Patch,0.2),Num_beetles/max(Num_beetles),ylab='Number of beetles',xlab='Patch number',cex.lab=1.2)
#plot(jitter(dispersal48$patch,0.2),dispersal48$final_count/max(dispersal48$final_count),col='red',pch=3)


mydatawide <- reshape(dispersal6[,c('landscape','patch','final_count')],
                      idvar=c("landscape"),
                      timevar="patch",direction="wide")
datamat <- as.matrix(mydatawide[,-(1)])
rownames(datamat) <- colnames(datamat) <- NULL
#datamat = cbind(datamat,rep(0,12))
datamat= round(datamat)
#datamat=rbind(datamat,datamat)

diffmultinompois.nll <- function(D){
  kbinompois <- Diff(c(1,rep(0,5)),1/8,D) 
  lnlik <- apply(datamat,1,dmultinom,prob=kbinompois,log=TRUE)
  return(-sum(lnlik))
}



llfit <- mle( diffmultinompois.nll,start=list(D=0.1) )
summary(llfit)
phat <- coef(llfit)
print(phat)
#logLik(llfit)
#ci <- confint(llfit)
#ci
#AIC(llfit)
matplot(1:6,t(datamat),type="n",xlab="Patch",ylab="Abundance",ylim=c(0,30))
#Overlay simulation
Ninit <- rowSums(datamat)
for (i in 1:100){
  n <- sample(Ninit,1,replace=TRUE) #Bootstrap n
  points(jitter((1:6),0.2),DiffSt(c(n,rep(0,5)),1/8,phat,theta_D=0),col="red")
}
matpoints(1:6,t(datamat),pch=16,col=1)
#lines(1:6,Diff(c(200,rep(0,5)),1,phat))

#Likelihood profile
grid <- seq(0.1,0.99,length.out=200)
profile <- grid*NA
for (i in 1:length(grid)) {
  profile[i] <- diffmultinompois.nll(grid[i])
}
plot(grid,profile,xlab="D",ylab="-Log likelihood")
grid()
cfint <- min(profile) + 1.92
abline(h=cfint,col="red") #Asymptotic confidence interval


####################################

# Poisson with overdispersion
#Likelihood for Poisson diffusion kernel with fixed initial
# abundance, allowing for overdispersion using a Polya
# distribution. The following parameterization is used for the
# Polya:
#   alpha = s * p
# where s is the "dispersion parameter" and p are the mean
# probabilities in each class (here p = kbinompois).
#
diffpolyapois.nll <- function(D,s){
  #  s <- exp(lns)  #Sometimes helpful to constrain positive
  kbinompois <- Diff(c(1,rep(0,5)),1,D) 
  lnlik <- apply(datamat,1,dpolya,alpha=s*kbinompois,log=TRUE)
  return(-sum(lnlik))
}

llfit <- mle( diffpolyapois.nll,start=list(D=0.3,s=30),lower=0.01,method='L-BFGS-B')
summary(llfit)
phat <- coef(llfit)
#logLik(llfit)
#ci <- confint(llfit) #The profile fails on s. Need to do by hand.
#ci
#AIC(llfit)
set.seed(4) #4 for science paper
matplot(1:6,t(datamat),type="n",xlab="Patch",ylab="Abundance")
#Overlay simulation
Ninit <- rowSums(datamat)
kbinompois <- Diff(c(1,rep(0,5)),1,phat[1]) #Kernel can be outside loop
for (i in 1:24){
  n <- sample(Ninit,1,replace=TRUE) #Bootstrap n
  tmp <- rdirichlet(1,phat[2]*kbinompois)
  sim <- rmultinom(1,n,tmp)
  lines(jitter((1:6),0.2),sim,col="red")
}
matpoints(1:6,t(datamat),pch=16,col="blue")
lines(1:6,Diff(c(182,rep(0,5)),1,phat[1]),col="black",lwd=4)
legend('topright',c("Data","Kernel","Simulation"),lty=1,
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



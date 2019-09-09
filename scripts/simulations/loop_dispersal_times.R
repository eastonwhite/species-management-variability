
par(mfrow=c(1,3))
#all_data=NULL
for (dispersal_times in c(24,6,1)){
  
  
  
dTime=dispersal_times  

dispersal_subset <- get(paste('dispersal',dTime,sep=''))

mydatawide <- reshape(dispersal_subset[,c('landscape','patch','count')],
                      idvar=c("landscape"),
                      timevar="patch",direction="wide")
datamat <- as.matrix(mydatawide[,-(1)])
rownames(datamat) <- colnames(datamat) <- NULL
#datamat = cbind(datamat,rep(0,12))
datamat= round(datamat)

set.seed(4) #4 for science paper
matplot(1:6,t(datamat),type="n",xlab="Patch",ylab="Abundance",las=1,main=paste('dispersal time = ',dTime,sep=''),ylim=c(0,150))
#Overlay simulation
Ninit <- rowSums(datamat)
kbinompois <- Diff(c(1,rep(0,5)),dTime/48,phat[1]) #Kernel can be outside loop
for (i in 1:24){
  n <- sample(Ninit,1,replace=TRUE) #Bootstrap n
  tmp <- rdirichlet(1,phat[2]*kbinompois)
  sim <- rmultinom(1,n,tmp)
  lines(jitter((1:6),0.2),sim,col="red")
}
matpoints(1:6,t(datamat),pch=1,col="blue")
#lines(1:6,Diff(c(colMeans(datamat)[1],rep(0,5)),1,phat[1]),col="black",lwd=4)
legend('topright',c("Data","Simulation"),
       col=c("blue","red"),pch=c(0,15),lty=1)


#all_data = rbind(all_data,datamat)
}

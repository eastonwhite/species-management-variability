source("DiffvarD_LAint.r")

N <- c(1,rep(0,11))
N0.5 <- DiffvarD_LAint(ls=N,tt=0.5,Dbar=0.5,k=1e10)
N1 <- DiffvarD_LAint(ls=N,tt=1,Dbar=0.5,k=1e10)
N2 <- DiffvarD_LAint(ls=N,tt=2,Dbar=0.5,k=1e10)

plot(1:12,N,xlab="Patch")
points(1:12,N0.5, col=2)
points(1:12,N1, col=3)
points(1:12,N2, col=4)
legend("topright",legend=paste("Time=",c(0,0.5,1,2)),pch=1,col=1:4)


source("Diff.r") #The "Poisson" diffusion model

N <- c(1,rep(0,11))
N0.5 <- Diff(N,t=0.5,D=0.5)
N1 <- Diff(N,t=1,D=0.5)
N2 <- Diff(N,t=2,D=0.5)

plot(1:12,N,xlab="Patch")
points(1:12,N0.5, col=2)
points(1:12,N1, col=3)
points(1:12,N2, col=4)
legend("topright",legend=paste("Time=",c(0,0.5,1,2)),pch=1,col=1:4)


source("DiffSt.r") #Poisson simulator

N <- c(100,rep(0,11))
N0.5 <- DiffSt(N,tt=0.5,D=0.5)
N1 <- DiffSt(N,tt=1,D=0.5)
N2 <- DiffSt(N,tt=2,D=0.5)

plot(1:12,N,xlab="Patch")
points(1:12,N0.5, col=2)
points(1:12,N1, col=3)
points(1:12,N2, col=4)
legend("topright",legend=paste("Time=",c(0,0.5,1,2)),pch=1,col=1:4)





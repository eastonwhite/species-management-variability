# Dispersal by diffusion (deterministic)
# dNx/dt = D(Nx+1 - 2Nx + Nx-1)
# or DD version:
# dNx/dt = f(Nx+1).Nx+1 - 2f(Nx).Nx + f(Nx-1).Nx-1
#
library(deSolve)

Diff <- function(y,t,D,theta_D=0) {
  p <- c(D,theta_D)
  diffuse <- function( t, y, p ) {
    np <- length(y)
    Dx <- p[1:(length(p)-1)]
    thetaD <- p[length(p)]
    fN <- Dx * y^thetaD #Density dependent dispersal function, theta_D = 0 is DI        
    dN <- fN * y  #Number dispersing per hole
    holes <- c(1,rep(2,(np-2)),1) #Number of holes
    dNfrR <- c(dN[2:np],0)        #Number disperse from right
    dNfrL <- c(0,dN[1:(np-1)])    #Number disperse from left
    #cbind(y,dN,dNfrR,dNfrL)
    list(dNfrR + dNfrL - holes * dN)
  }
  out <- lsoda(y,c(0,t),diffuse, p, rtol=1e-4, atol=1e-6)
  out[2,2:ncol(out)]
}


#There is a faster way using a LA solution. See DiffvarD_LAint.r.
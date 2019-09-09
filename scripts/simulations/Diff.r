# 9 Oct 2013
# Added the Diff_LA function, moved the original Diff function to Diff_solver
# and cleaned it up a bit.

# Solves the diffusion equation for a dispersal period (e.g. 24 hours)

# Dispersal by diffusion (deterministic)
# dNx/dt = D(N[x+1] - hN[x] + N[x-1])
# or density dependent (DD) version:
# dNx/dt = f(N[x+1])N[x+1] - hf(N[x])N[x] + f(N[x-1])N[x-1]
# or patch specific dispersal version
# dNx/dt = D[x+1]N[x+1] - hD[x]N[x] + D[x-1]N[x-1]
# where h is the number of holes (usually 2, except for end patches = 1)

# Arguments
# N - a vector of abundances (N[x]) representing the landscape
# tt - the total length of time of the dispersal period (e.g. 2 hr)
# D  - the diffusion coefficient, either a scalar or vector. The vector version
#      allows D to vary by patch, i.e. D[x].
# theta_D is not properly implemented yet

#If D is scalar and the landscape is smallish (<97), use the linear algebra
#solution, otherwise use the solver.
Diff <- function( N, tt, D ) {
    if ( length(D) == 1 & length(N) < 97 ) {
        return( Diff_LA(N,tt,D) )
    } else {
        return( Diff_solver(N,tt,D) )
    }
}


# Diff_solver()
# This version of diffusion dispersal uses the ODE solver in R to solve the 
# system of ODEs and was used in the Science 2009 paper. This is fairly slow. It
# might be speeded up by coding the diffuse function in C (see help for lsoda).
# It might also be speeded up using a banded form for the solver. See PDE code
# for Berestycki model.

library(deSolve)

Diff_solver <- function( N, tt, D ) {
    diffuse <- function( t, N, D ) {
        np <- length(N)
        dN <- N * D  #Number dispersing per hole
        holes <- c(1,rep(2,(np-2)),1) #Number of holes
        dNfrR <- c(dN[2:np],0)        #Number disperse from right
        dNfrL <- c(0,dN[1:(np-1)])    #Number disperse from left
        list(dNfrR + dNfrL - holes * dN)
    }
    out <- lsoda( N, c(0,tt), diffuse, D, rtol=1e-4, atol=1e-6 )
    out[2,2:ncol(out)]
}

# Diff_LA() We can formulate this as a linear ODE system and solve using linear 
# algebra. The matrix A specifies the system. Currently this only does scalar D.
# I think we could do this for vector D - need to look at the math. This is
# faster than Diff_solver when landscape is small (e.g. 10) but becomes slower
# with increasing landscape size (e.g. 100).
# 
Diff_LA <- function( N, tt, D ) {
    np <- length(N) #number of patches
    A <- matrix(c(-1,1,rep(0,np-2),rep(c(1,-2,1,rep(0,np-2)),np-2),1,-1),np,np)
    eig <- eigen(A);
    const <- solve(eig$vectors,N); #N is initial conditions
    return( eig$vectors %*% (const*exp(tt*eig$values*D)) )
}



#This is a prototype for density dependent dispersal that may or may not work
#properly.
DD_Diff <- function( N, tt, D, theta_D=0 ) {
    parms <- c(D,theta_D)
    diffuse <- function( t, N, parms ) {
        np <- length(N)
        Dx <- parms[1:(length(parms)-1)]
        thetaD <- parms[length(parms)]
        fN <- Dx * N^thetaD #Density dependent dispersal function, theta_D = 0 is DI        
        dN <- fN * N  #Number dispersing per hole
        holes <- c(1,rep(2,(np-2)),1) #Number of holes
        dNfrR <- c(dN[2:np],0)        #Number disperse from right
        dNfrL <- c(0,dN[1:(np-1)])    #Number disperse from left
        list(dNfrR + dNfrL - holes * dN)
    }
    out <- lsoda( N, c(0,t), diffuse, parms, rtol=1e-4, atol=1e-6 )
    out[2,2:ncol(out)]
}

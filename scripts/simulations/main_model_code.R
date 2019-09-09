#Default is NBBg Ricker model ("RickerStSexBS_DEhB") described in Nature 2008 paper
#together with stochastic model of dispersal ("DiffSt_Polya") described
#in Science 2009 paper.
#
# Adapted from Brett Melbourne 19 April 2009
# Created by Easton White
# Last updated 2-Sep-2019


source("Ricker.r")              #Deterministic Ricker
source("RickerStBS.r")          #Stochastic Rickers:
source("RickerStBS_DhB.r")
source("RickerStBS_EhB.r")
source("RickerStBS_DEhB.r")
source("RickerStSexBS.r")
source("RickerStSexBS_DhB.r")
source("RickerStSexBS_EhB.r")
source("RickerStSexBS_DEhB.r")
source("Diff2.r")            #Deterministic dispersal, diffusion
source("DiffSt.r")          #Stochastic dispersal, diffusion
source("DiffvarD_LAint.r")  #Deterministic dispersal with individual variance in D
source("DiffvarDSt_IBM.r")  #Stochastic dispersal with individual variance in D
source("DiffSt_Polya.r")    #As for DiffSt + Dirichlet overdispersion
source("Distrav.r")
source("DiffSt_single_patch.R")
source("harvesting_functions.R")


# This creates a function where we can use different levels of dispersal and harvesting rate
simulate_model = function(dispersal_rate,harvesting_rate,s_param,disp_time,R,alpha,kRd,kRe,kDd,np,nl,ng){
  
  Ninit <- c(20,rep(0,times=np-1))            #Number of adults to start in all patches  

D <- dispersal_rate         # dispersal rate (fraction of beetles leaving patch).
s <- s_param             #dispersion parameter for Dirichlet.
H = harvesting_rate             # Harvesting rate outside of reserves   

spatial_structure =  rep(1,times=np)# MPA spatial structure (vector of length np). 0 designates reserve. 1 is not reserve

#Stochastic growth settings
GrType <- "RickerStSexBS_DEhB"     
                       #"Ricker" = deterministic
                       #"RickerStBS" = stochastic birth (B) and survival (S)
                       #"RickerStBS_DhB" = stochastic B + S + demographic heterogeneity in birth (DhB)
                       #"RickerStBS_EhB" = stochastic B + S + environmental heterogeneity in birth (DhB)
                       #"RickerStBS_DEhB" = stochastic B + S + dem & env heterogeneity in birth (DhB)
                       #"RickerStSexBS" = stochastic B + S + sex ratio
                       #"RickerStSexBS_DhB" = stochastic B + S + sex ratio + DhB
                       #"RickerStSexBS_EhB" = stochastic B + S + sex ratio + EhB
                       #"RickerStSexBS_DEhB" = stochastic B + S + sex ratio + DhB + EhB

#Stochastic dispersal settings
DisType <- "DiffSt_Polya"#"DiffSt_single_patch"    
                       #"Diff" = deterministic diffusion.
                       #"DiffSt" = stochastic (Poisson) diffusion.
                       #"DiffvarD_LAint" = deterministic diffusion + demographic heterogeneity in D.
                       #"DiffvarDSt_IBM" = stochastic diffusion + DhD (negative binomial)
                       #"DiffSt_Polya" = stochastic (Poisson) diffusion with overdispersion.

HarvestType <- "fixed_proportion" #fixed_quota, #fixed_proportion_target_high, #fixed_proportion_target_high

#----Main Program--------------------------------------------------------------


output=NULL

#Replicate landscapes
for (ll in 1:nl) {
  LS=NULL
  Yield=NULL
  

# Initialize landscapes  
  ls <- 0 * 1:np 
  ls[1:np] <- Ninit #Seed with adults (male+female)

# Spread through a landscape
  for (g in 1:ng){
    dtt <- proc.time() #Record time of loop (for graphic display)

  # Grow locally with environmental and demographic stochasticity
    ls <- 
    switch( GrType,
      Ricker = Ricker(ls,R,alpha),
      RickerStBS = RickerStBS(ls,R,alpha),
      RickerStBS_DhB = RickerStBS_DhB(ls,R,alpha,kRd),
      RickerStBS_EhB = RickerStBS_EhB(ls,R,alpha,kRe),
      RickerStBS_DEhB = RickerStBS_DEhB(ls,R,alpha,kRd,kRe),
      RickerStSexBS = RickerStSexBS(ls,R,alpha),
      RickerStSexBS_DhB = RickerStSexBS_DhB(ls,R,alpha,kRd),
      RickerStSexBS_EhB = RickerStSexBS_EhB(ls,R,alpha,kRe),
      RickerStSexBS_DEhB = RickerStSexBS_DEhB(ls,R,alpha,kRd,kRe),
    )

  # Dispersal with demographic stochasticity
    ls <- 
      switch( DisType,
              Diff = Diff(ls,disp_time/48,D),
              DiffSt = DiffSt(ls,disp_time/48,D),
              DiffvarD_LAint = DiffvarD_LAint(ls,disp_time/48,D,kD),
              DiffvarDSt_IBM =DiffvarDSt_IBM(ls,disp_time/48,D,kD),
              DiffSt_Polya = DiffSt_Polya(ls,disp_time/48,D,s_param)
      )
    
    
    #harvest_rates=spatial_structure*H
    #yield = harvest_rates*ls
    yield <- switch(HarvestType,
                   fixed_proportion = fixed_proportion(ls,H),
                   fixed_quota = fixed_quota(ls,H),
                   fixed_proportion_target_high = fixed_proportion_target_high(ls,H),
                   fixed_proportion_target_low = fixed_proportion_target_low(ls,H))
    
    ls <- ls - yield
    
  
    
    Yield = cbind(Yield,yield)
    LS = cbind(LS,ls)

 # Model outputs
    output = rbind(output,cbind(rep(H,np),rep(disp_time,np),rep(ll,np),1:6,g,ls,yield))
  }#End landscape

# output[ll,1:5]= c(ll,sum(LS),sum(Yield),H,D)


 #if (ShowSim == 1) {print(paste(ll,sum(ls),sum(yield)))}
 
 
}#End replicates. All landscapes done.

return(as.data.frame(output))
}

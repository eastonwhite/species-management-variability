# Created by Easton R. White (eastonrwhite@gmail.com)
# Last edited: 2-Sep-2019
# BeetleControl project

# This code runs the main model script, but allows for changes to D (dispersal) and H (harvesting)

# Default model parameters
  #Simulation setup
  np <- 6            #Number of patches in a landscape
  nl <- 1000               #Number of landscapes
  ng <- 6               #Number of generations
  
  #Population parameters
  R <- 3.53#3.53             #R for Ricker (births)
  alpha <- 0.0020#0.0020      #Competition/cannibalism parameter 
  kRd <- 2.93#2.93            #k (shape parameter of gamma) for demographic heterogeneity in R
  kRe <- 26.74#26.74           #k (shape parameter of gamma) for environmental heterogeneity in R
  kDd <- 100000               #k variation in D (shape parameter of gamma) between individuals
  

# Define parameter vectors to try
#D_vector = 0.24#c(0.01,0.1,0.3)#seq(0.01,0.9,by=0.3)
disp_time_vector = c(1,6,24)#seq(1,24,1)#c(24,6,1) #seq(1,20,2)#
H_vector=seq(0,0.8,0.05)#0.75#c(0,0.25,0.5,0.75)#seq(0,0.89,by=0.3)
parameters <- expand.grid(H_vector,disp_time_vector)
names(parameters) = c('hRate','dTime')
parameters$D_value = 0.3844#c(rep(0.27,4),rep(0.42,4),rep(0.72,4),rep(0.22,4),rep(0.95,4)) #rep(0.24,nrow(parameters))#
parameters$s_value = 46.3253#c(rep(127,4),rep(37,4),rep(17,4),rep(39000,4),rep(306,4))#rep(8.00,nrow(parameters))#

model_output <- NULL

source('main_model_code.R')
source('harvesting_functions.R')
for (sim_num in 1:nrow(parameters)){
#  model_output = rbind(model_output,simulate_model(parameters$D_value[sim_num],parameters$hRate[sim_num],parameters$s_value[sim_num],parameters$dTime[sim_num],R,alpha,kRd,kRe,kDd,np,nl,ng))
  model_output = rbind(model_output,cbind(simulate_model(parameters$D_value[sim_num],parameters[sim_num,1],parameters$s_value[sim_num],parameters[sim_num,2],R,alpha,kRd,kRe,kDd,np,nl,ng),R))
  print(sim_num)
}

names(model_output) <- c('hRate','dTime','landscape','patch','gen','count','hcount','param')

#write.csv(x = model_output,file = "../../model_outputs/simulation_results2.csv",row.names = F)
write.csv(x = model_output,file = "../../model_outputs/simulation_results_spatial_spread2.csv",row.names = F)


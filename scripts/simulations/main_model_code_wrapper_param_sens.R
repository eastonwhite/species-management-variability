# Created by Easton R. White (eastonrwhite@gmail.com)
# Last edited: 2-Sep-2019
# BeetleControl project

# This code runs the main model script, but allows for changes to D (dispersal) and H (harvesting)

# Default model parameters
#Simulation setup
np <- 6            #Number of patches in a landscape
nl <- 50               #Number of landscapes
ng <- 6               #Number of generations

param_to_test = 'alpha'
vec_length = 100

#Population parameters
R <- rep(3.53,vec_length)#seq(1,4,length.out = vec_length)#3.5#3.5             #R for Ricker (births)
alpha <- seq(0.0001,0.01,length.out = vec_length)#rep(0.0020,vec_length)#0.003731      #Competition/cannibalism parameter 0.003
kRd <- rep(2.9303,vec_length)#1.07            #k (shape parameter of gamma) for demographic heterogeneity in R
kRe <- rep(26.7391,vec_length)#17.62           #k (shape parameter of gamma) for environmental heterogeneity in R
kDd <- rep(10000,vec_length)               #k variation in D (shape parameter of gamma) between individuals

H <- rep(0.5,vec_length)
D <- rep(0.3844,vec_length)
s <- rep(46.3253,vec_length)
dTime=rep(24,vec_length)

model_output <- NULL

source('main_model_code.R')
source('harvesting_functions.R')
for (sim_num in 1:vec_length){
  model_output = rbind(model_output,cbind(simulate_model(D[sim_num],H[sim_num],s[sim_num],dTime[sim_num],R[sim_num],alpha[sim_num],kRd[sim_num],kRe[sim_num],kDd[sim_num],np,nl,ng),rep(alpha[sim_num],np*nl*ng)))
  print(sim_num)
}

names(model_output) <- c('hRate','dTime','landscape','patch','gen','count','hcount','param')

#write.csv(x = model_output,file = "../../model_outputs/parameter_sensitivity_alpha.csv",row.names = F)


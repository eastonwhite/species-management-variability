# Created by Easton R White
# Last edited 5-Oct-2017
# Harvesting functions


# This is a set of different harvesting functions that could be used in simulations. The default function is the fixed_proportion function


fixed_proportion = function(abundances, harvest_rates){
  yield = round(abundances*harvest_rates)
  #post_harvest_abundances = abundances - yield
  return(yield)
}


fixed_quota = function(abundances, harvest_quota){
  if(length(harvest_quota)==1){yield = rep(harvest_quota,times=length(abundances))}
  yield[which(abundances-harvest_quota<0)] = abundances[which(abundances-harvest_quota<0)]
  return(yield)
}


fixed_proportion_target_high = function(abundances, harvest_rates){
  yield = rep(0,times=length(abundances)) 
  target_patch = which(abundances==max(abundances))
  if(length(target_patch)>1){target_patch = sample(target_patch,1)}
  yield[target_patch] = round(abundances[target_patch]*harvest_rates)
  
  return(yield)
}



fixed_proportion_target_low = function(abundances, harvest_rates){
  yield = rep(0,times=length(abundances)) 
  target_patch = which(abundances==min(abundances[abundances>0]))
  if(length(target_patch)>1){target_patch = sample(target_patch,1)}
  yield[target_patch] = round(abundances[target_patch]*harvest_rates)
  
  return(yield)
}

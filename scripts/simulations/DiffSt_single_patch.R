
# This script defines the 'DiffSt_single_patch' function. The function takes a vector of population sizes and a dispersal parameter D. D is the fraction of individuals dispersing from each patch into nearby patches. Then, the function draws from a multinomial function to distribute dispersing beetles to nearby patches. Importantly, the current iteration of the model only allows beetles to move one patch right or left, or to remain in their natal patch. 


# Created by Easton R. White
# Last edited 14-Jun-2017

DiffSt_single_patch <- function(N,D) {
  
  # define matrix of dispersal probabilities (will be sparse in simple case)
  DM = matrix(0,nrow=length(N),ncol=length(N))
  diag(DM) = 1-D
  diag(DM[-nrow(DM),-1]) = D
  diag(DM[-1,-ncol(DM)]) = D
  
  # Create object which tracks number of dispersers from each patch
  dispersers = matrix(0,nrow=length(N),ncol=length(N))
  for (q in 1:length(N)){
    dispersers[,q]=rmultinom(1,size=N[q],prob=DM[,q])
  }
  
  return(rowSums(dispersers)) #total numbers of beetles in each patch after dispersal
}






get_mean_edge_weights <- function(mcmc_res_var, start_index = NULL) {
  
  
  ## Get the list with the edge weight samples.
  betas = mcmc_res_var$MCMC.chain$betas
  nr_samples = length(betas)
  
  if(is.null(start_index))
    ## Specify from what index position to take the samples from.
    start_index = floor(nr_samples / 2)
  
  
  ## Extract number of elements in each edge weight record.
  ##  This contains the segment identifier (1st element), the edge weights
  ##  of the parents and the bias and SAC node (the last two elements). 
  nr_elements = length(betas[[1]][1,])
  
  
  ## 
  ## Each element in 'betas' contains the edge weight samples for 
  ##  the different segments. First average over these segments, and
  ##  later average over all samples.
  ## 
  
  beta_tmp = rep(0, nr_elements)
  
  
  ## Loop over the chain samples
  for(i in start_index:nr_samples) {
    
    ## Sum up the weights aver taking the mean over the segments.
    beta_tmp = beta_tmp + apply(betas[[i]], 2, mean)
    
  }
  
  ## Exclude the last one or two elements, which correspond to the bias and SAC node (if included).
  beta_tmp = beta_tmp[1:(nr_elements - mcmc_res_var$Grid.obj$additional.parents)]


  ## Exclude the first element, which correspond to the segment ID.
  beta_tmp = beta_tmp[-1]
  
  
  ## Normalize.
  mean_beta = beta_tmp / (nr_samples - start_index + 1)
  
  
  return(mean_beta)
  

  
}


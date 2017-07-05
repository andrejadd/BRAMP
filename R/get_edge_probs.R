
##
## Calculate the edge probabilities for each parent. This is done by averaging 
##  over the binary edge indicators of the chain samples, starting from the middle
##  of the chain. 
##
get_edge_probs <- function(mcmc_res_var, start_index = NULL) {

  
  ## Get the edge samples and number of samples
  edge_samples = mcmc_res_var$MCMC.chain$Structsamples$struct
  nr_samples = length(edge_samples)
  
  if(is.null(start_index))
    ## Specify from what index position to take the samples from.
    start_index = floor(nr_samples / 2)
  
  
  ## Number of edges in each samples. 
  nr_edges = length(edge_samples[[start_index]])
  
  
  ## Sum up all the binary indicators.
  sum_edges = rep(0, nr_edges)
  
  
  ## Loop over chain iterations.
  for (i in start_index:nr_samples) {
    
    sum_edges = sum_edges + edge_samples[[i]]
  }
  
  
  ## Average.
  edge_probs = sum_edges / (nr_samples - start_index + 1)

  
  ## Exclude the last one or two edges, which correspond to the bias and SAC node (if included).
  edge_probs = edge_probs[1:(nr_edges - mcmc_res_var$Grid.obj$additional.parents)]
  
  
  return(edge_probs)
  
}
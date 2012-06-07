# l = length of time series
# min_phase_length = minimum length of each phase
# k_bar = maximum number of phases
# q = number of nodes
simulate_network <- function(l=100, min_phase_length=10, k_bar=10, q=11, 
  lambda_2=0.45, noise=1, net=NULL) {

  if(is.null(net)) {
    net = generate_network(lambda_2, q, min_phase_length, k_bar)
  }
  
  network = net$network;
  epsilon = net$epsilon;
  k = net$k;
  
  # Simulate data from network
  sim_data = matrix(0, q, l) 
  sim_data[,1] = matrix(rnorm(q), q, 1);
  
  matrix_offset = 0;
  begin = 2;
  
  for (i in 1:(k+1)) {
	
	parent_set = network[[i]];
	change = epsilon[[i]];
	
    for (j in begin:change) {
	  new_pt = t(parent_set) %*% sim_data[,j-1]; 	
	  
	  # Nodes without parent are drawn from Gaussian
	  new_pt[new_pt == 0] = rnorm(sum(new_pt == 0));
	  
	  # Add noise
	  #new_pt = new_pt + matrix(rnorm(q, 0, 0.25), q, 1);
	  #new_pt = matrix(rnorm(q, 0, 1), q, 1);
	  new_pt = new_pt + matrix(rnorm(q, 0, noise), q, 1);
	  
	  # Scale to preserve variance = 1
	  sim_data[, j] = new_pt / sqrt(1 + noise*noise);  
	}
	
	begin = change + 1;
  }
  
  return (list(sim_data=sim_data, epsilon=epsilon, k=k, network=network));
  
}

generate_network <- function(lambda_2, q, min_phase_length, k_bar) { 
  # Draw hyperparameters
  #lambda_1 = rgamma(1, 0.6, 0.1)
  lambda_1 = 6  
  #lambda_2 = rgamma(1, 0.75, 0.25)
  #lambda_2 = 0.45


  k = k_bar + 1;
  epsilon = c();

  # Choose number and location of change points
  while(k > k_bar || any(c(epsilon, l) - c(0, epsilon) < min_phase_length)) {
    k = rpois(1, lambda_1);
    epsilon = sort(sample(1:l, k, replace=FALSE))
  }

  # Draw initial network from prior
  parents = matrix(FALSE, q, q);

  for(i in 1:q) {

  num_parents = q+1;

  while(num_parents > q) {
    num_parents = rpois(1, lambda_2); 
  }

  parents_i = 1:q %in% sample(1:q, num_parents, replace=FALSE) 
	  parents[i, parents_i] = TRUE; 
  }

  parents = t(parents);

  network = list();

  # For each phase, draw number of changes and changed edges with 
  # respect to previous phase, then record new network  
  for(i in 1:(k+1)) {  
    change_num = 0;

  while(change_num == 0) {
    change_num = rpois(1, lambda_2);
  }
  

    add_edges = sum(runif(change_num) > 0.5);
    
       edge_num = sum(parents);
    non_edge_num = q*q - edge_num;  
    
    new_parents = parents;
    
    if(add_edges > 0) {
      add_changes = 1:non_edge_num %in% sample(1:non_edge_num, add_edges, replace=FALSE)
         edges = parents[parents == 0];
      edges[add_changes] = 1;
      new_parents[parents == 0] = edges;
    }
    
    if(change_num - add_edges > 0) {
         remove_changes = 1:edge_num %in% sample(1:edge_num, change_num - add_edges, replace=FALSE)
         non_edges = parents[parents > 0];
      non_edges[remove_changes] = 0;
      new_parents[parents>0] = non_edges;
    }
             
    parents = new_parents;
    parent_num = colSums(parents)
    parent_num[parent_num == 0] = 1;
    
    # Scale weights to preserve variance = 1
       network[[i]] = parents / t(matrix(sqrt(parent_num), q, q))
    
    e_v = eigen(network[[i]]);
    
    while(any(abs(e_v$values) > 1)) {
      edge_num = sum(parents); 
      remove_changes = 1:edge_num %in% 
        sample(1:edge_num, 1, replace=FALSE)
      non_edges = parents[parents > 0];
      non_edges[remove_changes] = 0;
      new_parents[parents>0] = non_edges;
      
      parents = new_parents;
      parent_num = colSums(parents)
      parent_num[parent_num == 0] = 1;
    
      # Scale weights to preserve variance = 1
         network[[i]] = parents / t(matrix(sqrt(parent_num), q, q))
      e_v = eigen(network[[i]]);
      
    }
  }
  
  epsilon = c(epsilon, l);
 
  return(list(network=network, epsilon=epsilon, k=k))
}

  
simulate_network_simple <- function(l=50) {  
  parent_set1 = matrix(0, 10, 10);
  parent_set1[c(1, 2, 3), 5] = 1/sqrt(3); 
  parent_set1[9, 6] = 1;
  parent_set1[6, 10] = 1;
  parent_set1[10, 9] = 1;
  
  #min_weight = 0.5;
  #weighting = matrix(rnorm(parent_set1), 10, 10);
  #weighting[abs(weighting) < min_weight] = sign(weighting[abs(weighting) < min_weight]) * min_weight; 
  #parent_set1 = parent_set1 * weighting;
  
  parent_set2 = matrix(0, 10, 10);
  parent_set2[c(1, 2), 5] = 1/sqrt(2); 
  parent_set2[6, 10] = 1;
  parent_set2[10, 9] = 1;
  parent_set2[9, 4] = 1;
  
  #parent_set2 = parent_set2 * weighting;
  
  network = list(parent_set1, parent_set2);
  
  sim_data = matrix(0, 10, l) 
  sim_data[,1] = matrix(rnorm(10), 10, 1);
  
  matrix_offset = 0;
  begin = 2;
  
  half_l = round(l/2);
  
  for (parent_set in network) {
	
    for (i in begin:half_l) {
	  new_pt = t(parent_set) %*% sim_data[,matrix_offset + i-1]; 	    
	  new_pt[new_pt == 0] = rnorm(sum(new_pt == 0));
	  new_pt = new_pt + matrix(rnorm(10, 0, 1), 10, 1);
	  sim_data[, matrix_offset + i] = new_pt / sqrt(2);
	    
	}
	
	matrix_offset = matrix_offset + half_l;
	begin = 1;
  }
  return (sim_data);

}

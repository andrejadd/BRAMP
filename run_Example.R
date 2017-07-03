
source("BRAMP.R")
source("Code/get_edge_probs.R")

## Target node for which to calculate the parent probabilities
target_node = 1

## Input data file and result file
input_data = "Data/Example_Data.Rdata"
result_file = paste("Results/Example_output_target", target_node, ".Rdata", sep="")


## 
## Load input data:
##
##   data_mat  : A m-by-n matrix where m is the number of observations (rows) and n is the number of nodes (columns).
##   xlocs     : Number of locations (observations) along the x-axis.
##   ylocs     : Number of locations (observations) along the y-axis.
##               (Note: xlocs * ylocs must be equal to 'm')
## Optional:
##
##   SAC.nodes : a n-by-m matrix that contains a spatial autocorrelation node for each node (n) and observation (m)  
##               If this is not provided but ENABLE.SAC is True, the function spatAutoCorrelation() below will attempt to calculate it.
##
load(input_data)


## Extract the target values and scale it.
y = as.vector(scale(data_mat[,target_node]))


## Create design matrix, excluding the target node (no self-loop).
X = scale(data_mat[,-target_node])


## 
## Run BRAMP given the input data file, for a target node 'target' 
##  and for 'niter' number of iterations. 
##
mcmc_result = BRAMP(y, X, xlocs, ylocs, target = target_node, nr_iterations = 1000)
save(file=result_file, "mcmc_result")


cat("\n")
##
## Run again but longer. By passing 'result.file' the simulation will
##  pick up the MCMC chain where it last stopped (for this 'niter' must 
##  larger than the previous run). 
##
mcmc_result = BRAMP(y, X, xlocs, ylocs, target = target_node, nr_iterations = 2000, result.file=result_file)
save(file=result_file, "mcmc_result")


## Calculate edge probabilities from the chain samples.
edge_probs = get_edge_probs(mcmc_result)
cat("\nEdge probabilities:\n")
print(edge_probs)

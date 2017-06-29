
source("BRAMP.R")

## Target node for which to calculate the parent probabilities
target_node = 1

## Input data file and result file
input_data = "Data/Data_id10.Rdata"
result_file = paste("Data_id10_target", target_node, ".Rdata", sep="")



## This loads the R data file. 
##
## Currently the data file contains all kinds of information from the synthetic data 
## creation tool:
##
## [1] "networks"           "nrnodes"           
## [3] "Ymatrix"            "Alist"             
## [5] "Ematrix"            "Bmatrix"           
## [7] "cycles"             "XE"                
## [9] "YE"                 "minPhase"          
## [11] "HOMOGENEOUS_STRUCT" "Description"       
## [13] "node.names"         "modelid"           
## [15] "xlocs"              "ylocs"             
## [17] "totallocs"          "nrphases"          
## [19] "type"              
## 
## However, only the following variables are required: 
##
##   Ymatrix  : A n-by-m matrix where n is the number of nodes (rows) and m is the number of observations.
##   xlocs    : Number of locations (observations) along the x-axis).
##   ylocs    : Number of locations (observations along the y-axis).
##
## Optional:
##   SAC.nodes : a n-by-m matrix that contains a spatial autocorrelation node for each node (n) and observation (m)  
##               If this is not provided but ENABLE.SAC is True, the function spatAutoCorrelation() below will attempt to calculate it.
##
load(input_data)

## extract the target values and scale it
Y = as.vector(scale(Model$Ymatrix[target_node,]))

## Temporary helper vector of potential parent nodes, Default: allow all excluding the target itself (no self-loop)
posTF=c(1:(nrow(Model$Ymatrix)))[-c(target_node)]

## Extract all the parent data, transpose and scale it.
X = scale(t(Model$Ymatrix[posTF,]))




## 
## Run BRAMP given the input data file, for a target node 'target' 
##  and for 'niter' number of iterations. 
##
mcmc_result = BRAMP(Y, X, Model$xlocs, Model$ylocs, target = target_node, nr_iterations = 1000)
save(file=result_file, "mcmc_result")


cat("\n")
##
## Run again but longer. By passing 'result.file' the simulation will
##  pick up the MCMC chain where it last stopped (for this 'niter' must 
##  larger than the previous run). 
##
mcmc_result = BRAMP(Y, X, Model$xlocs, Model$ylocs, target = target_node, nr_iterations = 2000, result.file=result_file)
save(file=result_file, "mcmc_result")

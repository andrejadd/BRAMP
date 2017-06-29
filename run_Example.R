
source("BRAMP.R")

input_data = "Data/Data_id10.Rdata"
result_file = "Results/Data_id10_target1_test.Rdata"

## 
## Run BRAMP given the input data file, for a target node 'target' 
##  and for 'niter' number of iterations. 
##
mcmc_result = BRAMP(input_data, target = 1, nr_iterations = 1000)
save(file=result_file, "mcmc_result")


cat("\n")
##
## Run again but longer. By passing 'result.file' the simulation will
##  pick up the MCMC chain where it last stopped (for this 'niter' must 
##  larger than the previous run). 
##
mcmc_result = BRAMP(input_data, target = 1, nr_iterations = 2000, result.file=result_file)
save(file=result_file, "mcmc_result")

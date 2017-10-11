

run_Example <- function() {
    
    ## Target node for which to calculate the parent probabilities
    target_node = 1
    
    ## Input data file and result file
    input_data = "data/Example_Data.Rdata"
    #result_file = paste("Results/Example_output_target", target_node, ".Rdata", sep="")
    
    
    ## 
    ## Load input data:
    ##
    ##   data_mat   : A m-by-n matrix where m is the number of observations (rows) and n is the number of nodes (columns).
    ##   xlocs      : Number of locations (observations) along the x-axis.
    ##   ylocs      : Number of locations (observations) along the y-axis.
    ##               (Note: xlocs * ylocs must be equal to 'm')
    ## Optional:
    ##
    ##   y_SAC_node : A m-length vector that contains the spatial autocorrelation data for the target.
    ##                Can be calculated with function spatAuroCorrelation(), see below.
    ##
    load(input_data)
    
    
    ## Extract the target values and scale it.
    y = as.vector(scale(data_mat[,target_node]))
    
    ## Calculate SAC node for the target.
    ## Define as 'NULL' if not needed.
    y_SAC_node = spatAutoCorrelation(y, xlocs, ylocs)
    
    
    ## Create design matrix, excluding the target node (no self-loop).
    X = scale(data_mat[,-target_node])
    
    
    
    ## 
    ## Run BRAMP given the input data file, for a target node 'target' 
    ##  and for 'niter' number of iterations. 
    ##
    mcmc_result = BRAMP(y, X, y_SAC_node, xlocs, ylocs, nr_iterations = 1000)
    print(mcmc_result)
    
    #save(file=result_file, "mcmc_result")
    
    
    ## Run without additional spatial autocorrelation node (y_SAC_node == NULL).
    BRAMP(y, X, NULL, xlocs, ylocs, nr_iterations = 1000)
    
    
    ## Run with specific random generator seed, will be saved for later reference.
    ##
    ## Note: BRAMP() will not set the seed if the parameter mcmc_rnd_seed is ignored.
    ##       The user is responsible to generate and set a seed, e.g. with 
    ##             set.seed(as.numeric(Sys.time())) 
    ##
    BRAMP(y, X, NULL, xlocs, ylocs, nr_iterations = 1000, mcmc_rnd_seed = 1)
    cat("\n")
    
    
    ##
    ## Continue a previous MCMC simulation by specifying the file name of an old result file.
    ##  This will continue the chain if 'nr_iterations' is greater than in the result file.
    ##
    mcmc_result = BRAMP(y, X, y_SAC_node, xlocs, ylocs, nr_iterations = 2000, result.file=result_file)
    #save(file=result_file, "mcmc_result")
    
    
    ## Calculate edge probabilities from the chain samples.
    edge_probs = get_edge_probs(mcmc_result)
    cat("\nEdge probabilities:\n")
    print(edge_probs)
    
    
    ## Calculate the mean edge weights from the chain samples.
    ##  Note that edge weights between different segments will differ and 
    ##  segment number of locations will change along chain iteration.
    ##
    mean_betas = get_mean_edge_weights(mcmc_result)
    cat("\nMean edge weights (regression coefficients):\n")
    print(mean_betas)
    
    
    ##
    ## Define certain edges to be fixed and run the simulation. "Fixed" means they will
    ##   be always included as incoming nodes and never be removed in one of the edge moves.
    ##
    ##   The fixed edges are defined with the function parameter vecotr 'fixed_edges' that 
    ##   contains the index values of the nodes in the design matrix X. 
    ##   E.g. in order to set the 2nd and 4th node fixed, define the following:
    ##
    cat("\nRun with fixed edges..\n")
    fixed_edges = c(2,4)
    mcmc_result = BRAMP(y, X, y_SAC_node, xlocs, ylocs, nr_iterations = 1000, fixed_edges = fixed_edges)
    
    edge_probs = get_edge_probs(mcmc_result)
    
    cat("\nEdge probabilities:\n")
    print(edge_probs)
    
}    


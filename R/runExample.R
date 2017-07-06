
# source("BRAMP.R")

runExample <- function() {

    ## Target node for which to calculate the parent probabilities
    target_node = 1

    ## Load the example data included in this package.
    ##
    ##  The data set has following attributes:
    ##
    ##   data_mat   : A m-by-n matrix where m is the number of observations (rows) and n is the number of nodes (columns).
    ##   xlocs      : Number of locations (observations) along the x-axis.
    ##   ylocs      : Number of locations (observations) along the y-axis.
    ##               (Note: xlocs * ylocs must be equal to 'm')
    ##
    data(example_data)

    data_mat = example_data$data_mat
    xlocs = example_data$xlocs
    ylocs = example_data$ylocs
    

##   y_SAC_node : A m-length vector that contains the spatial autocorrelation data for the target.
##                Can be calculated with function spatAuroCorrelation(), see below.


    if(!exists("example_data")) {
        cat("could see data\n")
        return(-1)
    }


    ## Extract the target values and scale it.
    y = as.vector(scale(data_mat[,target_node]))

    
    ## Calculate the spatial autocorrelation data (SAC) node for the target.
    ##  This is an additional covariate that reflects the influence of surrounding target covariate quantities.
    ##
    ##  y_SAC_node : A m-length vector that contains the spatial autocorrelation data for the target.
    ##
    ##  Note: can be set to NULL or ignored in call to BRAMP()
    y_SAC_node = spatAutoCorrelation(y, xlocs, ylocs)


    ## Create design matrix, excluding the target node (no self-loop).
    X = scale(data_mat[,-target_node])


     
    cat("\n\nRun BRAMP given the input data file, for a target node 'target' and for 'niter' number of iterations.\n\n")
    mcmc_result = BRAMP(y, X, y_SAC_node, xlocs, ylocs, nr_iterations = 1000)
    print(mcmc_result)

    
    
    cat("\n\nContinue a previous MCMC simulation by passing the return variable of BRAMP()\n back to BRAMP() into the parameter 'mcmc_obj'. This will continue the chain if 'nr_iterations' is greater than in 'mcmc_obj'.\n\n")
    mcmc_result = BRAMP(mcmc_obj = mcmc_result, nr_iterations = 2000)


    cat("\n\nRun without additional spatial autocorrelation node (y_SAC_node == NULL).\n\n")
    BRAMP(y, X, NULL, xlocs, ylocs, nr_iterations = 1000)
   

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
    cat("\n\nRun with fixed edges.\n\n")
    fixed_edges = c(2,4)
    mcmc_result = BRAMP(y, X, y_SAC_node, xlocs, ylocs, nr_iterations = 1000, fixed_edges = fixed_edges)

    edge_probs = get_edge_probs(mcmc_result)

    cat("\nEdge probabilities:\n")
    print(edge_probs)


}

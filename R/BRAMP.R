
codePath=paste(getwd(),"/Code/",sep="")

source(paste(codePath,"mcmc_main.R",sep=""))
source(paste(codePath,"edge_moves.R",sep=""))
source(paste(codePath,"mondrian_moves.R",sep=""))
source(paste(codePath,"initEngine.R",sep=""))
source(paste(codePath,"segments.R", sep=""))
source(paste(codePath,"util.R",sep=""))   #requires pseudoinverse
source(paste(codePath,"invGamma.R",sep=""))
source(paste(codePath,"mvrnorm.R",sep=""))
source(paste(codePath,"ginv.R",sep=""))
source(paste(codePath,"convert.R",sep="")) ## need this ??
source(paste(codePath,"Tree.R",sep=""))

source(paste(codePath,"print_bramp.R",sep=""))
source(paste(codePath,"get_mean_edge_weights.R",sep=""))
source(paste(codePath,"get_edge_probs.R",sep=""))
source(paste(codePath,"spatAutoCorrelation.R",sep=""))


##
##
## This is the main BRAMP function. The function arguments are
##
##   Y             : a m-length vector with the number of 'm' observations for the target node.
##   X             : a n-by-m matrix with 'n' nodes and 'm' observations.
##   y_SAC_node    : a m-length vector with the spatial autocorrelation data. Can be left NULL
##                   if not used.
##   xlocs         : an integer defining the number of observations along the x-axis.
##   ylocs         : an integer defining the number of observations along the x-axis.
##                   [Note that (xlocs * ylocs) == m]. 
##   nr_iterations : the number of MCMC iterations to run.
##   chain_thinout : save samples from the chain every 'chain_thinout' iteration.
##   result_file   : A result file from a previous MCMC simulation. If 'nr_iterations' is 
##                   greater than the iterations in the result file, the simulation is continued
##                   until 'nr_iterations'. 
##
BRAMP <- function(y = NULL, 
                  X = NULL, 
                  y_SAC_node = NULL, 
                  xlocs = 0, 
                  ylocs = 0, 
                  nr_iterations = 0, 
                  chain_thinout = 10, 
                  result_file = NULL,
                  mcmc_rnd_seed = NULL,
                  fixed_edges = NULL,
                  edge_fanin = 5
                  ) {  
  
  ##
  ## Make some initial checks.
  ##
  
  
  ## if no result file was specified, check the input data
  if(is.null(result_file)) {
  
      if((xlocs * ylocs) != nrow(X))
      stop("Number of total observation has to match (xlocs * ylocs).")
  
    if(length(y) != nrow(X))
      stop("Length of response vector y has to match row number in design matrix X.")
    
    if(!is.null(y_SAC_node))
      if(length(y) != length(y_SAC_node))
        stop("Length of response vector y has to match length of SAC vector y_SAC_node.")

    if(chain_thinout < 1) 
      stop("Specify valid chain_thinout greater or equal to 1.")
  }
  
  
  if(nr_iterations <= 0)
    stop("Specify valid iteration number greater than 0.")
  
  
  
  
  
  method.name = "BRAMP"
  
  ## Start with this iteration (required when a simulation is continued).
  end.iter = nr_iterations

  
  ## The Mondrian process start budget.
  start.budget = 1
  
  
  ## The start of iteration can be changed through loading an already existing result file.
  start.iter = 1
 
  
  ## Flag that tells us if to proceed the MCMC chain from a previous simulation.
  ##  Do not change, this is done by checks below.
  PROCEED_CHAIN = F
 
 
  ## Reset the main data structures.
  Grid.obj = NULL
  MCMC.chain = NULL
  HYPERvar = NULL

 
  ##
  ## Check if a result file was passed and if it exists.
  ##  If yes, extract the MCMC chain and the corresponding data (which means no data file is needed)
  ##  and proceed running the chain, if the number of iterations ('nr_iteration') passed to BRAMP() is larger 
  ##  than the chain in the result file.
  ##
  
  if(!is.null(result_file)) {
  
    ## Check if 'mcmc_result' was 
    ## Try to load file with 'mcmc_result' data.
    tryCatch({    
      
      ## Load previous MCMC chain, includes the 'mcmc_result' data structure.
    	load(result_file)
      
    }, error = function(e) {
        
      cat("[", method.name, "] Failed to load result file ", result_file, ". Exiting ...\n")
      return(-1)
    })

    
    ## Check if the data structure was loaded.
    ## Do NOT use exists() because it also looks into the global user environment
    ## i.e. outside this function.
    if(!any("mcmc_result" == ls())) {
      cat("[", method.name, "] Failed to load 'mcmc_result' data structure from result file ", result_file, ". Exiting ...\n")
      return(-1)
    }
  
    
    ## Check if the data structure is of class 'bramp'.
    if(class(mcmc_result) != "bramp") {
      cat("[", method.name, "] Loaded 'mcmc_result' data structure is not of class 'bramp' ", result_file, ". Exiting ...\n")
      return(-1)
    }
    

    ## 
    ## Extract necessary data from the previous MCMC run
    ## This also includes the data itself (X and Y) and the hyper-parameters
    ## 
    MCMC.chain = mcmc_result$MCMC.chain
    Grid.obj = mcmc_result$Grid.obj
    HYPERvar = mcmc_result$HYPERvar
    
    
    ## Load the data.
    X = mcmc_result$X
    y = mcmc_result$y
    
    
    ## Extract the number of iterations this chain was run before.  
    iters = MCMC.chain$Structsamples$iter
    last.iter = iters[[length(iters)]]

    
    ## Check if less iterations are in the chain than requested
    ## (The added value of 10 is to make sure we don't start senseless for small amounts of iterations)
    if( (last.iter + 10) < end.iter) {

      ## set the iteration from which to continue
      start.iter = last.iter + 1
        
      ## flag that we will proceed with this MCMC chain. This skips all the initialization stuff below.
      PROCEED_CHAIN = T
	
      cat("[", method.name, "] Proceeding with chain from ", start.iter, " -> ", end.iter, "\n")
        
    } else {
       
      cat("[", method.name, "] Chain in result file has complete iteration count: ", last.iter, " and requested were ", end.iter, ".",sep="")
      return(-1)
      
    } 

  } ## End of reading the result file.

  
  ## 
  ## This is executed if we do not proceed from an old chain.
  ##  It will initialize some variables and load the input data.
  ##
  if(!PROCEED_CHAIN) {

    ## Set random generator seed if provided. Variable 'mcmc_rnd_seed' is 
    ##  saved into MCMC.chain$mcmc_rnd_seed for later reference.
    if(!is.null(mcmc_rnd_seed)) {
      set.seed(mcmc_rnd_seed)
    } 
    
    
    ## Get number of available parent nodes.
    nr.parents = ncol(X)
    
    
    ## Add the bias (intercept) node - this is just a '1' at the last column.
    X = cbind(X,array(1,length(X[,1])))
    
    
    ## Maximum number of parent nodes (fan-in restriction). A low limit is needed 
    ## for birth proposals based on precomupted posterior distribution (method 4),
    ## otherwise you can set smax = q
    smax = min(edge_fanin,nr.parents);
    
    
    ## Minimum number of locations in a segment. Don't make this too small, otherwise there is 
    ## not enough data to make proper calculations. 
    minSegsize = 9
    
    
    cat("[", method.name, "] Number locations: " , xlocs * ylocs, " with x: ",  xlocs, " , y: ",  ylocs, ", start.budget: ", start.budget, ", parent nodes: ", nr.parents, ", fan-in: ", smax, ", minimum segment size: ", minSegsize, "\n")

    
    ## Defines the number of additional entries in the design matrix
    ## A value of 1 stands for the additional bias node. 
    ## If the SAC node, below, is added, it will be added as another additional node.
    additional.parents = 1

    
    ##
    ## Check if the SAC node data was provided to the BRAMP() function call.
    ##   If Yes, add it to the design matrix.
    ##
    if(!is.null(y_SAC_node)) {

        ## Add the SAC node data after it was scaled.
        X = cbind(X, as.vector(scale(y_SAC_node)))
     
        ## Remember the SAC node as additional parent node.
        additional.parents = additional.parents + 1  
      
        cat("[", method.name, "] Spatial autocorrelation (SAC) node added.\n")
      
    } else {
      
      cat("[", method.name, "] No Spatial autocorrelation (SAC) node provided.\n")
      
    }

    
    ## Initialize the hyper-parameters
    ## For the signal-to-noise (SNR) ratio: delta2 ~ IG(alpha.snr,beta.snr)
    alpha.snr = 2
    beta.snr = 0.2
    v0 = 1 
    
    
    ## Inverse Gamma for the SNR. 
    delta.snr = (1 / rgamma(1, shape=alpha.snr, scale=1/beta.snr))

    
    ## Initialize the hyper-parameter variable list.                   
    HYPERvar = list( c = 0.5,              # for edge move proposals.
      alpha.var = v0/2, beta.var = v0/2,   # for weight variance sigma2, v0/2 as in Marcos AISTATs paper
      alpha.snr = alpha.snr,               # alpha parameter for SNR sampling.
      beta.snr = beta.snr,                 # beta parameter for SNR sampling.
      delta.snr = delta.snr,               # the sampled SNR value.
      alphalbd = 1, betalbd = 0.5          # for the parent nodes.
      )
    
    
    ## Initialize MCMC.chain data structure.
    MCMC.chain = list(Structsamples = list(struct = list(), 
                                           iter=list(),
                                           mondrian.tree=list()),
			                  segment_map = list(),             ## Samples of the the Mondrian map of segments.
			                  betas = list(),                   ## Samples of the edge weights for each segment.
                        counters=list(),                  ## Keep track of acceptance and rejection moves.
                        delta.snr = c(),                  ## Samples of the SNR parameter
                        params=matrix(0,nrow=0, ncol=6),
                        chain_thinout = chain_thinout,    ## save the chain samples in every 'chain_thinout' iteration (saves memory) 
			                  iteration_save_betas = 1,         ## start to save betas (edge weights) from iteration 'iteration_save_betas (saves memory)
			                  mcmc_rnd_seed = mcmc_rnd_seed
			)
    

    cat("[", method.name, "] Init.. ")
    
 
    ## Initialisation of first iteration parameters
    Grid.obj = initEngine(
      X = X,
      HYPERvar = HYPERvar,
      additional.parents = additional.parents,
      nr.parents = nr.parents,
      start.budget = start.budget,
      xlocs = xlocs,
      ylocs = ylocs,
      minSegsize = minSegsize,
      smax = smax,
      fixed_edges = fixed_edges
      )
    
    cat("ok\n")
    
  } ## end of !PROCEED_CHAIN, i.e. a new MCMC simulation initialization

  
  ##
  ##  
  ## Execute the main MCMC procedure in main()
  ##
  ##
  cat("[", method.name, "] Starting MCMC chain.. \n")
  
  result = mcmc_main(y, X, MCMC.chain, Grid.obj, HYPERvar, start.iter, end.iter)

  cat("\n---------------------------------------------------\n")
  cat("End of MCMC iterations\n")
  
  
  ##
  ## Evaluate the chain
  ##
  
  result$edge_probs = get_edge_probs(result)
  
  result$mean_coefficients = get_mean_edge_weights(result)
  
  result$call = match.call()
  
  class(result) = "bramp"
  
  return(result)
  
}

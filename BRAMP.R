
codePath=paste(getwd(),"/Code/",sep="")

source(paste(codePath,"moves.R",sep=""))
source(paste(codePath,"main.R",sep=""))
source(paste(codePath,"initEngine.R",sep=""))
source(paste(codePath,"segments.R", sep=""))
source(paste(codePath,"util.R",sep=""))   #requires pseudoinverse
source(paste(codePath,"spatAutoCorrelation.R", sep=""))
source(paste(codePath,"invGamma.R",sep=""))
source(paste(codePath,"mvrnorm.R",sep=""))
source(paste(codePath,"ginv.R",sep=""))
source(paste(codePath,"convert.R",sep="")) ## need this ??
source(paste(codePath,"Tree.R",sep=""))



BRAMP <- function(Y, X, xlocs, ylocs, target, nr_iterations = NULL, ENABLE.SAC=T, result.file = 'none_sense_file', SAC.nodes=NULL) {  
  
  ## remove all but arguments  
  ## rm(list= ls()[ls() != "X"ls()!="indata" && ls()!= "target" && ls()!="nr_iterations" && ls()!="ENABLE.SAC"  && ls()!="result.file"])

  
  method.name = "BRAMP"
  
  ## Start with this iteration (required when a simulation is continued).
  end.iter = nr_iterations

  ## The Mondrian process start budget.
  start.budget = 1
  
  ## The start of iteration can be changed through loading an already existing result file.
  start.iter = 1
 
  ## Unix time (seconds since 1970), should be fine for the cluster, might add/subtract milliseconds 
  set.seed(as.numeric(Sys.time()))

  ## flag that tells if to proceed MCMC chain from previous run
  ## ! Do not change, this is done by checks below !
  PROCEED.CHAIN = F
  RESULT.EXISTS = F

  ## flags result file existence
  result.exists.str = ""

  ## set these here, so check below will not fail (otherwise I would use exists() but even then I have to check for NULL)
  Grid.obj = NULL
  MCMC.chain = NULL

  cat("[", method.name, "] Target node: ", target, "\n")
 
  ##
  ## Check if a result file was passed and if it exists.
  ##  If yes, extract the MCMC chain and the corresponding data (which means no data file is needed)
  ##  and proceed running the chain, if the number of iterations ('nr_iteration') passed to BRAMP() is larger 
  ##  than the chain in the result file.
  ##
  
  if(file.exists(result.file)) {

    ## Try to open file
    tryCatch({    
      
      ## Load previous MCMC chain, includes the 'mcmc_result' data structure.
    	load(result.file)
      
    }, error = function(e) {
        RESULT.EXISTS = F	
	cat("[", method.name, "] Something wrong while trying to open existing results file, starting simulation from scratch.\n.")
    })

    ##
    ## Check if the 'mcmc_result' data structure was loaded
    ##
    if(!is.null(mcmc_result)) {
      
      RESULT.EXISTS = T
    
      ## 
      ## Extract necessary data from the previous MCMC run
      ## This also includes the data itself (X and Y) and the hyper-parameters
      ## 
      MCMC.chain = mcmc_result$MCMC.chain
      Grid.obj = mcmc_result$Grid.obj
      X = mcmc_result$X
      Y = mcmc_result$Y
      HYPERvar = mcmc_result$HYPERvar
      
      ## extract the number of iterations this chain was run before  
      iters = MCMC.chain$Structsamples$iter

      ## this is the last iteration number
      last.iter = iters[[length(iters)]]

      ## check if less iterations are in the chain than requested
      ## (The added value of 100 is to make sure we don't start senseless for small amounts of iterations)
      if( (last.iter + 100) < end.iter) {

        ## set the iteration from which to continue
        start.iter = last.iter + 1
        
        ## flag that we will proceed with this MCMC chain. This skips all the initialization stuff below.
        PROCEED.CHAIN = T
	
        cat("[", method.name, "] Proceeding with chain from ", start.iter, " -> ", end.iter, "\n")
        
      } else {
        
        result.exists.str = paste("and chain has all iterations (", last.iter, " and requested were ", end.iter, ").",sep="")
      
        } 
      
    } 
    
  }


  ## If the result file exists and the chain is complete given the requested number of iterations, exit..
  if(RESULT.EXISTS && !PROCEED.CHAIN) {
  
      cat("[", method.name, "] Exiting simulation: result file exists ", result.exists.str, "\n")
     
      return(1)
      
  }
  
  ## 
  ## This is executed if we do not proceed from an old chain.
  ##  It will initialize some variables and load the input data.
  if(!PROCEED.CHAIN) {

    ## This will hold the actual values of the chain states.
    MCMC.chain = NULL

    ## number of putative parents
    nr.parents = ncol(X)
    
    ## Add the bias (intercept) node; this is just a '1' for the first column.
    X = cbind(X,array(1,length(X[,1])))
    
    ## This variable is used to tell the method not to alter some specific initially set edges during run
    ## It should be NULL if no fixed edges are assumed (except of course the bias and optional the SAC edge)
    ##    otherwise it is a vector with the parent node indices that define the fixed incoming edges
    
    ## check if this makes sense, seems confusing: FIXED.INIT.EDGES will be set to this value and INIT.EDGES.FROM.FILE will be set to a file from where to read these edges - which causes the initilization routine to read the fixed edges from a file.
    FIXED.INIT.EDGES=NULL

    # FIXED.INIT.EDGES=seq(1,12) ## use for the Outer Hebrides data when soil attributes (node 1..12) shall be fixed
   
    
    ## Maximum number of parent nodes (fan-in restriction). A low limit is needed 
    ## for birth proposals based on precomupted posterior distribution (method 4),
    ## otherwise you can set smax = q
    smax = min(8,nr.parents);
    
    ## Minimum number of locations in a segment. Don't make this too small, otherwise there is 
    ## not enough data to make proper calculations. 
    minSegsize = 9
    
    cat("[", method.name, "] Number locations: " , xlocs * ylocs, " with x: ",  xlocs, " , y: ",  ylocs, ", start.budget: ", start.budget, ", parent nodes: ", nr.parents, ", fan-in: ", smax, ", minimum segment size: ", minSegsize, "\n")

  
   

    ## Check if we want spatial autocorrelation (SAC Nodes).
    if(ENABLE.SAC) {
      
      ## Check if SAC nodes already provided.
      if(!is.null(SAC.nodes)) {
        
        ## Add a constant vector to predictor data (representing the bias nodes) and the spatial autocorrelation data of the target node.
        X = cbind(X, as.vector(scale(SAC.nodes[target,])))
        
      } else {
        ## .. otherwise, try to calculate the SAC nodes.
        spatAC = spatAutoCorrelation(Y, xlocs, ylocs)
        
        ## Scale the SAC node and append to the design matrix 'X'.
        X = cbind(X, scale(spatAC))
      }
      
      ## Define the additional incoming edges: the bias and the SAC node
      additional.parents = 2  
      
      cat("[", method.name, "] Enabled spatial autocorrelation (SAC) node\n")
      
    } else {
      
      ## Only define the bias as additional edge, since the SAC node is not included.
      additional.parents = 1 

      cat("[", method.name, "] Disabled spatial autocorrelation (SAC) node\n")

    }

    ## Initialize the hyper-parameters
    ## For the signal-to-noise (SNR) ratio: delta2 ~ IG(alpha.snr,beta.snr)
    alpha.snr = 2
    beta.snr = 0.2
    v0 = 1 
    
    ## Inverse Gamma for the SNR
    delta.snr = (1 / rgamma(1, shape=alpha.snr, scale=1/beta.snr))

                       
    HYPERvar = list( c = 0.5, # for edge move proposals
      alpha.var = v0/2, beta.var = v0/2,   # for weight variance sigma2, v0/2 as in Marcos AISTATs paper
      alpha.snr=alpha.snr, beta.snr=beta.snr, delta.snr = delta.snr,  
      alphalbd = 1, betalbd = 0.5  # for parent nodes
      )
    
    
    
    ## Note, this is just a placeholder, define the file below and look into initEngine to make use
    ## of predefined edges (e.g. from another method like Lasso)
    INIT.EDGES.FROM.FILE = NULL
    
    ## if this is not null it means we can set the file from where to read it
    #if(!is.null(FIXED.INIT.EDGES)) {
    #  INIT.EDGES.FROM.FILE = paste("../Data/", DATA.TYPE, "/Data_id", dataid, "_edgeprobs.Rdata", sep="")
    #}
    
    cat("[", method.name, "] Init.. ")
 
    ## Initialisation of first iteration parameters
     Grid.obj = initEngine(X=X,
      HYPERvar=HYPERvar,
      additional.parents=additional.parents,
      nr.parents=nr.parents,
      start.budget=start.budget,
      xlocs=xlocs,
      ylocs=ylocs,
      minSegsize=minSegsize,
      smax=smax,
      target=target,
      INIT.EDGES.FROM.FILE=INIT.EDGES.FROM.FILE,
      FIXED.INIT.EDGES=FIXED.INIT.EDGES
    )
    
    cat("ok\n")
    
    ## check if this is a fixed set (which is not altered), in this case we need to increase the nr. of max. edges
    ## FIXME: move to initEngine()
    if(!is.null(FIXED.INIT.EDGES)) {
      Grid.obj$smax =  Grid.obj$smax + sum(Grid.obj$edge.struct) - Grid.obj$additional.parents
      Grid.obj$FIXED.INIT.EDGES = which(Grid.obj$edge.struct[1:(length(Grid.obj$edge.struct) - Grid.obj$additional.parents )] == 1)
    }
  }
  
    
  
  ##
  ##  
  ## Execute the main MCMC procedure in main()
  ##
  ##
  cat("[", method.name, "] Starting MCMC chain.. \n")
  
  result = main(X, Y, start.iter, end.iter, MCMC.chain, Grid.obj, HYPERvar)
  
  
  ## Add the hyper variables for further reference
  ## Will be used when chain is continued in another run. 
  result$HYPERvar = HYPERvar
  
  ## Add the predictor data and target node data for further reference
  result$X = X
  result$Y = Y
  
  cat("\n---------------------------------------------------\n")
  cat("End of MCMC iterations")
  
  return(result)
  
}

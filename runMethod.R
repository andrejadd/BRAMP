
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



runMethod <- function(dataid=NULL, target=NULL, runid=NULL, niter=NULL, data.prefix=NULL, ENABLE.SAC=T) {  
  
  ## remove all but arguments  
  rm(list= ls()[ls()!="dataid" && ls()!= "target" && ls()!="runid" && ls()!="niter" && ls()!="data.prefix" && ls()!="ENABLE.SAC"])

  method.name = "BRAMPi"
  
  end.iter = niter

  ## the Mondrian process start budget
  start.budget = 1
  
  ## the start of iteration can be changed through loading an already existing result file
  start.iter = 1
 
  ## Unix time (seconds since 1970), should be fine for the cluster, might add/subtract milliseconds 
  set.seed(as.numeric(Sys.time()))

  ## the file we write to
  result.file = paste("./Results/Result_", data.prefix, "_id", dataid, "_n", target, "_run", runid, sep="")

  
  ## flag that tells if to proceed MCMC chain from previous run
  PROCEED.CHAIN = F

  ## flags result file existence
  result.exists.str = ""

  ## set these here, so check below will not fail (otherwise I would use exists() but even then I have to check for NULL)
  Grid.obj = NULL
  MCMC.chain = NULL
  X = NULL
  Y = NULL

  cat("[", method.name, "] data id", dataid, ", target: ", target, ", runid: ", runid, "\n")

  ## check if a previous result file exists
  if(file.exists(result.file)) {
    
    ## open to see if I can generate more iterations
    load(result.file)

    ## are all objects there to proceed chain?
    if(!is.null(Grid.obj) && !is.null(MCMC.chain) && !is.null(X) && !is.null(Y) && !is.null(HYPERvar)) {

      iters = MCMC.chain$Structsamples$iter

      last.iter = iters[[length(iters)]]

      ## check if less iters are in the chain (the added value just makes sure we don't start senseless for small iteration amount)
      if( (last.iter + 100) < end.iter) {

        ## set the iteration from which to continue
        start.iter = last.iter + 1
        
        ## flag that we proceed this file, skips all the init stuff below!
        PROCEED.CHAIN = T

        cat("[", method.name, "] Proceeding with chain from ", start.iter, " -> ", end.iter, "\n")
        
      } else {
        result.exists.str = paste("and chain has all iterations (", last.iter, " and required were ", end.iter, ").",sep="")
      } 
      
    } else {
      result.exists.str = "but chain data missing."
    }
    
  }


  ## if a result exists and the chain is ok or some data does not exist (old version), do nothing
  #if(RESULT.EXISTS && !PROCEED.CHAIN) {
  if(nchar(result.exists.str) > 0) {
    ## nothing to do
    cat("[", method.name, "] Skipping run: result file exists (",result.file,") ",result.exists.str, "\n")
    cat("[", method.name, "] Exit.\n")
    return(1)
      
  }
  
  ## look if we start a chain
  if(!PROCEED.CHAIN) {

    ## this will hold the actual values of the chain states, I declare it here because it is passed and might be filled
    ## when I load a already existing chain
    MCMC.chain = NULL

    ## This variable is used to tell the method not to alter some specific initially set edges during run
    ## It should be NULL if no fixed edges are assumed (except of course the bias and optional the SAC edge)
    ##    otherwise it is a vector with the parent node indices that define the fixed incoming edges
    
    ## check if this makes sense, seems confusing: FIXED.INIT.EDGES will be set to this value and INIT.EDGES.FROM.FILE will be set to a file from where to read these edges - which causes the initilization routine to read the fixed edges from a file.
    FIXED.INIT.EDGES=NULL

    # FIXED.INIT.EDGES=seq(1,12) ## use for the Outer Hebrides data when soil attributes (node 1..12) shall be fixed

    indata = paste("../Data/", data.prefix, "/Data_", data.prefix, "_id", dataid, ".Rdata",sep="")

    cat("[", method.name, "] Input: ", indata,"\n")
    
    ## the main sample data that is loaded is a target/predictor matrix [nodes x locations] and a SAC (spatial autocorrelation) matrix with
    ## dame size. The SAC matrix has for each target and location a SAC node which is included into the linear regression 
    load(indata)
    
    ## number of putative parents (sign q)
    nr.parents=dim(Model$Ymatrix)[1]-1
    
    ## Maximum number of parent nodes (fan-in restriction). A low limit is needed 
    ## for birth proposals based on precomupted posterior distribution (method 4),
    ## otherwise you can set smax = q
    smax = min(8,nr.parents);
    
    ## minimum number of locations in a segment
    minSegsize = 9
    
    cat("[", method.name, "] locs: " , Model$xlocs*Model$ylocs, " with x: ",  Model$xlocs, " , y: ",  Model$ylocs, ", start.budget: ", start.budget, ", parents: ", nr.parents, ", fan-in: ", smax, ", minSegsize: ", minSegsize, "\n")


    ## NOTE: fullData and sacData are already scaled for the hebrides data, doing it again via scale() does not bring changes
    ##       just keeping it for other maybe not scaled data
    ## Build response Y and predictor X
    ## extract the target values and scale
    Y = as.vector(scale(Model$Ymatrix[target,]))
    
    ## Helper vector of putative predictor nodes (default all) except the target itself
    posTF=c(1:(nr.parents + 1))[-c(target)]
    
    ## extract only predictors (excluding the target defined in posTF) and scale
    X = scale(t(Model$Ymatrix[posTF,]))
    
    ## add the bias node
    X = cbind(X,array(1,length(X[,1])))

    ## check if we want spatial autocorrelation (SAC Nodes)
    if(ENABLE.SAC) {

      ## check if SAC nodes already provided
      if(!is.null(Model$SAC.nodes)) {
        
        ## add a constant vector to predictor data (representing the bias nodes) and the spatial autocorrelation data of the target node
        X = cbind(X, as.vector(scale(Model$SAC.nodes[target,])))
        
      } else {
        ## otherwise, try to calculate the SAC nodes
        spatAC = spatAutoCorrelation(Y,Model$xlocs, Model$ylocs)
        
        ## scale the SAC node and append to X
        X = cbind(X, scale(spatAC))
      }
      
      additional.parents = 2  # bias and SAC parent
      
      cat("[", method.name, "] enabled spatial autocorrelation (SAC) node\n")
      
    } else {
      
      additional.parents = 1 # the bias parent

      cat("[", method.name, "] disabled spatial autocorrelation (SAC) node\n")

    }

    ## INIT hyper parameters
    ## for the signal-to-noise ratio, delta2 ~ IG(alphad2,betad2)
    alpha.snr = 2
    beta.snr = 0.2
    v0 = 1 # as in BRAM
    
    HYPERvar = list( c = 0.5, # for edge move proposals
      alpha.var = v0/2, beta.var = v0/2,   # for weight variance sigma2, v0/2 as in Marcos AISTATs paper
      alpha.snr=alpha.snr, beta.snr=beta.snr, delta.snr = (1 / rgamma(1, shape=alpha.snr, scale=1/beta.snr)),  # for signal-to-noise (snr) delta2 ~ IG(alphad2, betad2)
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
      xlocs=Model$xlocs,
      ylocs=Model$ylocs,
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
  
    
    
  ## run niter iterations
  cat("[", method.name, "] Starting MCMC chain.. \n")
  result = main(X, Y, start.iter, end.iter, MCMC.chain, Grid.obj, HYPERvar)
  cat("\n---------------------------------------------------\n")
  cat("End of iterations")
  
  cat("\n[END] attempting to write MCMC.chain and Grid.obj to: ", result.file, "\n")

  ## retrieve chain and last state
  MCMC.chain = result$MCMC.chain
  Grid.obj = result$Grid.obj
  
  save(MCMC.chain, Grid.obj, X, Y, HYPERvar, file = result.file)
  
}

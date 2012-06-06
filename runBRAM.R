
codePath=paste(getwd(),"/Code/",sep="")

source(paste(codePath,"moves.R",sep=""))
source(paste(codePath,"main.R",sep=""))
source(paste(codePath,"initEngine.R",sep=""))
source(paste(codePath,"segments.R", sep=""))
source(paste(codePath,"hyperParms.R", sep=""))
source(paste(codePath,"util.R",sep=""))   #requires pseudoinverse
source(paste(codePath,"sample_and_update.R",sep=""))  #requires pseudoinverse  
source(paste(codePath, "spatAutoCorrelation.R", sep=""))
source(paste(codePath,"invGamma.R",sep=""))
source(paste(codePath,"mvrnorm.R",sep=""))
source(paste(codePath,"ginv.R",sep=""))
source(paste(codePath,"convert.R",sep="")) ## need this ??
source(paste(codePath,"Tree.R",sep=""))



runBRAM <- function(dataid=NULL, target=NULL, runid=NULL, niter=NULL, start.budget=1, ENABLE.SAC=T){

  ## remove all but arguments  
  rm(list= ls()[ls()!="dataid" && ls()!= "target" && ls()!="runid" && ls()!="niter" && ls()!="start.budget" && ls()!="ENABLE.SAC"])

  end.iter = niter
  
  ## the start of iteration can be changed through loading an already existing result file
  start.iter = 1
 
  ## Unix time (seconds since 1970), should be fine for the cluster, might add/subtract milliseconds 
  set.seed(as.numeric(Sys.time()))

  ## the file we write to
  result.file = paste("./Results/SC2D_m", dataid, "_i", target, "_run", runid, sep="")

  ## flag that tells if to proceed MCMC chain from previous run
  PROCEED.CHAIN = F

  ## flags result file existence
  result.exists.str = ""

  ## set these here, so check below will not fail (otherwise I would use exists() but even then I have to check for NULL)
  Grid.obj = NULL
  MCMC.chain = NULL
  X = NULL
  Y = NULL

  cat("[runBRAM] data id", dataid, ", target: ", target, ", runid: ", runid, "\n")

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

        cat("[runBRAM] Proceeding with chain from ", start.iter, " -> ", end.iter, "\n")
        
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
    cat("[runBRAMflex] Skipping run: result file exists (",result.file,") ",result.exists.str, "\n")
    cat("[runBRAMflex] Exit.\n")
    return(1)
      
  }
  
  ## look if we start a chain
  if(!PROCEED.CHAIN) {

    ## this will hold the actual values of the chain states, I declare it here because it is passed and might be filled
    ## when I load a already existing chain
    MCMC.chain = NULL

    ## This variable is used to tell BRAM to not alter initially set edges during run (which are defined in a Data_id file, see initEdges..()
    ## It can be NULL (no fixed edges) or a vector of edge indicies (e.g., c(1) but it is just a placeholder and set later)
    ## FIXED.INIT.EDGES will be set to this value and INIT.EDGES.FROM.FILE will be set to a file from where to read these edges - which causes the initilization routine to read the fixed edges from a file.
                                        #FIXED.INIT.EDGES=c(1)

    FIXED.INIT.EDGES=seq(1,12) ## hebrides soil edges

    ## setup Data directory
    DATA.TYPE = "LOTKA.VOLTERRA"
    if(dataid < 100) DATA.TYPE = "SYNTHETIC"
    if(dataid < 10) DATA.TYPE = "HEBRIDES.PLANTS"
    
    indata = paste("../Data/", DATA.TYPE, "/Data_id", dataid, ".Rdata",sep="")

    cat("[runBRAMflex] Input: ", indata,"\n")
    
    
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
    
    cat("[runBRAMflex] locs: " , Model$xlocs*Model$ylocs, " with x: ",  Model$xlocs, " , y: ",  Model$ylocs, ", start.budget: ", start.budget, ", parents: ", nr.parents, ", fan-in: ", smax, ", minSegsize: ", minSegsize, "\n")


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
      
      cat("[runBRAMflex] with SAC node\n")
      
    } else {
      
      additional.parents = 1 # the bias parent
      
      cat("[runBRAMflex] no SAC node\n")
    }

    ## Create HyperParms Variables used in all functions
    HYPERvar = HyperParms(alphaCP = 1, betaCP = 0.5, alphaTF = 1, betaTF = 0.5)

    ## Note, this is just a placeholder, define the file below and look into initEngine to make use
    ## of predefined edges (e.g. from another method like Lasso)
    INIT.EDGES.FROM.FILE = NULL
    
    ## if this is not null it means we can set the file from where to read it
    #if(!is.null(FIXED.INIT.EDGES)) {
    #  INIT.EDGES.FROM.FILE = paste("../Data/", DATA.TYPE, "/Data_id", dataid, "_edgeprobs.Rdata", sep="")
    #}
    
    cat("[runBRAMflex] Init.. ")
  
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
  cat("[runBRAMflex] Starting BRAMflex chain.. \n")
  result = main(X, Y, start.iter, end.iter, MCMC.chain, Grid.obj, HYPERvar)
  cat("\n---------------------------------------------------\n")
  cat("End of iterations")
  
  cat("\n[END] attempting to write MCMC.chain and Grid.obj to: ", result.file, "\n")

  ## retrieve chain and last state
  MCMC.chain = result$MCMC.chain
  Grid.obj = result$Grid.obj
  
  save(MCMC.chain, Grid.obj, X, Y, HYPERvar, file = result.file)
  
}

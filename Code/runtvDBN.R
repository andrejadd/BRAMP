
codePath=paste(getwd(),"/Code/",sep="")
print(codePath)

## +  main functions
source(paste(codePath,"init.R", sep=""))
source(paste(codePath,"moves.R",sep=""))
source(paste(codePath,"main.R",sep=""))
source(paste(codePath,"output_functions.R", sep=""))

## +  useful tools
source(paste(codePath,"hyperParms.R", sep=""))
source(paste(codePath,"util.R",sep=""))#requires pseudoinverse
source(paste(codePath,"sample_and_update.R",sep=""))#requires pseudoinverse  
source(paste(codePath,"extractData.R",sep=""))
source(paste(codePath, "spatAutoCorrelation.R", sep=""))

## +  external functions
source(paste(codePath,"invGamma.R",sep=""))
source(paste(codePath,"mvrnorm.R",sep=""))
# source(paste(codePath,"fast.svd.R",sep="")) 
# source(paste(codePath,"pseudoinverse.R",sep=""))
source(paste(codePath,"ginv.R",sep=""))
source(paste(codePath,"simulate_network.R",sep=""))
source(paste(codePath,"convert.R",sep=""))
  



##modifie par Sophie 01/03/09: ajout de parametres
##modifie par Sophie 01/03/09: ajout de targetNamesFile
##### Main function to use to run the whole program
runtvDBN <- function(fullData,      ## all the data, response and predictors
                     sacData,       ## data of spatial autocorrelation nodes
                     q,             ## nr. of existing (putative) edges
                     nr.locations,  ## total nr. of locations on the grid
                     xlocs, ylocs,
                     minPhase=2,
                     smax,               ## max. number incoming edges (fan-in)
                     kmax.x, kmax.y,     ## max. number of CP along x and y axis  
                     alphaCP=1, betaCP=0.5, ## hyperparms for sampling the number k of CP : k ~ Gamma(alphaCP,betaCP)  (default: alphaCP=0.5, betaCP =1). You can use function choosePriors to set alphaCP and  betaCP according to the desired dimension penalisation.
                     alphaTF=1, betaTF=0.5, ## hyperparms for sampling the number l of TF : l ~ Gamma(alphaTF,betaTF) (default: alphaTF=0.5, betaTF =1).  You can use function choosePriors to set alphaTF and  betaTF according to the desired dimension penalisation.
                     niter=0,             ## number of iterations of the MCMC
                     modelid=NULL,        ## id of Data aka simulation model
                     runid=NULL,          ## id of run, can differ with same settings
                     target=NULL,         ## the actual target node for which the inference is done
                     FIXED.INIT.EDGES=NULL)  ## specifies if there are edges, which are not changed during run, see GLOBvar
{


  # fixme: m = number of samples per location (default 1) , should be nr. of locations along the y axis
  
  # Which method to use for proposing a new network structure in the CP birth 
  # move:
  #      1 - Old (incorrect) method
  #      2 - Sample from prior (correctly)
  #      3 - Sample based on Hamming Distance
  #      4 - Sample based on posterior distribution (calculated with no 
  #          changepoints)
  #      5 - Sample based on mixturee of prior and Hamming distance (2+3)
  birth_proposals = 2

  ## nr. of repeated measurements per location
  m = 1  
  dyn=0

  ## NEEDED? TRUE when a specific variance is estimated for each segment, FALSE otherwise (default: TRUE).
  multipleVar=TRUE
  
  ##  Mphase=seq(1,n*m+1,by=m)-dyn*m
  XMphase=seq(1,xlocs * m+1,by=m) - dyn * m
  YMphase=seq(1,ylocs * m+1,by=m) - dyn * m
 
  
  ## Create Global Variables used in all functions
  GLOBvar = list(n=nr.locations, xlocs=xlocs, ylocs=ylocs,m=m, p=1, q=q, smax=smax, kmax.x=kmax.x, kmax.y=kmax.y, dyn=dyn, 
    minPhase=minPhase, XMphase=XMphase, YMphase=YMphase, niter=niter, birth_proposals=birth_proposals,
    modelid=modelid, runid=runid, target=target, INIT.EDGES.FROM.FILE=NULL, FIXED.INIT.EDGES=FIXED.INIT.EDGES)

  ## if this is not null it means we can set the file from where to read it
  if(!is.null(FIXED.INIT.EDGES)) {
    GLOBvar$INIT.EDGES.FROM.FILE = paste("../Data/Data_id", modelid, "_edgeprobs.Rdata", sep="")
  }
  
  ## Create HyperParms Variables used in all functions
  HYPERvar = HyperParms(alphaCP, betaCP, alphaTF, betaTF)

  ## NOTE: fullData and sacData are already scaled for the hebrides data, doing it again via scale() does not bring changes
  ##       just keeping it for other maybe not scaled data
  ## Build response Y and predictor X
  ## extract the target values and scale
  Y = as.vector(scale(fullData[target,]))

  ## Helper vector of putative predictor nodes (default all) except the target itself
  posTF=c(1:(q+1))[-c(target)]

  ## extract only predictors (excluding the target defined in posTF) and scale
  X = scale(t(fullData[posTF,]))

  ## add a constant vector to predictor data (representing the bias nodes) and the spatial autocorrelation data of the target node
  X = cbind(X,array(1,length(X[,1])), as.vector(scale(sacData[target,])))

  ## NOTE: use this block if to calculate the SAC right here (based on the grid if there is no sacData vector)
  # spatAC = spatAutoCorrelation(Y,GLOBvar$xlocs,GLOBvar$ylocs)
  # scale the SAC node and add to X
  # X = cbind(X, scale(spatAC))

  print("Extract data ok")
    
  ## initialize system
  initiation = init(X, Y, GLOBvar, HYPERvar)
  print("Initialisation ok")

  ## check if this is a fixed set (which is not altered), in this case we need to increase the nr. of max. edges
  if(!is.null(GLOBvar$FIXED.INIT.EDGES)) {
    GLOBvar$smax =  GLOBvar$smax + sum(initiation$initState$S2Dall) - 2
    GLOBvar$FIXED.INIT.EDGES = which(initiation$initState$S2Dall[1:(length(initiation$initState$S2Dall) - 2)] == 1)
  }

  
  ## run niter iterations
  print("Starting tvDBN iterations...")
  runiteration = main(X, Y, initiation, GLOBvar, HYPERvar)
  print("---------------------------------------------------")
  print("End of iterations")
  
  
}

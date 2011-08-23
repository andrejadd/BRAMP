  # Import code :


# AA: hä?? warum ändert sich codepath wieder zu dem alten standard pfad
codePath=paste(getwd(),"/Code/",sep="")
print(codePath)

  # +  main functions
  source(paste(codePath,"buildXY.R",sep=""))  
  source(paste(codePath,"init.R",sep=""))
  source(paste(codePath,"moves.R",sep=""))
  source(paste(codePath,"main.R",sep=""))
  source(paste(codePath,"output_main.R",sep=""))
  source(paste(codePath,"output_functions.R", sep=""))
  # +  useful tools
  source(paste(codePath,"hyperParms.R", sep=""))
  source(paste(codePath,"util.R",sep=""))#requires pseudoinverse
  source(paste(codePath,"sample_and_update.R",sep=""))#requires pseudoinverse  
  source(paste(codePath,"extractData.R",sep=""))

## +  external functions
  source(paste(codePath,"invGamma.R",sep=""))
  source(paste(codePath,"mvrnorm.R",sep=""))
  source(paste(codePath,"fast.svd.R",sep=""))
  source(paste(codePath,"pseudoinverse.R",sep=""))
  source(paste(codePath,"simulate_network.R",sep=""))
  source(paste(codePath,"convert.R",sep=""))
  
  ## AA: added - my PSRF edge functions
  source(paste(codePath,"Helper/psrf.R",sep=""))
  source(paste(codePath, "spatAutoCorrelation.R", sep=""))

##modifie par Sophie 01/03/09: ajout de parametres
##modifie par Sophie 01/03/09: ajout de targetNamesFile
##### Main function to use to run the whole program
runtvDBN <- function(fullData, sacData, q, n, xlocs, ylocs, m,
                     multipleVar=TRUE, minPhase=2, smax, kmax, 
                     alphaCP=1, betaCP=0.5, alphaTF=1, betaTF=0.5,
                     pkCP=NULL, pkTF=NULL, bestPosMatFile=NULL, 
                     niter=0, predNamesFile=NULL, 
                     targetNamesFile=NULL,
                     kpriorsfile=paste(codePath,"k_priors.txt",sep=""), 
                     outputFile="tvDBNoutput",
                     modelid=NULL, runid=NULL, target=NULL){
  # runtvDBN: run whole program
  # Description:
  #
  # Arguments:
  # targetdata = target data (either the name of a file, or directly a matrix)
  # preddata = optional, file with predictor data when differing from the target data (either the name of a file, or directly a matrix), default=NULL.
  # q = number of existing parents
  # n = number of locations along the x axis
  # fixme: m = number of samples per location (default 1) , should be nr. of locations along the y axis

  ## ?? brauche ich? p = number of response variables (default=1)
                                        
  # multipleVar = TRUE when a specific variance is estimated for each phase, FALSE otherwise (default: TRUE).
  # minPhase = minimal length of a phase (>1 if no repetition, default: 2)
  # kmax = maximal number of CP 
  # smax = maximal number of TF 
  # alphaCP, betaCP = hyperparms for sampling the number k of CP : k ~ Gamma(alphaCP,betaCP)  (default: alphaCP=0.5, betaCP =1). You can use function choosePriors to set alphaCP and  betaCP according to the desired dimension penalisation.
  # alphaTF, betaTF= hyperparms for sampling the number l of TF : l ~ Gamma(alphaTF,betaTF) (default: alphaTF=0.5, betaTF =1).  You can use function choosePriors to set alphaTF and  betaTF according to the desired dimension penalisation.
  # pkCP, pkTF = prior distribution for the nulmber of CP or TF (default=NULL, necessary when BFOut = TRUE and the hyperparameters alpha, beta are note among the one provided by the package, see choosePriors for the list of available hyperparameters).
  # posResponse = row position (in targetdata) of targets to be analyzed (default: analize all genes)
  # bestPosMatFile = file containing row position of predictors for each gene (see default below)
  # niter = number of iterations (default: 20000)
  # predNamesFile = file containing the names of the predictor (matrix format) (by default rownames of preddata will be used)
  # outputFile = name of output file (default: tvDBNoutput)
  
  # fixed parameters
  # Position of each time point in the data (designed for the algorithm)
  # AA: a help vector that increments from Mphase[1] = 0 to Mphase[n+1]= n, where 'n' is the nr. of timepoints in the serie
  #     Mphase is used to map the indices from the data in R matrices to the data in the actual data (which uses 0 to start the first element)

  # Do you want a 'BF' Result analysis ?
  BFOut = FALSE

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

  ##  Mphase=seq(1,n*m+1,by=m)-dyn*m
  XMphase=seq(1,xlocs * m+1,by=m) - dyn * m
  YMphase=seq(1,ylocs * m+1,by=m) - dyn * m
  
  # nbVarMax = maximal number of variances (default: kmax, put it to 1 if the variance is the same for all phases)
  if(multipleVar){
	nbVarMax=kmax+1
  }else{ 
	nbVarMax=1
  }

  
  ## The predictor nodes for the target (default all) except the target itself
  ## posTF is an array with the predictor indices, simple from 1 to q-1 (with excluded target) 
  posTF=c(1:(q+1))[-c(target)]

  ## in case a predictor file is specified, try to read it and set posTF again 
  if(!is.null(bestPosMatFile)) {

    ## check if valid, read it and set posTF
    if(is.character(bestPosMatFile) & file.exists(bestPosMatFile)){

      cat("reading a bestPosMatFile..\n")
      bestPosMat = read.table(bestPosMatFile)
      posTF=as.matrix(bestPosMat)[target,1:q]
      
    } else {
      stop(paste("File", bestPosMatFile,"is incorrect or does not exist\n"))
    }
  }

  
  ### Create Global Variables used in all functions
  ### modifie par Sophie 02/07/09 :Ajout de PredNames et TargetNames
  GLOBvar = list(n=n, xlocs=xlocs, ylocs=ylocs,m=m, p=1, q=q, smax=smax, kmax=kmax, dyn=dyn, 
    minPhase=minPhase, nbVarMax=nbVarMax, XMphase=XMphase, YMphase=YMphase, posTF=posTF, 
    niter=niter, birth_proposals=birth_proposals, modelid=modelid, runid=runid, target=target)

  ## modifie par Sophie 01/03/09
  ### Create HyperParms Variables used in all functions
  HYPERvar = HyperParms(alphaCP, betaCP, alphaTF, betaTF, pkCP, pkTF, kpriorsfile, n, q, kmax, smax, dyn, BFOut)

  ## Build response Y and predictor X
  ## extract the target values and scale
  Y = as.vector(scale(fullData[target,]))

  ## extract only predictors (excluding the target defined in posTF) and scale
  X = scale(t(fullData[posTF,]))
  
  ## add a constant vector to predictor data (representing the bias nodes) and the spatial autocorrelation data of the target node
  X = cbind(X,array(1,length(X[,1])), as.vector(scale(sacData[target,])))

  ## NOTE: use this block if to calculate the SAC right here (based on the grid if there is no sacData vector)
  # spatAC = spatAutoCorrelation(Y,GLOBvar$xlocs,GLOBvar$ylocs)
  # scale the SAC node and add to X
  # X = cbind(X, scale(spatAC))
  browser()
  print("Extract data ok")
    
  ## initialize system
  initiation = init(X, Y, GLOBvar, HYPERvar)
  print("Initialisation ok")
  
  ## run niter iterations
  print("Starting tvDBN iterations...")
  runiteration = main(X, Y, initiation, GLOBvar, HYPERvar)
  print("---------------------------------------------------")
  print("End of iterations")
  
  
}

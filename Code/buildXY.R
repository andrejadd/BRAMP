
readDataTS <- function(data, posI, t0, tf, m, n){
  # data = matrix to read
  # pos = position of interest
  # m = # of repetitions
  # n = # of timepoints
  # input data is order by repetitions
 
  ### WARNING ###
  # when targetData is read: t0 = dyn and tf = n
  # when predData is read: t0 = 0 and tf = n-dyn

  # sort positions per time (tps1#1 tps1#2 .. tps1#M tps2#1 ... tpsN#M, (sapply(): apply fct over a list )
  posT = c(sapply(1:n,seq,m*n,by=n))
  
  # not all timepoints are considered if dyn>0
  # perdictor : from  T0 to TF-dyn
  # target : from T0+dyn to TF
  posT = posT[(m*t0+1):(m*tf)]

  # output matrix (positions sorted and truncated)
  Y=data[posI,posT]
  return(Y)
}



buildXY <- function(targetData, predData, GLOBvar){
  ### Build response Y and predictor X

  ### assignement of global variables used here ###
  n = GLOBvar$n
  m = GLOBvar$m
  q = GLOBvar$q
  dyn = GLOBvar$dyn
  target = GLOBvar$target
  posTF= GLOBvar$posTF
  ### end assignement ###

  ##
  ## This is SAC 1, first calculate SAC and then standardize (scale) all data
  ## 

  # read target data
  Y = as.array(readDataTS(data=targetData, posI=target, t0=dyn, tf=n, m=m, n=n))

  # build the predictor matrix X with only allowed predictors (parents) in it, so target itself is not included
  X = t(readDataTS(data=predData, posI=posTF, t0=0, tf=n-dyn, m=m, n=n))

  # scale the predictors, do this here because we dont want to scale the constant vector
  X = scale(X)

  ## add a constant vector to predictor data
  X = cbind(X,array(1,length(X[,1])))

  ## AA, 22.02.2011, add the spatial autocorrelation vector BEFORE scaling the target data (Y)
  spatAC = spatAutoCorrelation(Y,GLOBvar$xlocs,GLOBvar$ylocs)

  # scale the SAC node and add to X
  X = cbind(X, scale(spatAC))

  # finally scale the target data
  Y = as.vector(scale(Y))
  
  # return formatted data
  return(list(X=X,Y=Y))
}
  
 



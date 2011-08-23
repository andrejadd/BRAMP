

buildXY <- function(fullData, sacData, GLOBvar){
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

  ## extract the target values
  Y = as.array(fullData[target,])

  ## extract only predictors (excluding the target defined in posTF)
  X = t(fullData[posTF,])
  
  # scale the predictors, do this here because we dont want to scale the constant vector
  X = scale(X)

  ## add a constant vector to predictor data and the spatial autocorrelation Data (sacData)
  X = cbind(X,array(1,length(X[,1])), scale(sacData))

  ## AA: use this block if to calculate the SAC right here (based on the grid)
  # spatAC = spatAutoCorrelation(Y,GLOBvar$xlocs,GLOBvar$ylocs)
  # scale the SAC node and add to X
  # X = cbind(X, scale(spatAC))

  # finally scale the target data
  Y = as.vector(scale(Y))
  
  # return formatted data
  return(list(X=X,Y=Y))
}
  
 



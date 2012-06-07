##########################################################
###        sample something .. what else is needed ?? 
##########################################################



## Sample k from a truncated Poisson distribution, (eq 4.5)
sampleK <- function(mini, maxi, lambda, nb){
  # AA: mini: minimal nr. of endges (set to 0)
  #     maxi: max nr. edges (=smax)
  #     lambda: nr. of expected edges
  #     nb: 1
  # it basically says, sample the probability of each k in [mini..maxi] following a poission distribution with given lambda
  # this is the part with prob, then select one k in [mini..maxi] dependent on the probability
  if( mini == maxi) { print("Error with sampling from a truncated Poisson: mini = maxi") }
  out = sample( mini:maxi, nb, replace=TRUE, prob=lambda^(mini:maxi)/apply(matrix(mini:maxi, 1, maxi-mini+1), 2, factorial))
  return(out)
  
}



sampledelta2Global <- function(seg.set, X, Y, HYPERvar, DEBUGLVL2 =F) {

  xlocs = seg.set$xlocs
  ylocs = seg.set$ylocs
  
  scaleForHomogeneousStruct = HYPERvar$betad2

  segids = getSegmentIDs(seg.set)
  
  ## loop over segments  
  for(segid in segids) {

    ## get regr. coeff. for this segment
    B = seg.set$edge.weights[which( seg.set$edge.weights[,1] == segid),2:ncol(seg.set$edge.weights)]

    ## get the predictor and target data
    x = extractData(seg.set, X, segid)
    y = extractData(seg.set, Y, segid)
    
    scaleForHomogeneousStruct = scaleForHomogeneousStruct + B[which(seg.set$edge.struct == 1)] %*% t(x[,which(seg.set$edge.struct == 1)]) %*% x[,which(seg.set$edge.struct==1)] %*% B[which(seg.set$edge.struct == 1)] / ( 2 * seg.set$sigma2) 
    
  }
  
  nrsegs = length(segids)

  ## Updating hyperparameters
  ## AA: do not 'sum(S2Dall) - 1' because the bias also has to be considered
  ## the shape is the count for the dimension of the projection matrix (nr. of counts), it must match in the following way
  ##     sum(S2Dall) * nrsegs = sum(dim((x[,which(S2Dall==1)]))[,2])  , where dim()[,2] gives the column count which is nr. parents + bias
  delta2 = rinvgamma(1, shape= sum(seg.set$edge.struct) * nrsegs + HYPERvar$alphad2, scale=scaleForHomogeneousStruct)

  return(delta2)

}

    


## AA, change arguments, instead of extracting from Sall, Ball, X, Y - directly pass from moves.R since they are already used there
updateSigGlobal <- function(seg.set, X, Y, delta2, v0, gamma0){

  S = seg.set$edge.struct
  
  sumP = 0

  ## loop over each segment
  for(segid in getSegmentIDs(seg.set)) {

    ## get the predictor and target data
    x = extractData(seg.set, X, segid)
    y = extractData(seg.set, Y, segid)
    
    matPx = computeProjection(as.matrix(x[, which(S == 1)]), delta2)
    
    sumP = sumP + (t(y) %*% matPx %*% y)
  }

  total.locs = getNrElements(seg.set)
  
  ## Inverse Gamma
  out = rinvgamma(1, shape=(v0 + total.locs) / 2, scale= (gamma0 + sumP) / 2)
  return(out)
  
}




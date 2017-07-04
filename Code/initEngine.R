
initEngine <- function(X, HYPERvar, additional.parents, nr.parents, start.budget,
                       xlocs, ylocs, minSegsize, smax, target, fixed_edges){

  ### assignement of hyperparameters variables used here ###
  alphalbd = HYPERvar$alphalbd 
  betalbd = HYPERvar$betalbd

  ## create segment data type, has single segment in it
  Seg.set = createGrid(xlocs = xlocs, 
                       ylocs = ylocs,
                       minSegsize = minSegsize, 
                       start.budget = start.budget,
                       additional.parents = additional.parents,
                       nr.parents = nr.parents,
                       smax = smax,
                       fixed_edges = fixed_edges
    )

  ## set the valid locations in the segment map
  valid.vec = rep(1,  xlocs * ylocs)  ## prepare a vector
  valid.vec[which(is.nan(X[,1]))] = 0  ## rows in the data matrix X with NaN are invalid positions
  
  valid.map = matrix(valid.vec, ncol=ylocs,  nrow=xlocs)  ## create transpose matrix so that the vector fits
  Seg.set$valid.locs.map = t(valid.map)  ## transpose to original and assign

  total.nr.parents = nr.parents + additional.parents

  
  ## Sample model structure, only one structure because this is homogeneous
  structure = matrix(0, 1, total.nr.parents)


  ## Set the edges that are fixed to 1.
  structure[1,fixed_edges] = 1
   

  ## the bias and SAC edges are the last. Set them here to 1. Its possible the SAC edge is not present and only the bias edge is set.
  structure[, c((nr.parents+1):total.nr.parents)] = 1

  
  ## Sample the sigma variance - one global value.
  ##   ~ IG(alpha/2,beta/2) 
  ##   Note, this is not the real sigma prior because we lack the knowledge of regression coefficients at this point.
  ##
  sigma.var = 1 / rgamma(n=1, shape=HYPERvar$alpha.var/2, scale=1/(HYPERvar$beta.var/2))
  
  
  ## create the regr. coefficient matrix: c(seg_id, regr. coeff,.., bias, [SAC] )
  weights = matrix(0,nrow=0,ncol=(1 + total.nr.parents))
  
  
  ## for all the segments created on createSegments (above)
  for(seg.id in getSegmentIDs(Seg.set)) {
  
    ## get the predictor data
    X_tmp = extractData(Seg.set, X, seg.id)

    ## sample new regr. coeff. vector
    newB <- array(0, nr.parents + additional.parents)
  
    for(l in which(structure[1,] == 1)){
      newB[l] <- rnorm(1, mean=0, sd = sqrt(HYPERvar$delta.snr * sigma.var * t(X_tmp[,l]) %*% X_tmp[,l]))
    }

    ## create entry for segment in the regression coefficient matrix
    weights = rbind(weights, c(seg.id, newB))

  }

  ## Make final variable assignment:
  Seg.set$edge.weights = weights  
  Seg.set$edge.struct = structure[1,]
  Seg.set$sigma.var = sigma.var

  
  ## Init the mean on the weight means prior as zero.
  Seg.set$mu.on.prior = rep(0, total.nr.parents)

  
  ## Init covariance matrix for weight sampling to the diagonal identity matrix I.
  Seg.set$Cov.mat.on.prior = diag(total.nr.parents) 
  
  
  ## Set number of available predictor nodes including bias and SAC node.
  Seg.set$total.nr.parents = total.nr.parents
  
  
  ## Increase the edge fan-in by the number of fixed_edges.
  #if(!is.null(fixed_edges))
  #  Seg.set$smax = Seg.set$smax + length(fixed_edges)
  
  
  return(Seg.set)
}



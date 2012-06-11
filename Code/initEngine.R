
initEngine <- function(X, HYPERvar, additional.parents, nr.parents, start.budget,
                       xlocs, ylocs, minSegsize, smax, target,
                       INIT.EDGES.FROM.FILE,
                       FIXED.INIT.EDGES){

  ### assignement of hyperparameters variables used here ###
  alphalbd = HYPERvar$alphalbd 
  betalbd = HYPERvar$betalbd

  ## create segment data type, has single segment in it
  Seg.set = createGrid(xlocs=xlocs, ylocs=ylocs,
    minSegsize=minSegsize,
    start.budget = start.budget,
    additional.parents=additional.parents,
    nr.parents=nr.parents,
    smax=smax,
    target=target,
    FIXED.INIT.EDGES=FIXED.INIT.EDGES
    )

  ## set the valid locations in the segment map
  valid.vec = rep(1,  xlocs * ylocs)  ## prepare a vector
  valid.vec[which(is.nan(X[,1]))] = 0  ## rows in the data matrix X with NaN are invalid positions
  
  valid.map = matrix(valid.vec, ncol=ylocs,  nrow=xlocs)  ## create transpose matrix so that the vector fits
  Seg.set$valid.locs.map = t(valid.map)  ## transpose to original and assign

  total.nr.parents = nr.parents + additional.parents

  ## sample model structure, only one structure because this is homogeneous
  structure = matrix(0, 1, total.nr.parents)

  ## check if to sample from prior or set edges 
  if(is.null(FIXED.INIT.EDGES)) {
  
    ## sample mean nr. edges (Lambda) with rgamma and sampleK (is the same inverse gamma for CPs and edges) 
    sPred = sampleK(0, smax, rgamma(1, shape=alphalbd, rate = betalbd), 1)    # scale s= 1/rate => f(x)= 1/(s^a Gamma(a)) x^(a-1) e^-(x/s)

    ## if there are any edges..
    if(sPred>0){
      ## set random position in Structure S to edge (1)
      structure[1, sample(1:nr.parents, sPred, replace=FALSE)] = array(1, sPred) 
    }
  } else {

    ## set the edges that are fixed to 1
    structure[1,FIXED.INIT.EDGES] = 1
    
    ## If you want to initialize edges from a file ..maybe use this later
    #cat("[initEngine] initializing preset edges from file.. ", INIT.EDGES.FROM.FILE, "\n")

    ## at what probability we accept an edge (not too low because we want to have the best edges in but not too much)
    #prob.cutoff = 0.4
  
    ## attempt to read file, this read in a talbe
    ##  edge.pp.mat = read.table("../Data/Data_id2_edgeprobs.txt", header=F)

    ## load from Rdata, post.edge.pp
    #load(INIT.EDGES.FROM.FILE)
    #edge.pp.mat = post.edge.pp
  
    ## get the edges for the target and remove target self-edge
    #S = edge.pp.mat[target,][-target]

    ## create the matrix with values of '1' where an edge is above the prob. cut-off
    #S = matrix((S > prob.cutoff) * rep(1, length(S)), nrow=1, ncol=length(S))

    ## append the bias and SAC position 
    #S = cbind(S, c(0))
    #S = cbind(S, c(0))
    
  }

  ## the bias and SAC edges are the last. Set them here to 1. Its possible the SAC edge is not present and only the bias edge is set.
  structure[, c((nr.parents+1):total.nr.parents)] = 1

  ## FIXME: set the initial variance according to BRAMPi prior
  ## create the variance sigma - one global value  : IG(alpha/2,beta/2), note this is not the real sigma prior because we lack the knowledge of regre. at this point
  sigma.var = 1 / rgamma(n=1, shape=HYPERvar$alpha.var/2, scale=1/(HYPERvar$beta.var/2))
  
  ## create the regr. coefficient matrix: c(seg_id, regr. coeff,.., bias, [SAC] )
  weights = matrix(0,nrow=0,ncol=(1 + total.nr.parents))
  
  ## for all the segments created on createSegments (above)
  for(seg.id in getSegmentIDs(Seg.set)) {
  
    ## get the predictor data
    x = extractData(Seg.set, X, seg.id)

    ## sample new regr. coeff. vector
    newB <- array(0, nr.parents + additional.parents)
  
    for(l in which(structure[1,] == 1)){
      newB[l] <- rnorm(1, mean=0, sd=sqrt(HYPERvar$delta.snr * sigma.var * t(x[,l]) %*% x[,l]))
    }

    ## create entry for segment in the regression coefficient matrix
    weights = rbind(weights, c(seg.id, newB))

  }

  ## make final assignments
  Seg.set$edge.weights = weights  

  Seg.set$edge.struct = structure[1,]

  Seg.set$sigma.var = sigma.var

  ## init the mean on the weight means prior: mu_n ~ N(mu.on.prior, Cov.mat.on.prior)
  Seg.set$mu.on.prior = rep(0, total.nr.parents)

  ## init covariance matrix for weight sampling to small values
  Seg.set$Cov.mat.on.prior = diag(total.nr.parents) # matrix(runif(total.nr.parents, 0, 0.001), nrow=total.nr.parents, ncol=total.nr.parents)

  ## incoming species + bias + SAC (if present)
  Seg.set$total.nr.parents = total.nr.parents
  
  return(Seg.set)
}



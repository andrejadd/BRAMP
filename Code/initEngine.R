
initEngine <- function(X, HYPERvar, additional.parents, nr.parents, start.budget,
                       xlocs, ylocs, minSegsize, smax, target,
                       INIT.EDGES.FROM.FILE,
                       FIXED.INIT.EDGES){

  ### assignement of hyperparameters variables used here ###
  ##AA needed here ???
  alphaD = HYPERvar$alphaD
  betaD = HYPERvar$betaD
  alphalbd = HYPERvar$alphalbd
  betalbd = HYPERvar$betalbd
  v0 = HYPERvar$v0
  gamma0 = HYPERvar$gamma0
  alphad2 = HYPERvar$alphad2
  betad2 = HYPERvar$betad2
  delta2 = HYPERvar$delta2

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

  ## sample model structure, only one structure because this is homogeneous
  S = matrix(0, 1, nr.parents + additional.parents)

  ## check if to sample from prior or set edges 
  if(is.null(FIXED.INIT.EDGES)) {
  
    ## sample mean nr. edges (Lambda) with rgamma and sampleK (is the same inverse gamma for CPs and edges) 
    sPred = sampleK(0, smax, rgamma(1, shape=alphalbd, rate = betalbd), 1)    # scale s= 1/rate => f(x)= 1/(s^a Gamma(a)) x^(a-1) e^-(x/s)

    ## if there are any edges..
    if(sPred>0){
      ## set random position in Structure S to edge (1)
      S[1, sample(1:nr.parents, sPred, replace=FALSE)] = array(1, sPred) # structure du model (1 si pred in the model)
    }
  } else {

    ## set the edges that are fixed to 1
    S[1,FIXED.INIT.EDGES] = 1
    
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
  S[, c((nr.parents+1):(nr.parents + additional.parents))] = 1

  ## create the variance sigma - one global value  : IG(v0/2,gamma0/2), note this is not the real sigma prior because we lack the knowledge of regre. at this point
  Sig2_2Dall = rinvgamma(n=1, shape=v0/2, scale=gamma0/2)
  
  ## create the regr. coefficient matrix: c(seg_id, regr. coeff,.., bias, [SAC] )
  B = matrix(0,nrow=0,ncol=(1 + nr.parents + additional.parents))

  ## for all the segments created on createSegments (above)
  for(seg.id in getSegmentIDs(Seg.set)) {
  
    ## get the predictor data
    x = extractData(Seg.set, X, seg.id)

    ## sample new regr. coeff. vector
    newB <- array(0, nr.parents + additional.parents)
  
    for(l in which(S[1,] == 1)){
      newB[l] <- rnorm(1, mean=0, sd=sqrt(delta2 * Sig2_2Dall * t(x[,l]) %*% x[,l]))
    }

    ## create entry for segment in the regression coefficient matrix
    B = rbind(B, c(seg.id, newB))
  }

  ## make final assignments
  Seg.set$edge.weights = B  

  Seg.set$edge.struct = S[1,]

  Seg.set$sigma2 = Sig2_2Dall
  
  return(Seg.set)
}



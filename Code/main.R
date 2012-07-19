####
####
####
#### main function

main <- function(data, Y, start.iter, end.iter, MCMC.chain, Grid.obj, HYPERvar){

  DEBUGLVL = 0

  ## used extensively for updates and end of MCMC
  total.nr.parents = Grid.obj$total.nr.parents
  
  ## counting Segment moves (1..3) and edge update moves (4..7)
  ##  counters = list() ## this will hold the move counters from below
  cptMove = array(0,7)
  acceptMove = array(0,7)

  if(is.null(MCMC.chain)) {
    MCMC.chain = list(Structsamples = list(segment.map = list(), struct = list(), iter=list(), regression.coeff=list(), mondrian.tree=list()), counters=list(), delta.snr = c() )
  }
  
  ## everything important stored here and saved at the end to a file
  ## Structsamples = list(segment.map = list(), struct = list(), iter=list(), regression.coeff=list(), mondrian.tree=list())
  
  ## when to start saving the regression coefficients (takes up memory)
  start.save.regr = 1 #end.iter - floor(end.iter * 1/5)

  # do main iteration
  r = start.iter

  while(r < end.iter) {

    r = r + 1
    
    ## ----------->   FIXME!!!!
    ## the move probabilities in vector rho1 - calculate in a different fashion or even set fixed
    ##    The move probabilities below come from the BRAM which used segment proposals in the probabilities below
    ##    We dont need this and setting to fixed portions should be ok, e.g. (add segment, remove segment, shift cut, edge move) = (0.2, 0.2, 0.2, 0.4)
    ## TEST the difference on synthetic data! Could this lead to worse convergence ?

    ## get number of segments
    k = getNrSegments(Grid.obj)
    
    ## mean nr. of changepoints (lambda)    
    mean.nr.segments = rgamma(1, shape= k + 1, rate = 1 + 0.5)

    rho=array(1,4)
    cD = 0.2
    rho[1] = cD * min(1, mean.nr.segments / (k+1) )                                            ## segment cut depends on mean nr. of segments  
    if(k == 0) { rho[2] = rho[1] } else { rho[2] = rho[1]+ cD * min(1, k / mean.nr.segments) } ## segment merge
    if(k > 0) { rho[3] = rho[2]+(1-rho[2])/3 } else{ rho[3] = rho[2] }        ## cp shift, will be edge move if ignored below

    ## Sample u to choose one of the 4 moves : CP birth, CP death, CP shift, Update phases.
    u1 = runif(1, 0, 1)
    
    ## FIXME, see above
    #    rho[1] = 0.5
    #    rho[2] = 0.6
    #    rho[3] = 1
    
    ## Run 1 out of the 4 moves (depending on the value of u)
    if (u1 < rho[1]){
      ## Segment split (birth) move: return the new model if the move is accepted, the previous model otherwise.
      out = cut.segment(Grid.obj, data, Y, HYPERvar, DEBUGLVL)

    } else if(u1 < rho[2]){
     
      ## Segment merge 
      out = merge.segment(Grid.obj, data, Y, HYPERvar, DEBUGLVL)

    } else if(u1 < rho[3]){

      ## shift cut 
      out = shift.cut(Grid.obj, data, Y, HYPERvar, DEBUGLVL, counter=r)


    } else {
      ## Update phases: return the new model if the move is accepted, the previous model otherwise.
      out = segment.update(Grid.obj, data, Y, HYPERvar, DEBUGLVL)
    }
    
    ## apply changes made in the move functions
    Grid.obj = out$Grid.obj
    cptMove[out$move] = cptMove[out$move]+1
    acceptMove[out$move] = acceptMove[out$move]+out$accept

    ## If some move was accepted or the edge moves were selected, we need to update the edge weights
    ## FIXME: it would be wiser to always update all segment weights, this is to gain performance but check
    ##        the draw back of this!
   if(out$accept == 1 || out$move > 3) {

## --------------- UPDATE PARAMETERS ----------------------------------------------------------------------

    
     parents.vec = which(Grid.obj$edge.struct == 1)
     segidx.vec = getSegmentIDs(Grid.obj)

    
     ## the means for the weights, updated below and used for the information sharing
     mu.vec = Grid.obj$mu.on.prior[parents.vec]
    
     ## the covariance matrix Sigma_n, updated at the very bottom, here fixed for testing
     Cov.mat = Grid.obj$Cov.mat.on.prior[parents.vec, parents.vec]
        
     ## signal to noise ratio, updated below
     delta.snr = HYPERvar$delta.snr                 

    
##-------- variance (sigma.var) for edge weight udpate -----------------------------------------------------------------------

    
     ## sigma (variance) hyper parameters
     alpha.var = HYPERvar$alpha.var
     beta.var = HYPERvar$beta.var
       
    
     ## loop segment
     for(segidx in segidx.vec) {
        
       ## get nr. of observations (locations)
       n.obs = getNrElements(Grid.obj, segidx)
       
       ## get the predictor and target data
       x = extractData(Grid.obj, data, segidx)
       X = t(as.matrix(x[, parents.vec]))   # took transpose to make it conform with Marcos Code for testing (might change, see below transposes of X)
       
       ## what is the projection matrix exactly doing, need it?
       ## matPx = computeProjection(as.matrix(x[, which(seg.set$edge.struct == 1)]), delta.snr)
        
       y = extractData(Grid.obj, Y, segidx)

       ## thats a [locs,parents] x [parents,1] ; since mu.vec is a vector with length parents, %*% tranforms it to a [parents,1] matrix
       mu.tilde = t(X) %*% mu.vec

       inv.Sigma.tilde = diag(n.obs) - t(X) %*% ginv( ginv(delta.snr * Cov.mat) + X %*% t(X) ) %*% X
       
       mahalanobis.dist = t(y - mu.tilde) %*% inv.Sigma.tilde %*% (y - mu.tilde)
       
       alpha.var = alpha.var + (n.obs/2)
       
       beta.var = beta.var + (mahalanobis.dist/2)
       
     }
     
     ## IF RGAMMA FAILS ("NaNs produced") shape or scale are negative, catch this above
     inv.variance = rgamma(1,shape=alpha.var, scale=(1/beta.var))
     sigma.var = 1 / inv.variance
     
     ## save (but if calculated befor every weight sampling , might be not needed)
     Grid.obj$sigma.var = sigma.var
    

## ------- weights and signal-to-noise ratio (delta) update -------------------------------------------------------------------

     ## sigma (variance) hyper parameters
     beta.snr  = HYPERvar$beta.snr

     ## reset weight matrix, the +1 is for the segment idx
     Grid.obj$edge.weights = matrix(0,nrow=0,ncol=(total.nr.parents+1))

     ## update the weights for the segments that changed
     for(segidx in segidx.vec) {

       ## get nr. of observations (locations)
       n.obs = getNrElements(Grid.obj, segidx)
      
       ## get the predictor and target data
       x = extractData(Grid.obj, data, segidx)
       X = t(as.matrix(x[, parents.vec]))   # took transpose to make it conform with Marcos Code for testing (might change, see below transposes of X)
       y = extractData(Grid.obj, Y, segidx)


       ## mean and covariance for weight
       Sigma.inv.star = ginv(delta.snr * Cov.mat) + X %*% t(X)
       mu.star = ginv(Sigma.inv.star) %*% (ginv(delta.snr * Cov.mat) %*% mu.vec + X %*% y)

       ## sample weights for existing edges
       weights = mvrnorm(mu=mu.star, Sigma= sigma.var * ginv(Sigma.inv.star))
       weights = t(weights)

       ## assign to proper edge idx
       full.weights = rep(0, total.nr.parents)
       full.weights[parents.vec] = weights

       ## append with seg id
       Grid.obj$edge.weights = rbind(Grid.obj$edge.weights, c(segidx, full.weights))

       ## calculate scale (beta) for signal-to-noise delta, Note vector transpose different than in equation:   (w - mu) inv(Cov) t(w - mu)
       ##                                                   because (w - mu) gives a 1x3 matrix, whereas we would expect a 3x1 to appear
       beta.snr = drop(beta.snr + 0.5 * inv.variance * (weights - mu.vec) %*% ginv(Cov.mat) %*% t(weights - mu.vec) )
     }
    
     ## alpha + (nr.segs * nr. active parents)/2 
     alpha.snr = HYPERvar$alpha.snr + (length(segidx.vec) * length(parents.vec))/2
     
     ## sample signal-to-noise delta 
     inv.delta.snr  = rgamma(1, shape=alpha.snr, scale=(1/beta.snr));
     delta.snr = 1 / inv.delta.snr;

    ## and update
    HYPERvar$delta.snr = delta.snr;
               

## ------------ resample mean and covariance for the weights ----------------------------------------------------------
      

    nr.segs = length(segidx.vec)

    # weight mean over segments, exclude first column value which is the segment id
    mean.weights = apply(Grid.obj$edge.weights[,parents.vec + 1, drop=F], 2, mean)  # only means of existing edges

    c.factor = Grid.obj$sigma.var * HYPERvar$delta.snr
      
    MU.star  = nr.segs/(nr.segs + c.factor)  * mean.weights;
    COV.star = c.factor/(nr.segs + c.factor) * diag(length(parents.vec))

    ## sample and udpate prior edge means 
    Grid.obj$mu.on.prior = rep(0,total.nr.parents)   ## reset (what if we leave '0' edges untouched? to use for later iterations when edge is selected, maybe helps performance?)
    Grid.obj$mu.on.prior[parents.vec] = mvrnorm(1, MU.star,COV.star)      ## set parents

    Grid.obj$Cov.mat.on.prior = diag(total.nr.parents) ## simple Covariance, option to do more sophisticated?
    
    } # end parameter updates


    
    ## Save Data
    ## every 10th iteration (to save memory and MCMC thin-out)
    if((r %% 10) == 0) {

      # save data, this is written to the major outfile
      MCMC.chain$Structsamples$struct[[length(MCMC.chain$Structsamples$struct) + 1]] = Grid.obj$edge.struct
      MCMC.chain$Structsamples$iter[[length(MCMC.chain$Structsamples$iter) + 1]] = r
      MCMC.chain$delta.snr = c(MCMC.chain$delta.snr, HYPERvar$delta.snr)
      MCMC.chain$counters[[length(MCMC.chain$counters) + 1]] = list(cptMove=cptMove, acceptMove=acceptMove)
 
      
      ## only save in later stage, because needs memory
      if(r > start.save.regr) {

        MCMC.chain$Structsamples$regression.coeff[[length(MCMC.chain$Structsamples$regression.coeff) + 1]] = Grid.obj$edge.weights


	## check if segmentation exists
      	if(max(Grid.obj$segment.map) > 1) {
        			     ## save segmentation matrix only if there are more than one segment
           MCMC.chain$Structsamples$segment.map[[length(MCMC.chain$Structsamples$segment.map) + 1]] = Grid.obj$segment.map
           #Structsamples$mondrian.tree[[length(Structsamples$mondrian.tree) + 1]] = Grid.obj$mondrian.tree
        
        } else {
          ## only a single segment: put placeholder to save memory
          MCMC.chain$Structsamples$segment.map[[length(MCMC.chain$Structsamples$segment.map) + 1]] = 1
        }


      } else {
        MCMC.chain$Structsamples$regression.coeff[[length(MCMC.chain$Structsamples$regression.coeff) + 1]] = 1
        MCMC.chain$Structsamples$segment.map[[length(MCMC.chain$Structsamples$segment.map) + 1]] = NaN ## means nothing was recorded
      }

      
      
    }

    #  print status every X and save data to disk iteration 
    if((r %% 500) == 0) {
       cat("\nr: ", r, "\t")
       cat("memory (MByte): ", object.size(MCMC.chain)/1048600)
       
     }
        
  } # end iteration

 
 
  return(list(MCMC.chain=MCMC.chain, Grid.obj=Grid.obj))
}


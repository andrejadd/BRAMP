##
##
##
## The main function that runs the MCMC simulation.
##
##

mcmc_main <- function(y, X, MCMC.chain, Grid.obj, HYPERvar, start.iter, end.iter){

  
  DEBUGLVL = 0

  
  ## used extensively for updates and end of MCMC
  total.nr.parents = Grid.obj$total.nr.parents

  
  ## Counting the segment moves (1..3) and edge update moves (4..7)
  cptMove = array(0,7)
  acceptMove = array(0,7)
  
  
  ##
  ## Use as running counter of the current MCMC iteration.
  ##
  r = start.iter

 
  ##
  ##
  ## Main MCMC Loop. 
  ##
  ##
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
    if(k > 0) { rho[3] = rho[2]+(1-rho[2])/3 } else{ rho[3] = rho[2] }                         ## cp shift, will be edge move if ignored below

    ## Sample a uniform number 'u1' to choose one of the 4 moves : CP birth, CP death, CP shift, Update phases.
    u1 = runif(1, 0, 1)
    
    ## FIXME, see above
    #    rho[1] = 0.5
    #    rho[2] = 0.6
    #    rho[3] = 1
    
    ## Run 1 out of the 4 moves (depending on the value of u)
    if (u1 < rho[1]){
      
      ## Make the Mondrian split movement.
      out = cut.segment(Grid.obj, X, y, HYPERvar, DEBUGLVL)

    } else if(u1 < rho[2]){
     
      ## Make the Mondrian merge movement.
      out = merge.segment(Grid.obj, X, y, HYPERvar, DEBUGLVL)

    } else if(u1 < rho[3]){

      ## Make the Mondrian shift movement of a boundary (cut). 
      out = shift.cut(Grid.obj, X, y, HYPERvar, DEBUGLVL, counter=r)


    } else {
     
      ## Execute edge moves (add, delete, flip):
      out = edge_moves(Grid.obj, X, y, HYPERvar, DEBUGLVL)
    
      }
    
    ## Apply changes made in the move functions.
    Grid.obj = out$Grid.obj
    cptMove[out$move] = cptMove[out$move]+1
    acceptMove[out$move] = acceptMove[out$move]+out$accept

    ## If some move was accepted or the edge moves were selected, we need to update the edge weights
    ##
    ## NOTE: It might be wiser to always update the edge weights of each segment. This here is done
    ##       to gain some computing performance. Check if inference quality get much better if we do 
    ##       the update more often (I doubt it).
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
    
    ## Loop over each segment.
    for(segidx in segidx.vec) {
        
      ## Get number of observations for this segment.
      n.obs = getNrElements(Grid.obj, segidx)
       
      ## Get the design matrix for this segment.
      X_tmp = extractData(Grid.obj, X, segidx)
      X_tmp = t(as.matrix(X_tmp[, parents.vec]))  

      ## Get response data for this segment.     
      y_tmp = extractData(Grid.obj, y, segidx)

       ## That is matrix multiplication with: [locs,parents] x [parents,1].
       ## Because mu.vec is a vector with length parents, %*% tranforms it to a [parents,1] matrix.
       mu.tilde = t(X_tmp) %*% mu.vec

       inv.Sigma.tilde = diag(n.obs) - t(X_tmp) %*% ginv( ginv(delta.snr * Cov.mat) + X_tmp %*% t(X_tmp) ) %*% X_tmp
       
       mahalanobis.dist = t(y_tmp - mu.tilde) %*% inv.Sigma.tilde %*% (y_tmp - mu.tilde)
       
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
      
       
       ## Get the predictor and response data for the segment.
       X_tmp = extractData(Grid.obj, X, segidx)
       X_tmp = t(as.matrix(X_tmp[, parents.vec]))  
       y_tmp = extractData(Grid.obj, y, segidx)


       ## mean and covariance for weight
       Sigma.inv.star = ginv(delta.snr * Cov.mat) + X_tmp %*% t(X_tmp)
       mu.star = ginv(Sigma.inv.star) %*% (ginv(delta.snr * Cov.mat) %*% mu.vec + X_tmp %*% y_tmp)

       ## sample weights for existing edges
       weights = mvrnorm(mu = mu.star, Sigma = sigma.var * ginv(Sigma.inv.star))
       weights = t(weights)

       ## assign to proper edge idx
       full.weights = rep(0, total.nr.parents)
       full.weights[parents.vec] = weights

       ## append with seg id
       Grid.obj$edge.weights = rbind(Grid.obj$edge.weights, c(segidx, full.weights))

       ## Calculate the scale (beta) parameter for the signal-to-noise delta parameter.
       ##  Note: The vector transpose is different than in equation: (w - mu) inv(Cov) t(w - mu).
       ##        This is because (w - mu) gives a 1x3 matrix, whereas we would expect a 3x1 matrix.
       beta.snr = drop(beta.snr + 0.5 * inv.variance * (weights - mu.vec) %*% ginv(Cov.mat) %*% t(weights - mu.vec) )
     }
    
     ## Update the alpha hyper-parameter for the SNR sample: 
     ##     alpha = alpha + ((nr.segs * nr. active parents)/2) 
     alpha.snr = HYPERvar$alpha.snr + (length(segidx.vec) * length(parents.vec))/2
     ## FIXME: save alpha.snr back to HYPERvar, is this necessary ?

     
     ## Sample the signal-to-noise delta parameter.
     inv.delta.snr  = rgamma(1, shape=alpha.snr, scale=(1/beta.snr));
     delta.snr = 1 / inv.delta.snr;

    ## and update
    HYPERvar$delta.snr = delta.snr;
               

## ------------ Resample mean and covariance for the regression coefficients (weights) ----------------------------------------------------------
      

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


    
    ## Save the MCMC chain with a thin-out of 10, i.e. only every 10th iteration to save memory.
    if((r %% (MCMC.chain$chain_thinout)) == 0) {

      ## Save the edge structure and the number of the iteration
      MCMC.chain$Structsamples$struct[[length(MCMC.chain$Structsamples$struct) + 1]] = Grid.obj$edge.struct
      MCMC.chain$Structsamples$iter[[length(MCMC.chain$Structsamples$iter) + 1]] = r
      
      ## Save some of the parameters - important for continuing a run of the chain.
      MCMC.chain$params = rbind(MCMC.chain$params, c(HYPERvar$alpha.snr, HYPERvar$beta.snr, HYPERvar$delta.snr, HYPERvar$alpha.var, HYPERvar$beta.var, Grid.obj$sigma.var))
      MCMC.chain$delta.snr = c(MCMC.chain$delta.snr, HYPERvar$delta.snr)
      
      ## Save the counters record the number of accepted moves - can be used for later analysis.
      MCMC.chain$counters[[length(MCMC.chain$counters) + 1]] = list(cptMove=cptMove, acceptMove=acceptMove)
 
      ##
      ## Only save in later stage, because needs memory
      ##

      if(r > MCMC.chain$iteration_save_betas) {

        ## Save the edge weights, i.e. regression coefficients.
        MCMC.chain$betas[[length(MCMC.chain$betas) + 1]] = Grid.obj$edge.weights


	      ## Check if segmentation exists.
      	if(max(Grid.obj$segment.map) > 1) {
      	  
          ## Save segmentation matrix only if there are more than one segment
          MCMC.chain$segment_map[[length(MCMC.chain$segment_map) + 1]] = Grid.obj$segment.map
          
        } else {
          
          ## This means that there is only a single segment, use a place holder '1', instead of saving the whole
          ## segment matrix
          MCMC.chain$segment_map[[length(MCMC.chain$segment_map) + 1]] = 1
        }


      } else {
        
        ## Do not save any valueable information, instead put a placeholder into the list for the regression coefficients
        ##  and the segmentation map. 
        ## The reason for this is to use the index of these lists as a lookup for the iteration in the chain.
        MCMC.chain$betas[[length(MCMC.chain$betas) + 1]] = 1
        MCMC.chain$segment_map[[length(MCMC.chain$segment_map) + 1]] = NaN ## means nothing was recorded
      
      }

    }

    
    ##  Prints the memory usage of the MCMC.chain data structure. 
    if((r %% 500) == 0) {
       cat("\nIteration: ", r, "\t")
       cat("Memory used (MByte): ", object.size(MCMC.chain)/1048600)
       
     }
        
  } ## End of main MCMC iteration loop


  return(list(MCMC.chain = MCMC.chain, Grid.obj = Grid.obj, HYPERvar = HYPERvar, y = y, X = X))

}


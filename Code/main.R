####
####
####
#### main function

main <- function(X, Y, start.iter, end.iter, MCMC.chain, Grid.obj, HYPERvar){

  DEBUGLVL = 0
  
  ## counting Segment moves (1..3) and edge update moves (4..7)
  ##  counters = list() ## this will hold the move counters from below
  cptMove = array(0,7)
  acceptMove = array(0,7)

  if(is.null(MCMC.chain)) {
    MCMC.chain = list(Structsamples = list(segment.map = list(), struct = list(), iter=list(), regression.coeff=list(), mondrian.tree=list()), counters=list() )
  }
  
  ## everything important stored here and saved at the end to a file
  ## Structsamples = list(segment.map = list(), struct = list(), iter=list(), regression.coeff=list(), mondrian.tree=list())
  
  ## when to start saving the regression coefficients (takes up memory)
  start.save.regr = end.iter - floor(end.iter * 1/5)

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
      out = cut.segment(Grid.obj, X, Y, HYPERvar, DEBUGLVL)

    } else if(u1 < rho[2]){
     
      ## Segment merge 
      out = merge.segment(Grid.obj, X, Y, HYPERvar, DEBUGLVL)

    } else if(u1 < rho[3]){

      ## shift cut 
      out = shift.cut(Grid.obj, X, Y, HYPERvar, DEBUGLVL, counter=r)


    } else {
      ## Update phases: return the new model if the move is accepted, the previous model otherwise.
      out = segment.update(Grid.obj, X, Y, HYPERvar, DEBUGLVL)
    }
    
    ## apply changes made in the move functions
    Grid.obj = out$Grid.obj
    cptMove[out$move] = cptMove[out$move]+1
    acceptMove[out$move] = acceptMove[out$move]+out$accept

    ## If some move was accepted or the edge moves were selected, we need to update the edge weights
    ## FIXME: it would be wiser to always update all segment weights, this is to gain performance but check
    ##        the draw back of this!
    if(out$accept == 1 || out$move > 3) {

#      cat("move: ", out$move)
#      if(!is.null(out$changed.segids)) { cat(", seg.ids: "); print.table(out$changed.segids) }
#      browser()

     
      ## upate variance, globally (over all segments) - and before calculating the regression coefficient
      Grid.obj$sigma2 = updateSigGlobal(Grid.obj, X,Y, HYPERvar$delta2, HYPERvar$v0, HYPERvar$gamma0)

      ## check if we should update the weights of all segments - this is the case after a segment update move (move > 3)
      if(out$move > 3) {
        out$changed.segids = getSegmentIDs(Grid.obj)
      }
      
      ## update the weights for the segments that changed
      for(seg.id in out$changed.segids) {
      
        ## create regression coefficient for new segment 
        x = extractData(Grid.obj, X, seg.id)
        y = extractData(Grid.obj, Y, seg.id)
      
        Pr = computeProjection(as.matrix(x[,which(Grid.obj$edge.struct == 1)]), HYPERvar$delta2)
      
        ## sample weight vector ( + 2 means Bias and SAC edge)
        weight.vec = array(0, Grid.obj$nr.parents + Grid.obj$additional.parents)      

        xi = x[,which(Grid.obj$edge.struct==1)]

        Ml = ( HYPERvar$delta2 / ( HYPERvar$delta2 + 1)) * ginv( t(xi) %*% xi)
      
        weight.vec[which(Grid.obj$edge.struct == 1)] = mvrnorm(1, mu=Ml %*% t(xi) %*% y, Sigma = (Grid.obj$sigma2 * Ml))

        
        ## check if weights for this seg id already exists:
        ## extract row in the Regr.Coeff. matrix for this seg.id
        seg.row = which( Grid.obj$edge.weights[,1] == seg.id)

        if(length(seg.row) == 1) {   # if the entry exists

          ## replace the old entry 
          Grid.obj$edge.weights[seg.row,] = c(seg.id, weight.vec)

        } else {

          ## append to matrix
          Grid.obj$edge.weights = rbind(Grid.obj$edge.weights, c(seg.id, weight.vec))
        }

       

        
      }

      ## Update of delta2 must come after update of edge weights because it needs the edge weights of all segment; in the case of a cut however
      ## the weights vector of the new segment is only created above. 
      HYPERvar$delta2 = sampledelta2Global(Grid.obj, X, Y, HYPERvar, F)
      
      
    }
   

        
    ## Save Data
    ## every 10th iteration (to save memory and MCMC thin-out)
    if((r %% 10) == 0) {

      # save data, this is written to the major outfile
      MCMC.chain$Structsamples$struct[[length(MCMC.chain$Structsamples$struct) + 1]] = Grid.obj$edge.struct
      MCMC.chain$Structsamples$iter[[length(MCMC.chain$Structsamples$iter) + 1]] = r

      
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
      
      MCMC.chain$counters[[length(MCMC.chain$counters) + 1]] = list(cptMove=cptMove, acceptMove=acceptMove)
 
    }

    #  print status every X and save data to disk iteration 
    if((r %% 500) == 0) {
       cat("\nr: ", r, "\t")
       cat("memory (MByte): ", object.size(MCMC.chain)/1048600)
       
     }
        
  } # end iteration

 
 
  return(list(MCMC.chain=MCMC.chain, Grid.obj=Grid.obj))
}


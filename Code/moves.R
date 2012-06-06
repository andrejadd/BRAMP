

#####################################################################################
## Cut a Segment
## 
#####################################################################################

cut.segment <- function(Grid.obj, X, Y, HYPERvar, DEBUGLVL = 0) {

  if(DEBUGLVL == 1) cat("START segment.split > ")
  
  ## get segment ids, all can be split which are in the segment map
  leaf.id = getSegmentIDs(Grid.obj)

  ## if there are many segment, chose one by uniform chance
  if(length(leaf.id) > 1) {
    leaf.id = sample(c(leaf.id, leaf.id), 1)
  }

  if(DEBUGLVL == 2) cat("[cut.segment] proposed leaf: ", leaf.id)

  
  ## get length and width
  parent.dim = getSegmentDimScale(Grid.obj, leaf.id)

  ## the cost of cutting this segment conditional on the half parameter
  cost = expdist(sum(parent.dim))     

  ## get the budget for this segment
  node = getLeafNode(Grid.obj$mondrian.tree, leaf.id)

  if(is.null(node)) {
   
    rnd.id = ceiling(runif(1,1,100000))
    save(file=paste("debug.out.", rnd.id, sep=""), Grid.obj, leaf.id) 
   
    stop("error getting the leaf node for id ", leaf.id, " ; node returned is NULL\n")
  }

  parent.budget = node$budget
 
  ## check if a cut can be done at all
  if(cost > parent.budget)  {
    if(DEBUGLVL == 2) cat("  -> cost too high\n")
    return(list(Grid.obj=Grid.obj, accept=0, move=1, alpha=0))
  }
  
  ## decide along which axis to cut
  axis="y"     ## cut along y

  if(runif(1) < (parent.dim[1]/(sum(parent.dim))))  axis="x"  # x axis, parent.dim[1] corresponds to x-length, as higher as more likely we cut along x

  ## get cut position 
  positions = getSplitPositions(Grid.obj, leaf.id, axis)

  if(length(positions) == 0) {
    ## no position found, exit cut
    if(DEBUGLVL == 2) cat("  -> no position found\n")
    return(list(Grid.obj=Grid.obj, accept=0, move=1, alpha=0))
  }
    
  ## extract random position
  position = sample(c(positions,positions), 1)

  ## the budget of the child segments
  child.budget = parent.budget - cost

  if(is.infinite(child.budget)) {
    stop("INFINITE BUDGET CAUGHT in cut\n")
  }

  ## do the cut (and save the budgets)
  proposed.set = splitSegment(data.obj=Grid.obj, parent.id=leaf.id, axis=axis, position=position)

  ## this is set inside splitSegment()
  new.id = proposed.set$last.added.id
  
  ## check if elements are large enough, if not return
  if( (getNrElements(proposed.set, leaf.id) < proposed.set$minSegsize) || (getNrElements(proposed.set, new.id) < proposed.set$minSegsize)) {
    if(DEBUGLVL == 2) cat("  ->  minSegSize exceeded\n")
    return(list(Grid.obj=Grid.obj, accept=0, move=1, alpha=0))
  }

  ## compute the Log Likelihoods for the new childs and parent
  log.alpha = cp.computeAlpha(HYPERvar, X, Y, Grid.obj, proposed.set, leaf.id, c(leaf.id, new.id), parent.budget, child.budget, DEBUGLVL) 

  ## transform from log to normal form
  alpha = exp(log.alpha)

  ## Sample u to conclude either to  acceptation or to rejection
  u = runif(1,0,1)
  
  ## Boolean for the acceptation of the CP birth move initially set to 0 (=1 if birth accepted, 0 otherwise)
  accept = 0

  if(!is.nan(alpha) & u <= alpha){

    if(DEBUGLVL == 2) cat("  -> accepted\n")
    
    if(DEBUGLVL == 3) cat("C")
    
    ## Acceptation of the birth of the new CP (boolean = 1)
    accept=1
    
    ## upate variance, globally (over all segments) - and before calculating the regression coefficient
    proposed.set$sigma2 = updateSigGlobal(proposed.set, X,Y, HYPERvar$delta2, HYPERvar$v0, HYPERvar$gamma0)
    
    ## create regression coefficient for new segment 
    x = extractData(proposed.set, X, new.id)
    y = extractData(proposed.set, Y, new.id)
    
    Pr = computeProjection(as.matrix(x[,which(Grid.obj$edge.struct == 1)]), HYPERvar$delta2)
    
    ## create regression coefficient vector ( + 2 means Bias and SAC edge)
    newB = array(0, Grid.obj$nr.parents + Grid.obj$additional.parents)      
    newB[which(Grid.obj$edge.struct == 1)] = sampleBxy(x[,which(Grid.obj$edge.struct==1)], y, proposed.set$sigma2, HYPERvar$delta2)

    ## append to edge weights matrix
    proposed.set$edge.weights = rbind(proposed.set$edge.weights, c(new.id, newB))
    
    Grid.obj = proposed.set

    ## do this here because alterations to the tree affect the pointer, if I would have made this on proposed.set before it would
    ## also have affected the tree inside Grid.obj - would require to copy the tree data not the pointer
    ## This way it is easier
    addLeafpair(Grid.obj$mondrian.tree, leaf.id, leaf.id, new.id, child.budget)

    ##
    ## start DEBUGGING
    ##

    
    if(DEBUGLVL == 2) {

      heatmap(Grid.obj$segment.map,Colv=NA,Rowv=NA,scale="none", col=grey.colors(20,start=1,end=0))
      drawTree(Grid.obj$mondrian.tree)

#      browser()
      
      z <- try(silent=TRUE, timeout(readline(prompt="Hit me: "), seconds=50))
      z = ""
      if(z == "q") { stop("USER EXIT") }
      if(z == "b") { browser() }
      cat("\n")
    }

    ##
    ## END DEBUGGING
    ##
  } else {
    if(DEBUGLVL == 2) cat("  -> not accepted\n")
  }

#  if(DEBUGLVL == 2) cat("  nr. segments: ", length(getSegmentIDs(Grid.obj)), "\n")
  
  ##  Return all variables
  ## (+ variable move describing the move type  (1= CP birth, 2= CP death, 3= CP shift, 4= Update phases)
  return(list(Grid.obj=Grid.obj, accept=accept, move=1, alpha=alpha))
  
}



#####################################################################################
# Segment merge
#####################################################################################

merge.segment <- function(Grid.obj, X, Y, HYPERvar, DEBUGLVL = 0) {

  if(DEBUGLVL == 1) { cat("START segment.merge > \n") }

  ## get the blocks that can be merged (leafs)
  leaf.pair = getRandomLeafPair(Grid.obj)
  
  ## check if leafs exist at all
  if(is.null(leaf.pair)) {
    if(DEBUGLVL == 2) cat("[merge.segment]  nothing to merge.\n")
    return( list(Grid.obj=Grid.obj, accept=0, move=2, alpha=0, estar=-1) )
  }

  ## identify the involved segments
  parent.id = leaf.pair$leaf1$parent$id
  child1.id = leaf.pair$leaf1$id
  child2.id = leaf.pair$leaf2$id

  if(DEBUGLVL == 2) cat("[merge.segment] proposed leaf pair: ", child1.id , " - ", child2.id)

  ## merges the segment in $segment.map, keeps Mondrian tree unaltered
  proposed.set = mergeSegment(Grid.obj, c(child1.id, child2.id))
  
  ## get the budget of the parent
  parent.budget = leaf.pair$leaf1$parent$budget

  ## get the child budget
  child.budget = leaf.pair$leaf1$budget

  ## compute the acceptance probabilityx
  log.alpha = cp.computeAlpha(HYPERvar, X, Y, proposed.set, Grid.obj, parent.id, c(child1.id, child2.id), parent.budget, child.budget ,DEBUGLVL) 

  ## invert to reverse the cut (-> a merge) and get rid of the log
  alpha = exp(-1 * log.alpha)
  
  ## Sample u to conclude either to  acceptation or to rejection
  u = runif(1,0,1)

  ## Boolean for the acceptation of the CP death move initially set to 0 (=1 if birth accepted, 0 otherwise)
  accept = 0

  ## check if to accept
  if(!is.nan(alpha) & u <= alpha){

    if(DEBUGLVL == 3) cat("M")

    ## Acceptation of the death of the selected CP, Move acceptation boolean =1
    accept=1

    ## delete the parameters for the deleted segment id
    ## FIXME: I could bind edge weights to a Mondrian tree node.. 
    proposed.set$edge.weights = proposed.set$edge.weights[which(proposed.set$edge.weights[,1] != max(child1.id, child2.id)),,drop=F]     # drop=F makes sure the matrix is not transformed to vector when single row is left

    ## !! FIXME: recalc the edge weights after merge!!!

    ## overwrite old state segment set
    Grid.obj = proposed.set

    ## delete the two childs

    if(!delLeaf(Grid.obj$mondrian.tree, child1.id)) stop("failed to delete leaf")
    if(!delLeaf(Grid.obj$mondrian.tree, child2.id)) stop("failed to delete leaf")
    
    # no need to update the variance because only needed for calculating the regression coefficient B
    # Sig2_2Dall = updateSigGlobal(XE, YE, X,Y, Grid.obj$edge.struct, delta2, HYPERvar$v0, HYPERvar$gamma0)
    
    if(DEBUGLVL == 2) {

      cat(" -> accepted\n")
      heatmap(Grid.obj$segment.map,Colv=NA,Rowv=NA,scale="none", col=grey.colors(20,start=1,end=0))
      drawTree(Grid.obj$mondrian.tree)

#      browser()
      z <- try(silent=TRUE, timeout(readline(prompt="Hit me: "), seconds=50))
      z=""
      if(z == "q") { stop("USER EXIT") }
      if(z == "b") { browser() }
    }

    ##
    ## END DEBUGGING
    ##
  } else {
    if(DEBUGLVL == 2) cat(" -> not accepted\n")
  }
    

#  if(DEBUGLVL == 2) cat("  nr. segments: ", length(getSegmentIDs(Grid.obj)), "\n")
  
  ##  Return all variables
  ## (+ variable move describing the move type  (1= Segment split, 2= Segm. merge, 3= Segment grow, 4= Update phases)
  return(list(Grid.obj=Grid.obj, accept=accept, move=2, alpha=alpha))

}


#################################################################################################
# Shift an existing cut in a leaf pair
#
# Description: this moves just repositions the cut, i.e. x or y and locations in the same manner
#               we would set a new cut.
#
# FIXME: Maybe we could create two moves: one that slightly shifts a cut along his orientation
#                                         another that changes the orientation: x -> y, y -> x        
#
# In all cases, the prior should not change because halfperimeter and cost of cut stay the same
#################################################################################################
shift.cut <- function(Grid.obj, X, Y, HYPERvar, DEBUGLVL = 0, counter=0) {

  if(DEBUGLVL == 1) { cat("START segment.merge > \n") }

  ## get the blocks that can be merge (leafs)
  leaf.pair = getRandomLeafPair(Grid.obj)

  ## check if leafs exist at all
  if(is.null(leaf.pair)) {
    if(DEBUGLVL == 2) cat("  nothing to merge, returning.\n")
    return( list(Grid.obj=Grid.obj, accept=0, move=3, alpha=0, estar=-1) )
  }

  ## identify the involved segments
  parent.id = leaf.pair$leaf1$parent$id
  child1.id = leaf.pair$leaf1$id
  child2.id = leaf.pair$leaf2$id

  ## prepare new segment, first delete the cut
  proposed.set = Grid.obj

  ## get the coordinates of both segments
  x1.ij = which(proposed.set$segment.map == child1.id, arr.ind=T)
  x2.ij = which(proposed.set$segment.map == child2.id, arr.ind=T)
  
  ## merge: assign the parent id (smallest) 
  proposed.set$segment.map[x1.ij] = parent.id
  proposed.set$segment.map[x2.ij] = parent.id

  ## get the halfperimeter 
  full.dim = getSegmentDimScale(proposed.set, parent.id)
  
  ## decide along which axis to cut
  axis="y"     ## cut along y

  if(runif(1) < (full.dim[1] / sum(full.dim) ))  axis="x"  # x axis, parent.dim[1] corresponds to x-length, as higher as more likely we cut along x

  ## get cut position 
  positions = getSplitPositions(proposed.set, parent.id, axis)

  if(length(positions) == 0) {
    ## no position found, exit cut
    return(list(Grid.obj=Grid.obj, accept=0, move=3, alpha=0))
  }

  position.mismatch = F
  
  ## FIXME: is there a fast way to exclude the position that already exists
  for(i in c(1:3)) {
    ## try to find a valid position
    
    ## extract random position
    position = sample(c(positions,positions), 1)

    ## do the cut 
    tmp.set = splitSegment(data.obj=proposed.set, parent.id=parent.id, new.id=max(child1.id, child2.id), axis=axis, position=position)
    
    no.match = sum(Grid.obj$segment.map != tmp.set$segment.map)

  ## check if elements are large enough, if not return
    if( (getNrElements(tmp.set, child1.id) < tmp.set$minSegsize) || (getNrElements(tmp.set, child2.id) < tmp.set$minSegsize)) {
      if(DEBUGLVL == 2) cat("  ->  minSegSize exceeded\n")
      no.match = 0
      next ## leads to exit from function when loop exhausted
    }
    
    ## see if something changed
    if(no.match > 0) {
      #cat("change! ")
      proposed.set = tmp.set 
      position.mismatch = T
      break
    }
  }

  if(DEBUGLVL == 2) cat("[shift.cut] poposed leaf pair: ", child1.id, " - ", child2.id, " to ", position, ", axis: ", axis, ", no.match: ", no.match)
 
  
  ## check if a valid cut position was found
  if(!position.mismatch) {
    if(DEBUGLVL == 2) cat("  no change from cut\n")
    return(list(Grid.obj=Grid.obj, accept=0, move=3, alpha=0))
  }

  ## calculate the probability of shift
  alpha = shift.computeAlpha(X, Y, Grid.obj, proposed.set, seg.ids=c(child1.id, child2.id), HYPERvar, DEBUG_BIRTH_EXT=0) 

  ## decide if to do it
  u = runif(1,0,1)
  
  ## Boolean for the acceptation of the CP birth move initially set to 0 (=1 if birth accepted, 0 otherwise)
  accept = 0

  if(!is.nan(alpha) & u <= alpha){
    
    if(DEBUGLVL == 3) cat("S")
    if(DEBUGLVL == 2) cat("  accepted\n")
    
    ## Acceptation of shift
    accept=1

    ## upate variance, globally (over all segments) - and before calculating the regression coefficient
    proposed.set$sigma2 = updateSigGlobal(proposed.set, X,Y, HYPERvar$delta2, HYPERvar$v0, HYPERvar$gamma0)
    
    ## update the regression coefficient for the new aligned segments
    for(seg.id in c(child1.id, child2.id)) {
      
      ## create regression coefficient for new segment 
      x = extractData(proposed.set, X, seg.id)
      y = extractData(proposed.set, Y, seg.id)
    
      ## sample edge weights
      newB = array(0, Grid.obj$nr.parents + Grid.obj$additional.parents)
      newB[which(proposed.set$edge.struct == 1)] = sampleBxy(x[, which(proposed.set$edge.struct == 1)], y, proposed.set$sigma2, HYPERvar$delta2)

      ## which row in the Regr.Coeff. matrix
      seg.row = which( proposed.set$edge.weights[,1] == seg.id)

      ## replace the old entry for this segment
      proposed.set$edge.weights[seg.row,] = c(seg.id, newB)
    }

    
    if(DEBUGLVL == 2) {
      
      #cat("    accept at ", counter, ", ",  leaf.pair[1], " - ", leaf.pair[2], ", position: ", position, ", axis: ", axis, "\n")
      heatmap(proposed.set$segment.map,Colv=NA,Rowv=NA,scale="none", col=grey.colors(20,start=1,end=0))
      
      z <- try(silent=TRUE, timeout(readline(prompt="Hit me: "), seconds=50))
        
      if(z == "q") { stop("USER EXIT") }
      if(z == "b") { browser() }
    }

    ## assign 
    Grid.obj = proposed.set


  }   else {     ## end accept
    if(DEBUGLVL == 2) cat("  not accepted\n")
  }
  
  
  ##  Return all variables
  ## (+ variable move describing the move type  (1= cut, 2= merge, 3= shift, 4= Update phases)
  return(list(Grid.obj=Grid.obj, accept=accept, move=3, alpha=alpha))

}




###################################################################
# Do the edge moves or just update of model parameters
###################################################################

segment.update <- function(Grid.obj, X, Y, HYPERvar,  DEBUGLVL = 0) {

  if(DEBUGLVL == 1) { cat("START segment.update >\n") }

  ## current number of edges (sign s, subtract bias and SAC)
  nr.edges = sum(Grid.obj$edge.struct) - Grid.obj$additional.parents

  ## expected nr of edges (mean, Lambda)
  mean.nr.edges = rgamma(1, shape= nr.edges + HYPERvar$alphalbd, rate= 1 + HYPERvar$betalbd)
  
  ## Compute acceptation probability vector rho
  rho3 = computeRho3(nr.edges, 0, Grid.obj$smax, HYPERvar$c, mean.nr.edges)
  
  ## Sample u
  u = runif(1, 0, 1)
  
  ## Compute the corresponding move (Edge birth, Edge death or Update the regression coefficient) 
  bduout = edge.move.homogeneousStructure(u, rho3, X, Y, Grid.obj, HYPERvar, DEBUGLVL)

  ## assignments
  Grid.obj$edge.struct = bduout$struct
  Grid.obj$edge.weights = bduout$weights
  Grid.obj$sigma2 = bduout$sigma2
  
  ## (+ variable move describing the move type  (1= CP birth, 2= CP death, 3= CP shift, 4= Update segments)
  return( list( Grid.obj=Grid.obj, move=bduout$move, accept=bduout$accept))

}

##
## edge move function for homogeneous structure
##

edge.move.homogeneousStructure <- function(u, rho3, X, Y, Grid.obj, HYPERvar, DEBUGLVL = 0){

  ## Variable move describing the move type  (4 = update coefficient, 5 = edge birth, 6 = edge death, 7 = edge flip)
  move = 4

  ## Boolean indicating whether the move is accepted or not (=1 if accepted, 0 otherwise, default=0)
  accept = 0

  ## New edges vector, to be returned at the end of the function
  newS = Grid.obj$edge.struct
  
  ## Current number of edges
  nr.edges = sum(Grid.obj$edge.struct) - Grid.obj$additional.parents 

  ## in the case there are fixed edges we need to ignore these (because we are not allowed to operate on them)
  ## This is of course not true for the edge birth move, where we look if the max. nr. of edges is reached
  if(!is.null(Grid.obj$FIXED.INIT.EDGES)) {
    nr.edges = nr.edges - length(Grid.obj$FIXED.INIT.EDGES)
  }
  
  ## Choose between flip move and other moves
  choice = runif(1, 0, 1)

  ## Flip Move
  if(nr.edges > 0 && choice > 0.75) {    

    if(DEBUGLVL == 1)  cat("\n[edges] flip ..") 

    ## flip move is 4
    move = 7

    ## sample from the existing edges
    sample.from = which(Grid.obj$edge.struct[1:Grid.obj$nr.parents]==1)

    ## in the case fixed edges exist..
    if(!is.null(Grid.obj$FIXED.INIT.EDGES)) {

      ## take them out the sample vector
      sample.from = setdiff(sample.from, Grid.obj$FIXED.INIT.EDGES)
    }
    
    ## Sample the original parent, the double vector in sample() is if there is only one edge
    parent.orig = sample(c(sample.from, sample.from),1)  
    
    ## Sample the new parent
    parent.new = sample(c(which(Grid.obj$edge.struct==0), which(Grid.obj$edge.struct==0)), 1) # needed when there is only one position  S==0
        
    stmp = Grid.obj$edge.struct
    stmp[parent.orig] = 0
    stmp[parent.new]  = 1

    rfliplog = 0    

    for(seg.id in getSegmentIDs(Grid.obj)) {
    
      ## create regression coefficient for new segment 
      x = extractData(Grid.obj, X, seg.id)
      y = extractData(Grid.obj, Y, seg.id)
      
      ## number of locations
      omega = length(y)
      
      ## calculate projection matrix
      Pr = computeProjection(as.matrix(x[,which(Grid.obj$edge.struct == 1)]), HYPERvar$delta2)       # current structure 
      Prstar = computeProjection(as.matrix(x[,which(stmp == 1)]), HYPERvar$delta2)     # proposed structure

      ## normal way
      ##	  rflip = rflip * ((gamma0 + t(y) %*% Prstar %*% y)/(gamma0 + t(y) %*% Pr %*% y))^(-(length(y) + v0)/2)
      
      ## log version (better)
      rfliplog = rfliplog + (-(length(y) + HYPERvar$v0)/2) * ( log(HYPERvar$gamma0 + t(y) %*% Prstar %*% y) - log(HYPERvar$gamma0 + t(y) %*% Pr %*% y))
        
    }   

    u = runif(1,0,1)
    
    if(u <= min(1,exp(rfliplog))) {
      accept = 1
      newS = stmp
      if(DEBUGLVL == 1)  cat("accept") 
      if(DEBUGLVL == 3)  cat("f") 

    }
    
  } else {

    not.max.edges = nr.edges < Grid.obj$smax
    
    if(!is.null(Grid.obj$FIXED.INIT.EDGES)) {
        not.max.edges = (length(Grid.obj$FIXED.INIT.EDGES) + nr.edges) < Grid.obj$smax
    }
    
    ##########################
    ## Birth of an edge move
    ##########################
    if(u < rho3[1] && not.max.edges && (length(which(Grid.obj$edge.struct == 0)) > 0 ) ){


      if(DEBUGLVL == 1)  cat("\n[edges] birth ..") 

      ## Variable move describing the move type  (1= Edge birth, 2= Edge death, 3= Update coefficient)
      move = 5
      
      ## Sample the additional edge
      sstar = sample(c(which(Grid.obj$edge.struct==0), which(Grid.obj$edge.struct==0)), 1) # needed when there is only one position  Grid.obj$edge.struct==0
      
      ## Proposed edges vector (with an additional edge)
      stmp = Grid.obj$edge.struct
      stmp[sstar] = 1
      
      ## product over posterior probs.       
      rbirth = 1

      for(seg.id in getSegmentIDs(Grid.obj)) {
    
        if(DEBUGLVL == 2) {  cat("[edge birth] seg.id: ", seg.id) }
          
        ## create regression coefficient for new segment 
        x = extractData(Grid.obj, X, seg.id)
        y = extractData(Grid.obj, Y, seg.id)

        ## calculate projection matrix
        Pr =     computeProjection(as.matrix(x[,which(Grid.obj$edge.struct == 1)]), HYPERvar$delta2)       # current structure 
        Prplus = computeProjection(as.matrix(x[,which(stmp == 1)]), HYPERvar$delta2)     # + 1 edge

        ## number of locations
        omega = length(y)

        ## Compute birth ratio, no log needed because no critical gamma() calculation
        ## orig 1D homog.:
        ##rbirth =    rbirth * ((gamma0 + t(y) %*% Pxlp1 %*% y) /(gamma0 + t(y) %*% Pxl %*% y))^(-(length(y) + v0)/2)/sqrt(1 + HYPERvar$delta2)
        rbirth =rbirth * ((HYPERvar$gamma0 + t(y) %*% Prplus %*% y) / (HYPERvar$gamma0 + t(y) %*% Pr %*% y))^(-(HYPERvar$v0 + omega)/2)/sqrt(1 + HYPERvar$delta2)
        
      }

      ## Sample u 
      u = runif(1,0,1)
      
      tryCatch({ ## BUG ID 3
        if(u <= min(1,rbirth)){
          accept = 1
          newS = stmp

          if(DEBUGLVL == 3)  cat("b") 
        }
      }, error = function(e) {
        write("Error with rbirth in bdu.homogeneous: saving data producing rbirth to DEBUG.DATA/error_rbirth.RData", stderr()) 
        print(e)
        save(HYPERvar$v0, HYPERvar$gamma0, HYPERvar$delta2, y, x, Grid.obj$edge.struct, Grid.obj, file="DEBUG.DATA/error_rbirth.RData")
      })
            	  
     } else {

       
       if(u < rho3[2] & nr.edges > 0 ){  ## makes sure at least one edge exists (despite the bias and SAC edge)

         ##########################
         ## Death of an edge move
         ##########################
         if(DEBUGLVL == 1) { cat("\n[edges] death ..") }

         ## Variable describing the move type  (1 for Edge birth, 2 for Edge death, 3 for Update coefficient)
         move=6

         ## sample from the existing edges
         sample.from = which(Grid.obj$edge.struct[1:Grid.obj$nr.parents]==1)

         ## in the case fixed edges exist..
         if(!is.null(Grid.obj$FIXED.INIT.EDGES)) {
           
           ## take them out the sample vector
           sample.from = setdiff(sample.from, Grid.obj$FIXED.INIT.EDGES)
           
         }
         
         ## Sample the edge to remove
         sstar = sample(c(sample.from, sample.from),1) # needed when there is only one position  S[1:q]==1
         
         ## Proposed edges vector (after taking away one edge)
         stmp = Grid.obj$edge.struct
         stmp[sstar] = 0
         
         rdeath =1
         
         for(seg.id in getSegmentIDs(Grid.obj)) {
           
           if(DEBUGLVL == 2) {  cat("[edge death] seg.id: ", seg.id) }
           
           ## create regression coefficient for new segment 
           x = extractData(Grid.obj, X, seg.id)
           y = extractData(Grid.obj, Y, seg.id)
           
           ## calculate projection matrix
           Pr =      computeProjection(as.matrix(x[,which(Grid.obj$edge.struct == 1)]), HYPERvar$delta2)       # current structure 
           Prminus = computeProjection(as.matrix(x[,which(stmp == 1)]), HYPERvar$delta2)         # - 1 edge
          
           ## complies to TVDBN_SH1D
           rdeath = rdeath * ( (HYPERvar$gamma0 + t(y) %*% Pr %*% y) / (HYPERvar$gamma0 + t(y) %*% Prminus %*% y))^((length(y) + HYPERvar$v0) / 2)*(sqrt(1 + HYPERvar$delta2))
         }
        

         ## Sample u 
         u<-runif(1,0,1)
	  
         if(u <= min(1,rdeath)){
           ## Boolean for the acceptation of the CP death move (=1 if birth accepted, 0 otherwise)
           accept = 1
           newS = stmp

           if(DEBUGLVL == 3)  cat("d") 
           
         }
       } 
     } 
  }

  
  ##
  ## Update Sigma2
  ##
  ## BUG ID 4
  tryCatch({
    
    Sig2_2Dall = updateSigGlobal(Grid.obj, X,Y, HYPERvar$delta2, HYPERvar$v0, HYPERvar$gamma0)
    
  }, error = function(e) {
    
    write("Caught with updateSigGlobal: saving data to DEBUG.DATA/error_bugid4.RData", stderr()) 
    print(e)
    save(Grid.obj, X, Y, Grid.obj$edge.struct, HYPERvar$delta2, file="DEBUG.DATA/error_bugid4.RData")
  })

  ##
  ## Updating weights of each segment (regardless if structure changed)
  ##

  ## create regression parameters from scratch, this will replace the original one fully
  B2Dall = matrix(0, 0, (1 + Grid.obj$nr.parents + Grid.obj$additional.parents)) # +1: seg.id; nr.parent: regr.coeff; additional.parents: bias, sac
  
  
  for(seg.id in getSegmentIDs(Grid.obj)) {
    
    if(DEBUGLVL == 2) {  cat("[regrcoeff update] seg.id: ", seg.id) }
    
    ## create regression coefficient for new segment 
    x = extractData(Grid.obj, X, seg.id)
    y = extractData(Grid.obj, Y, seg.id)
    
    ## sample edge weights
    ## AA22.02.2011, +2 because of additional bias and SAC edge
    newB = array(0, Grid.obj$nr.parents + Grid.obj$additional.parents)
    newB[which(newS == 1)] = sampleBxy(x[, which(newS == 1)], y, Sig2_2Dall, HYPERvar$delta2)
    
    ## set B
    B2Dall = rbind(B2Dall, c(seg.id, newB))
    
  }


  ##  Return all variables
  return(list( u=u, move=move, accept=accept, struct=newS, weights=B2Dall, sigma2=Sig2_2Dall)) 
}












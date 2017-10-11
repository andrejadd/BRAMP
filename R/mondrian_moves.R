



#' Cut a Segment
#' 
#'
#' @importFrom grDevices grey.colors
#' @importFrom stats heatmap
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

    ## assign proposed state to current state
    Grid.obj = proposed.set

    ## do this here because alterations to the tree affect the pointer, if I would have made this on proposed.set before it would
    ## also have affected the tree inside Grid.obj - would require to copy the tree data not the pointer
    ## This way it is easier
    addLeafpair(Grid.obj$mondrian.tree, leaf.id, leaf.id, new.id, child.budget)

    ##
    ## start DEBUGGING
    ##
    if(DEBUGLVL == 2) {

      heatmap(Grid.obj$segment.map,Colv=NA,Rowv=NA,scale="none", col = grey.colors(20,start=1,end=0))
      drawTree(Grid.obj$mondrian.tree)
      
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

  ##  Return all variables
  ## (+ variable move describing the move type  (1= CP birth, 2= CP death, 3= CP shift, 4= Update phases)
  return(list(Grid.obj=Grid.obj, accept=accept, move=1, alpha=alpha, changed.segids=c(leaf.id, new.id)))
  
}



#' Merge a two segments.
#'
#'
#' @importFrom grDevices grey.colors
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
    ## MAYBE FIXME: I could bind edge weights to a Mondrian tree node.. 
    proposed.set$edge.weights = proposed.set$edge.weights[which(proposed.set$edge.weights[,1] != max(child1.id, child2.id)),,drop=F]     # drop=F makes sure the matrix is not transformed to vector when single row is left

    ## overwrite old state segment set
    Grid.obj = proposed.set

    ## delete the two childs
    if(!delLeaf(Grid.obj$mondrian.tree, child1.id)) stop("failed to delete leaf")
    if(!delLeaf(Grid.obj$mondrian.tree, child2.id)) stop("failed to delete leaf")
        
    if(DEBUGLVL == 2) {

      cat(" -> accepted\n")
      heatmap(Grid.obj$segment.map,Colv=NA,Rowv=NA,scale="none", col=grey.colors(20,start=1,end=0))
      drawTree(Grid.obj$mondrian.tree)

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

  
  ##  Return all variables
  ## (+ variable move describing the move type  (1= Segment split, 2= Segm. merge, 3= Segment grow, 4= Update phases)
  return(list(Grid.obj = Grid.obj, accept = accept, move = 2, alpha = alpha, 
              changed.segids = c(min(child1.id, child2.id))))

}


#' Shift an existing cut in a leaf pair.
#' This moves just repositions the cut, i.e. x or y and locations in the same 
#' manner we would set a new cut. 
#' FIXME: Maybe we could create two moves: one that slightly shifts a cut along 
#' his orientation another that changes the orientation: x -> y, y -> x        
#' In all cases, the prior should not change because halfperimeter and cost of 
#' cut stay the same.
#'
#' @importFrom grDevices grey.colors
shift.cut <- function(Grid.obj, X, Y, HYPERvar, DEBUGLVL = 0, counter=0) {

  if (DEBUGLVL == 1) { cat("START segment.merge > \n") }

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
  return(list(Grid.obj=Grid.obj, accept=accept, move=3, alpha=alpha, changed.segids=c(child1.id, child2.id)))

}









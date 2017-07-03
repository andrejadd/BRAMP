

##
## Split and Merge are straight forward
##
## When growing a structure:
## 1. get adjacent segments
##           getAdjacentSegments(newobj, segid)
## 2. check which neighbors are too small to grow into
##          getNrElements(segset, segid)
##
## 3. growSegment(newobj, segid, segid.into)
##
## extractData() 
##
##

batch.test <- function() {

  newobj = createSegments(20,20)

  newobj = splitSegment(newobj, 1, "x", 8)

  newobj = splitSegment(newobj, 2, "y", 6)

  newobj = mergeSegment(newobj, 3, 1)

  newobj = splitSegment(newobj, 2, "x", 4)

  newobj = splitSegment(newobj, 2, "y", 13)


  ## get adjacent segment IDs
  newobj = getAdjacentSegments(newobj, 2)


  return(newobj)

}

createGrid <- function(xlocs, ylocs, minSegsize, start.budget = 1, additional.parents,
                       nr.parents, smax, FIXED.INIT.EDGES) {

  
  data.obj = list(
    xlocs=xlocs,
    ylocs=ylocs,
    minSegsize = minSegsize, 
    segment.map=NULL,
    mondrian.tree=NULL,
    valid.locs.map=NULL,
    edge.weights = NULL,
    additional.parents=additional.parents,
    nr.parents = nr.parents,
    smax = smax,
    FIXED.INIT.EDGES = FIXED.INIT.EDGES
    )

  ## the location matrix, where each location/element has a segment id assigned to it
  data.obj$segment.map = matrix(1, nrow=ylocs, ncol=xlocs)

  ## init the first root node of the Mondrian tree
  data.obj$mondrian.tree = newNode(1, 1, start.budget)

  ## later set locations to '0' if not valid/sampled, these will be ignored when counting segment members etc. (speeds up things)
  data.obj$valid.locs.map = matrix(1, nrow=ylocs, ncol=xlocs)
  
  return(data.obj)

}

initMondrian <- function(data.obj) {

  leaf.ids = c(1)

  cat("\n")

  ## do nothing
  return(data.obj)
  
  while(length(leaf.ids) > 0) {

    ## put new leafs in here, will overwrite leaf.ids at the end
    ## as long as new.leafs is filled below, the main while loop will not exit
    new.leafs = c()
    
    ## do cut operations as long as there are leaves and cuts are above the budget
    for(leaf.id in leaf.ids) {

      cat("\nCUT segment ", leaf.id, " .. ")
      
      ## try to split the leave
      ## get length and width
      seg.dim = getSegmentDim(data.obj, leaf.id)
      
      cost = expdist(sum(seg.dim))     ## is the random.expovariate() -> see random.py 

      ## lookup budget for the block
      budget = data.obj$budget.mat[ which(data.obj$budget.mat[,1] == leaf.id), 2]

#      cat(", budget: ", budget, ", length: ", seg.dim[2], ", width: ", seg.dim[1], "\n")
      
      
      ## check if a cut can be done at all
      if(cost <= budget) {

        ## decide along which axis to cut
        axis="y"                                              ## y axis
        if(runif(1) < (seg.dim[1]/(sum(seg.dim)))) axis="x"   ## x axis

        ## get cut position 
        positions = getSplitPositions(data.obj, leaf.id, axis)
        
        ## extract random position
        position = sample(c(positions,positions), 1)

        ## do the cut
        proposed.set = splitSegment(data.obj, parent.id=leaf.id, axis=axis, position=position, parent.budget=budget, child.budget = (budget - cost) )

    
        new.id = proposed.set$last.added.id
        
        ## check if elements are large enough
        if(getNrElements(proposed.set, leaf.id) < proposed.set$minSegsize) { cat(" too small"); next }
        if(getNrElements(proposed.set, new.id) < proposed.set$minSegsize)  { cat(" too small"); next }

        cat(" ok\n")
        data.obj = proposed.set

        ## add regression coeff. 
        
        
   #     cat("  accept split at ", axis, ":", position, " -> new.id: ", new.id, "\n")
                
        ## push the id into the leaf vector
        new.leafs = c(new.leafs, leaf.id,new.id)
      } else {
        cat(" bad budget")
      }

      
    } ## end for

    leaf.ids = new.leafs
    
  #  heatmap(data.obj$segment.map,Colv=NA,Rowv=NA,scale="none", col=grey.colors(20,start=1,end=0))

  #  cat("INSIDE NODES:\n")
  #  print.table(data.obj$inside.node.pairs)
  #  cat("LEAVE PAIRS:\n")
  #  print.table(data.obj$leave.pairs)

  #  cat("\n-> on leave.ids stack: ")
  #  print.table(leave.ids)

#    browser()
  } ## end while leaf.ids

  
  return(data.obj)
}


getRandomLeafPair <- function(data.obj) {  ## tree.ok returns a leaf pair that can be merged

  pairs = getLeafPairIDs(data.obj$mondrian.tree)

  ## if no pairs where found
  if(nrow(pairs) == 0) {
    return(NULL) 
  }

  pair.row = sample(1:nrow(pairs), 1:nrow(pairs), 1)

  leaf.id1 = pairs[pair.row, 1]
  leaf.id2 = pairs[pair.row, 2]

#  cat("[segments.R:getRandomLeafPair] leaf.id1: ", leaf.id1, ", leaf.id2: ", leaf.id2, "\n")
  
  return(list(leaf1 = getLeafNode(data.obj$mondrian.tree, leaf.id1), leaf2 = getLeafNode(data.obj$mondrian.tree, leaf.id2)) )
  
}
  
expdist <- function(mean) {   ## my exponential distr. sampler

  rnd.samp = runif(1)
  
  draw.exp = 1/mean * exp(-1/mean * rnd.samp)

  draw.log =  -log(1.0 - rnd.samp) / mean

  return(draw.log)
}


splitSegment <- function(data.obj=NULL, parent.id=NULL, new.id=NULL, axis=NULL, position=NULL) {

  ## get locations of elements belonging to the parent segment
  x.ij = which(data.obj$segment.map == parent.id, arr.ind=T)

  ## get lowest none-used segment id that is larger than the parent.id (the compare vector starts width the parent.id)
  ## the last constraint is important because the larger id will be discarded with a merge
  compare = c(parent.id:(max(data.obj$segment.map) + 1))

  ## if not set (default) find one, otherwise only work on the segment.map and do nothing on the tree
  ## new.id is only provided for the shift.cut move, where we want to make sure to keep the same child id
  if(is.null(new.id)) {
    new.id = min(setdiff(compare, unique(as.vector(data.obj$segment.map))))
  } 

  ## set one of the new segments create by the split to the new id
  if(axis == "x") {
    
    data.obj$segment.map[ x.ij[ which(x.ij[,2] < position),]] = new.id

  } else {

    data.obj$segment.map[ x.ij[ which(x.ij[,1] < position),]] = new.id
  }

  data.obj$last.added.id = new.id

  return(data.obj)
}

mergeSegment <- function(data.obj=NULL, seg.ids=c()) {

  ## get the coordinates of both segments
  x1.ij = which(data.obj$segment.map == seg.ids[1], arr.ind=T)
  x2.ij = which(data.obj$segment.map == seg.ids[2], arr.ind=T)

  ## extract smallest id, which is also the parents id
  parent.id = min(seg.ids)
  
  ## assign the same (smallest) segment id
  data.obj$segment.map[x1.ij] = parent.id
  data.obj$segment.map[x2.ij] = parent.id

  return(data.obj)
  
}

extractData <- function(data.obj=NULL, X=NULL, seg.id) { ## tree.ok

  ## get all locations
  ## FIXME: the transpose is a workaround to get indices to proper index into X and Y, must avoid too much operations somehow
  x.ij = which( t(data.obj$segment.map) == seg.id & t(data.obj$valid.locs.map) == 1 )

  ## if its a vector, it is the target data from which we extract
  if(is.vector(X)) return(X[x.ij])

  ## otherwise its the full data vector
  return( X[x.ij,] )

}


##
## note, it does not grow diagonal but it shouldn't be a big problem , watch this!
##

growSegment <- function(data.obj=NULL, seg.id1, seg.id2) {

  ## We need the data structure growed.coord that was already used in adjacent neighbor (grows one cell)
  ##   
  if(data.obj$growed.seg.id != seg.id1) {
    stop("FIXME: it appears getAdjacentSegments() was not called before for this segment id, I could do now but something is wrong in the code flow")
  }

  ## get the coordinates of segment to grow into
  seg2.coord = which(data.obj$segment.map == seg.id2, arr.ind=T)

  unity = rbind(seg2.coord, data.obj$growed.coord)

  ## these are the position to expand seg.id1 into
  expand.into = unity[duplicated(unity),]

  ## set the cells to the expanded segments id
  data.obj$segment.map[ expand.into] = seg.id1

  return(data.obj)
}


#
# return all possible positions that produce valid offspring segments
#


getSplitPositions <- function(data.obj, seg.id, axis) { ##tree.ok

  x.ij = which(data.obj$segment.map == seg.id, arr.ind=T)

  
  if(axis == "x") { ## these are the columns
    smallest = min(x.ij[,2])
    largest = max(x.ij[,2])
  } else {
    smallest = min(x.ij[,1])
    largest = max(x.ij[,1])
  }

  ## check to avoid overflow below
  if(smallest == largest) {
    return(c(smallest))
  }
  
  return((smallest + 1):largest)

}


getSegmentIDs <- function(data.obj=NULL) {

  return(unique(as.vector(data.obj$segment.map)))
}

getNrElements <- function(data.obj=NULL, seg.id=NULL) { ## tree.ok

  ## if no segment id specified, return total nr. of locations/elements in grid
  if(is.null(seg.id)) return(sum(data.obj$valid.locs.map == 1))
  
  ## the second term is the map indicating valid positions (1), both sites produce a TRUE/FALSE matrix which is compared logically
  return(sum( data.obj$segment.map == seg.id & data.obj$valid.locs.map == 1  ))
}

getNrSegments <- function(data.obj=NULL) { ##tree.ok

  return(length(unique(as.vector(data.obj$segment.map))))

}

getNrLeafs <- function(data.obj=NULL) {
  
  nrow(data.obj$leaf.pairs)
}

##
##FIXME, not used , should check use of xlocs and ylocs
##
expandSegmentBoundaries <- function(data.obj=NULL, seg.id) {

  seg = which(data.obj$segment.map == seg.id, arr.ind=T)

  grow.coord = matrix(0,nrow=0, ncol=2)
  
  shift = seg
  shift[,2] = seg[,2] + 1   # x right (col)
  grow.coord = rbind(grow.coord, shift)
  
  shift = seg
  shift[,2] = seg[,2] - 1   # x left  (col)
  grow.coord = rbind(grow.coord, shift)

  shift = seg
  shift[,1] = shift[,1] + 1   # y up (row)
  grow.coord = rbind(grow.coord, shift)

  shift = seg
  shift[,1] = shift[,1] - 1   # y down (row)
  grow.coord = rbind(grow.coord, shift)

  grow.coord = grow.coord[ grow.coord[,1] > 0 & grow.coord[,1] <= data.obj$ylocs & grow.coord[,2] > 0 & grow.coord[,2] <= data.obj$xlocs ,]

  ## remove duplicates
  grow.coord = grow.coord[!duplicated(grow.coord),]

  return(grow.coord)
}


##
## return only IDs of adjacent segments
##

getAdjacentIDs <- function(data.obj=NULL, seg.id) {

  ## need this to find out bordering segments
  growed.coord = expandSegmentBoundaries(data.obj, seg.id)

  ## extract actual ids 
  ids = data.obj$segment.map[ growed.coord ]

  ## these are the neigbors
  return( setdiff(unique(ids), seg.id))
}

##
## return object including adjacent IDs and coordinates of growed segment
##

getAdjacentSegments <- function(data.obj=NULL, seg.id) {

  growed.coord = expandSegmentBoundaries(data.obj, seg.id)

  ids = data.obj$segment.map[ growed.coord ]

  ## these are the actual neigbors
  data.obj$adjacent.ids = setdiff(unique(ids), seg.id)
    
  ## save, used for the subsequent step to grow the segment
  data.obj$growed.coord = growed.coord

  ## save the current id to be sure of consistence when using later $growed.coord
  data.obj$growed.seg.id = seg.id

  return(data.obj)
  
}


setSegmentBudget <- function(data.obj, parent.id, child.ids, parent.budget, childs.budget) {  # save the budget for the parent and child blocks

  ## move the budget of the parent segment to the inner node budgets
  ## first delete
  data.obj$budget.mat = data.obj$budget.mat[ -which(data.obj$budget.mat[,1] == parent.id), ]

  ## add to inside.node.budget
  data.obj$inside.node.budget = rbind(data.obj$inside.node.budget, c(parent.id, parent.budget))

  ## finally add the child budgets to the leaf budget matrix (budget.mat)
  data.obj$budget.mat = rbind(data.obj$budget.mat, c(child.ids[1], childs.budget))
  data.obj$budget.mat = rbind(data.obj$budget.mat, c(child.ids[2], childs.budget))

  return(data.obj)

}

removeSegmentBudget <- function(data.obj, parent.id, child.id) {  # save the budget for the parent and child blocks


  ## remove the budget information for the child ids
#  tryCatch({
    data.obj$budget.mat = data.obj$budget.mat[ -which(data.obj$budget.mat[,1] == parent.id), ,drop=F]
    data.obj$budget.mat = data.obj$budget.mat[ -which(data.obj$budget.mat[,1] == child.id), ,drop=F]
#  }, error = function(e) {
#    print(e)
#    browser()
#  } )

           ## Move the budget of the parent from the inside.node.budget matrix to the budget.mat (leaf node) matrix
  ## first copy parent budget, for this reason I need to get the lowest budget for the parent
  
  ## get lowest budget

    tryCatch( {
      parent.budget = min(data.obj$inside.node.budget[ which(data.obj$inside.node.budget[,1] == parent.id) , 2])
    }, warning = function(w) {

      cat("CAUGHT WARNING WHILE FETCHING INSIDE BUDGET\n")
      browser()

    }

             )
  ## I remove this entry with the lowest id (if there are higher entries it means these are innern nodes)
  data.obj$inside.node.budget = data.obj$inside.node.budget[ -which(data.obj$inside.node.budget[,2] == parent.budget & data.obj$inside.node.budget[,1] == parent.id), , drop=F]
  
  ## append to leaf budget
  data.obj$budget.mat = rbind(data.obj$budget.mat, c(parent.id, parent.budget))
  
  return(data.obj)

}


getSegmentDim <- function(data.obj, seg.id) {  ## return the length and weight of the segment

  seg.coords = which(data.obj$segment.map == seg.id, arr.ind=T)

  x.length = max(seg.coords[,2]) - min(seg.coords[,2]) + 1
  y.length = max(seg.coords[,1]) - min(seg.coords[,1]) + 1

  return(c(x.length, y.length))

}

getPairDim <- function(data.obj, seg.ids) {  ## return the length and weight of two segments

  ## extract for both segment ids
  seg.coords1 = which(data.obj$segment.map == seg.ids[1], arr.ind=T)
  seg.coords2 = which(data.obj$segment.map == seg.ids[2], arr.ind=T)

  ## combine into one matrix
  seg.coords = rbind(seg.coords1, seg.coords2)
  
  x.length = max(seg.coords[,2]) - min(seg.coords[,2]) + 1
  y.length = max(seg.coords[,1]) - min(seg.coords[,1]) + 1

  return(c(x.length, y.length))
}


getSegmentDimScale <- function(data.obj, seg.id) {  ## sale the length and weight into interval [0,1] for the continuous Mondrian

  ## extract the length and width
  if(length(seg.id) == 1) {
    ## get it for a single segment
    dims = getSegmentDim(data.obj, seg.id)
  } else {
    ## get it for two segments, did not have time to find out the perimeter of a combined rectangular (maybe a FIXME?)
    dims = getPairDim(data.obj, seg.id)
  }

  ## scale x and y
  ## FIXME: for irregular grids the smaller side scales also into [0,1] but need to add a side ratio that penalized the shorter side
  x.scaled = dims[1] / data.obj$xlocs
  y.scaled = dims[2] / data.obj$ylocs
  
  return(c(x.scaled, y.scaled))

}

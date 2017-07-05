



###################################################################
# Do the edge moves or just update of model parameters
###################################################################

edge_moves <- function(Grid.obj, X, Y, HYPERvar,  DEBUGLVL = 0) {
  
  if(DEBUGLVL == 1) { cat("START segment.update >\n") }
  
  ## current number of edges (sign s, subtract bias and SAC)
  nr.edges = sum(Grid.obj$edge.struct) - Grid.obj$additional.parents
  
  ## expected nr of edges (mean, Lambda)
  mean.nr.edges = rgamma(1, shape= nr.edges + HYPERvar$alphalbd, rate= 1 + HYPERvar$betalbd)
  
  ## Compute acceptation probability vector rho
  rho3 = computeRho3(nr.edges, 0, Grid.obj$smax, HYPERvar$c, mean.nr.edges)
  
  ## Sample u
  u = runif(1, 0, 1)
 
  ## Variable move describing the move type  (4 = update coefficient, 5 = edge birth, 6 = edge death, 7 = edge flip)
  move = 4
  
  
  ## Boolean indicating whether the move is accepted or not (=1 if accepted, 0 otherwise, default=0)
  accept = 0
  
  
  ## New edges vector, to be returned at the end of the function
  newS = Grid.obj$edge.struct
  
  
  ## Current number of edges
  nr.edges = sum(Grid.obj$edge.struct) - Grid.obj$additional.parents 
  
  
  ## In the case there are fixed edges we need to ignore these (because we are not allowed to operate on them)
  ## This is of course not true for the edge birth move, where we look if the max. nr. of edges is reached.
  if(!is.null(Grid.obj$fixed_edges)) {
    nr.edges = nr.edges - length(Grid.obj$fixed_edges)
  }
  
  
  ## Choose between flip move and other moves.
  choice = runif(1, 0, 1)
  
  
  ##
  ## Exchange two parents (add one, delete one at the same time), i.e. flip two edges.
  ##
  
  ## sample from the existing edges
  sample.from = which(Grid.obj$edge.struct[1:Grid.obj$nr.parents] == 1)
  
  
  ## in the case fixed edges exist..
  if(!is.null(Grid.obj$fixed_edges)) {
    
    ## Take out the fixed edges from the active edges.
    sample.from = setdiff(sample.from, Grid.obj$fixed_edges)
  }
  
  
  ## Make flip if at least one (non-fixed) edge exists
  if((choice > 0.75) && length(sample.from) > 0) {    
    
    if(DEBUGLVL == 1)  cat("\n[edges] flip ..") 
    
    ## Variable move describing the move type.  
    ##   5 : edge birth
    ##   6 : edge death
    ##   7 : edge flip
    move = 7
    
    
    ## Sample the edge that is turned off.
    delete_edge = sample(c(sample.from, sample.from),1)  
    
    
    ## Sample the new parent
    activate_edge = sample(c(which(Grid.obj$edge.struct==0), which(Grid.obj$edge.struct==0)), 1)
    
    
    ## Get structure and make the flip.
    stmp = Grid.obj$edge.struct
    stmp[delete_edge] = 0
    stmp[activate_edge]  = 1
    
    
    ## Acceptance probability of a flip aggregated for all segments.
    rfliplog = 0    
    
    for(seg.id in getSegmentIDs(Grid.obj)) {
      
      ## create regression coefficient for new segment 
      X_tmp = extractData(Grid.obj, X, seg.id)
      y_tmp = extractData(Grid.obj, Y, seg.id)
      
      ## The number of observations in this segment.
      omega = length(y_tmp)
      
      ## Calculate the projection matrices given the old and new edge structure.
      Pr = computeProjection(as.matrix(X_tmp[,which(Grid.obj$edge.struct == 1)]), HYPERvar$delta.snr)       # current structure 
      Prstar = computeProjection(as.matrix(X_tmp[,which(stmp == 1)]), HYPERvar$delta.snr)     # proposed structure
      
      ## normal way
      ##	  rflip = rflip * ((beta.var + t(y) %*% Prstar %*% y)/(beta.var + t(y) %*% Pr %*% y))^(-(length(y) + alpha.var)/2)
      
      ## Add up acceptance probability: log version (better, more robust)
      rfliplog = rfliplog + (-(length(y_tmp) + HYPERvar$alpha.var)/2) * ( log(HYPERvar$beta.var + t(y_tmp) %*% Prstar %*% y_tmp) - log(HYPERvar$beta.var + t(y_tmp) %*% Pr %*% y_tmp))
      
    }   
    
    u = runif(1,0,1)
    
    if(u <= min(1,exp(rfliplog))) {
      accept = 1
      newS = stmp
      if(DEBUGLVL == 1)  cat("accept") 
      if(DEBUGLVL == 3)  cat("f") 
      
    }
    
  } else {
    
    
    ##
    ## Add a parent, i.e. birth of an edge.
    ##  Check if the move is possible, i.e. 
    ##           1. The number of active edges is smaller than the fanin (smax)
    ##           2. At least one in-active edge exists      
    ##                            
    if( (u < rho3[1]) && 
       (nr.edges < Grid.obj$smax) && 
       (length(which(Grid.obj$edge.struct == 0)) > 0 ) 
       ){
      
      
      if(DEBUGLVL == 1)  cat("\n[edges] birth ..") 
      
      
      ## Variable move describing the move type  
      ##   5 : edge birth
      ##   6 : edge death
      ##   7 : edge flip
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
        Pr =     computeProjection(as.matrix(x[,which(Grid.obj$edge.struct == 1)]), HYPERvar$delta.snr)       # current structure 
        Prplus = computeProjection(as.matrix(x[,which(stmp == 1)]), HYPERvar$delta.snr)     # + 1 edge
        
        ## number of locations
        omega = length(y)
        
        ## Compute birth ratio, no log needed because no critical gamma() calculation
        ## orig 1D homog.:
        ##rbirth =    rbirth * ((beta.var + t(y) %*% Pxlp1 %*% y) /(beta.var + t(y) %*% Pxl %*% y))^(-(length(y) + alpha.var)/2)/sqrt(1 + HYPERvar$delta.snr)
        rbirth =rbirth * ((HYPERvar$beta.var + t(y) %*% Prplus %*% y) / (HYPERvar$beta.var + t(y) %*% Pr %*% y))^(-(HYPERvar$alpha.var + omega)/2)/sqrt(1 + HYPERvar$delta.snr)
        
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
        save(HYPERvar$alpha.var, HYPERvar$beta.var, HYPERvar$delta.snr, y, x, Grid.obj$edge.struct, Grid.obj, file="DEBUG.DATA/error_rbirth.RData")
      })
      
    } else {
      
      ##
      ## Delete a parent, i.e. death of an edge.
      ##
      
      ## Get the set of active edges.
      sample.from = which(Grid.obj$edge.struct[1:Grid.obj$nr.parents] == 1)
      
      
      ## Take out the fixed_edges from the available set of edges.
      if(!is.null(Grid.obj$fixed_edges)) {
        
        ## take them out the sample vector
        sample.from = setdiff(sample.from, Grid.obj$fixed_edges)
        
      }
      
      
      ##
      ## Execute death move if at least one edge exists in 'sample.from'
      ##
      if(u < rho3[2] && length(sample.from) > 0){  
        
        
        if(DEBUGLVL == 1) { cat("\n[edges] death ..") }
        
        
        ## Variable move describing the move type  
        ##   5 : edge birth
        ##   6 : edge death
        ##   7 : edge flip
        move=6
        
        
        ## Sample the edge to remove.
        sstar = sample(c(sample.from, sample.from),1) 
        
        
        ## Set the edge to zero, i.e. turn the edge off.
        stmp = Grid.obj$edge.struct
        stmp[sstar] = 0
        
        
        ## Acceptance probability - updated over the segments.
        rdeath = 1
        
        
        ## Loop over each segment.
        for(seg.id in getSegmentIDs(Grid.obj)) {
          
          if(DEBUGLVL == 2) {  cat("[edge death] seg.id: ", seg.id) }
          
          
          ## Get the design matrix and response vector for this segment. 
          X_tmp = extractData(Grid.obj, X, seg.id)
          y_tmp = extractData(Grid.obj, Y, seg.id)
          
          
          ## calculate projection matrix
          Pr =      computeProjection(as.matrix(X_tmp[,which(Grid.obj$edge.struct == 1)]), HYPERvar$delta.snr)       # current structure 
          Prminus = computeProjection(as.matrix(X_tmp[,which(stmp == 1)]), HYPERvar$delta.snr)                       # structure with removed edge
          
          
          ## Take product of acceptance probabilities.
          rdeath = rdeath * ( (HYPERvar$beta.var + t(y_tmp) %*% Pr %*% y_tmp) / (HYPERvar$beta.var + t(y_tmp) %*% Prminus %*% y_tmp))^((length(y_tmp) + HYPERvar$alpha.var) / 2)*(sqrt(1 + HYPERvar$delta.snr))
        }
        
        
        ## Sample u 
        u<-runif(1,0,1)
        
        ## If smaller, accept the move and save the new structure into newS.
        if(u <= min(1,rdeath)){
          
          accept = 1
          newS = stmp
          
          if(DEBUGLVL == 3)  cat("d") 
          
        }
      } 
    } 
  }
  

  ## Save the latest edge structure.
  Grid.obj$edge.struct = newS 
  
  return( list( Grid.obj = Grid.obj, move = move, accept = accept))
  
}







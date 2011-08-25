#####################################################################################
## BIRTH OF A CHANGEPOINT
#####################################################################################

cp.birth <- function(ALTERX, XE, YE, S2Dall, B2Dall, Sig2_2Dall, X, Y, D, GLOBvar, HYPERvar, DEBUGLVL1 = F, DEBUG_BIRTH_EXT = F) {
  
  ### assignement of global variables used here ###
  q = GLOBvar$q
  qmax = GLOBvar$qmax
  minPhase = GLOBvar$minPhase
  smax = GLOBvar$smax
  dyn = GLOBvar$dyn
  XMphase = GLOBvar$XMphase
  YMphase = GLOBvar$YMphase
  birth_proposal = GLOBvar$birth_proposals
  xlocs = GLOBvar$xlocs
  ylocs = GLOBvar$ylocs
  
  ### end assignement ###

  ### assignement of hyperparameters variables used here ###
  alphalbd = HYPERvar$alphalbd
  betalbd = HYPERvar$betalbd
  alphad2 = HYPERvar$alphad2
  betad2 = HYPERvar$betad2
  v0 = HYPERvar$v0
  gamma0 = HYPERvar$gamma0
  delta2 = HYPERvar$delta2

  ### end assignement ###

  if(DEBUG_BIRTH_EXT == TRUE) { cat("\n START CP.BIRTH\n\n") }

  ## assign the changepoint vector of interest
  if(ALTERX) { E = XE; E.other = YE } else { E = YE; E.other = XE } 
  
  ## search for possible CP, not in E and not close to E if minPhase (length of phase) is > than 1
  toremove = E
  if(minPhase>1) for(i in 1:(minPhase-1)) toremove = c(toremove, E-i, E+i)

  ## possible CPs are those not in 'toremove'
  possibleCP = setdiff((1+dyn):E[length(E)], toremove)

  ## check if a valid cp position exists or sample() will fail
  if(length(possibleCP)== 0) {
    return(list(XE=XE, YE=YE, S2Dall=S2Dall, B2Dall=B2Dall, Sig2_2Dall=Sig2_2Dall, accept=0, move=1, alpha=0, estar=-1))
  }

  ## sample uniformly new cp, returns NaN in the case that all cp produce invalid segments (constraint to min.seglocs)
  cp.new = sampleValidateCPs(candidateCPs=possibleCP, min.seglocs = minPhase*minPhase, E, E.other, Y, ALTERX, xlocs, type="birth") 
  
  ## if no valid changepoint could be found, return
  if(is.nan(cp.new)) {
    #cat("giving up new cp\n") 
    return(list(XE=XE, YE=YE, S2Dall=S2Dall, B2Dall=B2Dall, Sig2_2Dall=Sig2_2Dall, accept=0, move=1, alpha=0, estar=-1))
  }
 

  ## Create next state CP vector
  Eplus = sort(c(E,cp.new))
    
  ## Position of the phase containing the new CP
  poskl = sum(E < cp.new)

  alpha = cp.computeAlpha(1, X, Y, xlocs, ylocs, ALTERX, XMphase, YMphase, E, Eplus, E.other,  poskl, HYPERvar, S2Dall, D, DEBUG_BIRTH_EXT) 
  
  ## Sample u to conclude either to  acceptation or to rejection
  u = runif(1,0,1)
  
  ## Boolean for the acceptation of the CP birth move initially set to 0 (=1 if birth accepted, 0 otherwise)
  accept = 0
  
  if(!is.nan(alpha) & u <= alpha){
    ## Acceptation of the birth of the new CP
    ## Move acceptation boolean = 1
    accept=1

    ## I will assume that the new segment is the right one. this is simplified but if to follow the 1D approach one
    ## would have to choose among more than 2 neighbors for the new segment. This can be of course done in some sort of local shuffle
    ## However, the structure stays the same and the regress. coefficient and variance (sigma2) are updated in any case in later iterations
    ## We assume poskl to hold the "old" segment and poskl+1 the new

    ## First increment the segment ids of segments that are greater than the inserted segments
    if(ALTERX) { colmn = 1 } else { colmn = 2}  # select the column

    B2Dall[which(B2Dall[,colmn] > poskl),colmn] = B2Dall[which(B2Dall[,colmn] > poskl),colmn] + 1

    ## sets the new segments segment id, and updates the new changepoint vector
    if(ALTERX) { xsegid = poskl+1; XE = Eplus; nrOtherSegs = length(YE) - 1  } else { ysegid = poskl+1; YE = Eplus; nrOtherSegs = length(XE) - 1 }

    ## upate variance, globally (over all segments) - and before calculating the regression coefficient
    Sig2_2Dall = updateSigGlobal(xlocs,XMphase, YMphase, XE, YE, X,Y, S2Dall, delta2, HYPERvar$v0, HYPERvar$gamma0)
       

    ## update regression coefficient  
    for(i in 1:nrOtherSegs) {

      if(ALTERX) { ysegid = i; } else { xsegid = i; }
   
      ## get segment coordinates
      segcoord = c(XMphase[XE[xsegid]], YMphase[YE[ysegid]],XMphase[XE[xsegid+1]]-1, YMphase[YE[ysegid+1]]-1)

      if(DEBUG_BIRTH_EXT == TRUE) { cat("xsegid: ", xsegid, ", ysegid: ", ysegid, "\n");
                                    cat("prodPhiPlus: segcoord ")
                                    print.table(segcoord) }
      
      
      x = extractNodes(X, segcoord, xlocs,F)
      y = extractNodes(Y, segcoord, xlocs,F)
      Pr = computeProjection(as.matrix(x[,which(S2Dall == 1)]), delta2)

      ## regression coefficient
      ## AA22.02.2011, takes the additional SAC edge into account
      newB = array(0, q + 2)      

      newB[which(S2Dall == 1)] = sampleBxy(x[,which(S2Dall==1)], y, Sig2_2Dall, delta2)
      B2Dall = rbind(B2Dall, c(xsegid, ysegid, newB))
    }
  } 

  ##  Return all variables
  ## (+ variable move describing the move type  (1= CP birth, 2= CP death, 3= CP shift, 4= Update phases)
  return(list(XE=XE, YE=YE, S2Dall=S2Dall, B2Dall=B2Dall, Sig2_2Dall=Sig2_2Dall, accept=accept, move=1, alpha=alpha, estar=cp.new))

}


#####################################################################################
#DEATH OF A CHANGEPOINT
#####################################################################################

cp.death <- function(ALTERX, XE, YE, S2Dall, B2Dall, Sig2_2Dall, X, Y, D, GLOBvar, HYPERvar, DEBUGLVL1 = F, DEBUG_BIRTH_EXT = F) {

  ### assignement of global variables used here ###
  q = GLOBvar$q
  smax = GLOBvar$smax
  qmax = GLOBvar$qmax
  xlocs = GLOBvar$xlocs
  ylocs = GLOBvar$ylocs
  XMphase = GLOBvar$XMphase
  YMphase = GLOBvar$YMphase

### end assignement ###

  ### assignement of hyperparameters variables used here ###
  alphad2 = HYPERvar$alphad2
  betad2 = HYPERvar$betad2
  v0 = HYPERvar$v0
  gamma0 = HYPERvar$gamma0
  ### end assignement ###

  if(DEBUG_BIRTH_EXT == TRUE) { cat("\n START CP.DEATH\n\n") }

  if(ALTERX) { E = XE; E.other = YE } else { E = YE; E.other = XE }
  
  ## check if there are at least one possible CP (should never happen that this fct gets selected for run but just to prevent worsed case and crash)
  if(length(E) < 3) {
    return(list(XE=XE, YE=YE, S2Dall=S2Dall, B2Dall=B2Dall, Sig2_2Dall=Sig2_2Dall, accept=0, move=2, alpha=0, estar=-1))
  }
  
  ## Sample the CP to be removed
  estar = sample(c(E[2:(length(E)-1)], E[2:(length(E)-1)]), 1)

  ##  Position of the phase starting at the selected CP
  poskstar = sum(E <= estar)

  ## CP changepoint vector without the CP
  Eminus = E[-which(E == estar)]
  
  ## compute the acceptance probability
  alpha = cp.computeAlpha(-1, X, Y, xlocs, ylocs, ALTERX, XMphase, YMphase, E=Eminus, Eplus=E, E.other,  poskl=(poskstar-1), HYPERvar, S2Dall, D, DEBUG_BIRTH_EXT) 

  ## Sample u to conclude either to  acceptation or to rejection
  u = runif(1,0,1)

  ## Boolean for the acceptation of the CP death move initially set to 0 (=1 if birth accepted, 0 otherwise)
  accept = 0
  
  if(!is.nan(alpha) & u <= alpha){

    ## Acceptation of the death of the selected CP, Move acceptation boolean =1
    accept=1

    ## update the CP vector
    if(ALTERX) { XE = Eminus} else { YE = Eminus}

    ## AA: Need to decide which parameter segments to keep, naturally and simply would be (poskstar - 1) - this segment id would stay
    ## when merging [(poskstar - 1),poskstar], although, we could allow the parameters of the larger segment to stay.
    ## But anyways, update of the parameters will take place later. This might only influence convergence a little

    ## delete the parameters with segid poskstar
    if(ALTERX) { colmn = 1 } else { colmn = 2}  # select the column
    B2Dall =     B2Dall[which(B2Dall[,colmn] != poskstar),,drop=F]     # drop=F makes sure the matrix is not transformed to vector when single row is left
    
    ## decrement the segments that are greater poststar in order to fill the gap of the now morged segments
    B2Dall[which(B2Dall[,colmn] > poskstar),colmn] = B2Dall[which(B2Dall[,colmn] > poskstar),colmn] - 1

    # no need to update the variance because only needed for calculating the regression coefficient B
    # Sig2_2Dall = updateSigGlobal(XE, YE, X,Y, S2Dall, delta2, HYPERvar$v0, HYPERvar$gamma0)
  }

  ##  Return all variables
  ## (+ variable move describing the move type  (1= CP birth, 2= CP death, 3= CP shift, 4= Update phases)
  return(list(XE=XE, YE=YE, S2Dall=S2Dall, B2Dall=B2Dall, Sig2_2Dall=Sig2_2Dall, accept=accept, move=2, alpha=alpha, estar=estar))
}

                                  

#####################################################################################
#MOVE OF A CHANGEPOINT
#####################################################################################


cp.shift <- function(ALTERX, XE, YE, S2Dall, B2Dall, Sig2_2Dall, X, Y, GLOBvar, HYPERvar, DEBUGLVL1 = F, DEBUG_BIRTH_EXT = F) {
  
  ### assignement of global variables used here ###
  q = GLOBvar$q
  minPhase = GLOBvar$minPhase
  smax = GLOBvar$smax
  xlocs = GLOBvar$xlocs
  ylocs = GLOBvar$ylocs
  XMphase = GLOBvar$XMphase
  YMphase = GLOBvar$YMphase
  
  ### end assignement ###

  ### assignement of hyperparameters variables used here ###
  alphad2 = HYPERvar$alphad2
  betad2 = HYPERvar$betad2
  v0 = HYPERvar$v0
  gamma0 = HYPERvar$gamma0
  delta2 = HYPERvar$delta2
  ### end assignement ###
  
  
  ## Select two segments that have adjacent boundaries, this boundary corresponds to a CP that can be shifted
  ## Makes use of the segment pair search cp.death

  if(DEBUG_BIRTH_EXT == TRUE) { cat("\n START CP.SHIFT\n\n") }

  ## assign changepoint vectors depending on axis of interest (ALTERX)
  if(ALTERX) { E = XE; E.other = YE } else { E = YE; E.other = XE }
    
  ## check if there are at least one possible CP to shift 
  if(length(E) < 3) {
    if(DEBUG_BIRTH_EXT == TRUE) { cat("nothing to shift, returning..\n") }
    return(list(XE=XE, YE=YE, S2Dall=S2Dall, B2Dall=B2Dall, Sig2_2Dall=Sig2_2Dall, accept=0, move=2, alpha=0, estar=-1))
  }
  
  ## extract the CP to be shifted
  estar = sample(c(E[2:(length(E)-1)], E[2:(length(E)-1)]), 1)

  ## Position of the phase starting at the selected CP
  poskstar = sum(E <= estar)

  ## Possible new position for the selected CP (CP-1, CP+1)
  newCPs = c(E[poskstar]-1,E[poskstar]+1)
  
  ## remove positions that create too short segments (minPhase)
  newCPs = newCPs[which( !( newCPs %in% c(E[poskstar-1],E[poskstar-1]+minPhase-1, E[poskstar+1],E[poskstar+1]-minPhase+1,E)))]
  if(DEBUG_BIRTH_EXT == TRUE) { cat("newCPs: ", newCPs, "\n") }
  
  ## Boolean for the acceptation of the CP shift move initially set to 0 (=1 if birth accepted, 0 otherwise)
  accept = 0
  
  ## If there is at least one option to shift the selected CP 
  if(length(newCPs) > 0){

    ## sample new CP position and check if shift produces valid segments
    cp.new = sampleValidateCPs(candidateCPs=newCPs, min.seglocs = minPhase*minPhase, E, E.other, Y, ALTERX, xlocs,  type="shift", cp.pos = poskstar) 
    
    ## if no valid shift could be found, return
    if(is.nan(cp.new)) {
      #cat("giving up new cp - shift\n") 
      return(list(XE=XE, YE=YE, S2Dall=S2Dall, B2Dall=B2Dall, Sig2_2Dall=Sig2_2Dall, accept=0, move=2, alpha=0, estar=-1))
    }

    ## new CP vector
    Eshift = E
    Eshift[poskstar] = cp.new
        
    ###
    ## Calculate the current state posterior, unshifted
  #  prodPhi = 1     # product over current state
    sumPhi = 0
  
    ## Assign proper CP vector to temporary type (used in segcoord extraction)
    if(ALTERX) {tmpYE = YE; tmpXE = E } else { tmpXE = XE; tmpYE = E }
        
    for(j in (poskstar-1):poskstar) {
      
      for(i in 1:(length(E.other)-1)) {

        ## Assign proper segment id
        if(ALTERX) { xsegid = j; ysegid = i } else { xsegid = i; ysegid = j}
      
        ## get segment coordinates
        segcoord = c(XMphase[tmpXE[xsegid]], YMphase[tmpYE[ysegid]],XMphase[tmpXE[xsegid+1]]-1, YMphase[tmpYE[ysegid+1]]-1)

        if(DEBUG_BIRTH_EXT == TRUE) {  cat("[prodPhi] xsegid: ", xsegid, ", ysegid: ", ysegid, ",  segcoord ")
                                       print.table(segcoord) }
        ## get the predictor data
        x = extractNodes(X, segcoord, xlocs,F)
        
        ## get the target data
        y = extractNodes(Y, segcoord, xlocs,F)
        
        ## number of locations
        omega = length(y)
        
        ## calculate projection matrix
        Pr = computeProjection(as.matrix(x[,which(S2Dall == 1)]), delta2)

        if( (dim(Pr)[1] != length(y)) ) {
          #browser()
        }
        
#        prodPhi = prodPhi * gamma((v0+omega)/2) * ((gamma0+ t(y) %*% Pr %*% y)/2)^(-(v0+omega)/2)

        tryCatch({
          
          sumPhi  = sumPhi  + lgamma((v0+omega)/2) + (-(v0+omega)/2) * log( (gamma0+ t(y) %*% Pr %*% y)/2)

        }, error = function(e) {
          cat("Caught error \n ")
          print(e)
          #browser()
        })



      }
    }

    ##
    ## Calculate the next state posterior, shifted
    ## NOTE, this is exactly the same as with current state, only the Eshift assigned as tmpXE or tmpYE
 #   prodPhiPlus = 1     # product over current state
    sumPhiPlus = 0
  
    ## Assign proper CP vector to temporary type (used in segcoord extraction)
    if(ALTERX) {tmpXE = Eshift } else { tmpYE = Eshift }
    
    for(j in (poskstar-1):poskstar) {
      
      for(i in 1:(length(E.other)-1)) {

        ## Assign proper segment id
        if(ALTERX) { xsegid = j; ysegid = i } else { xsegid = i; ysegid = j}
      
        ## get segment coordinates
        segcoord = c(XMphase[tmpXE[xsegid]], YMphase[tmpYE[ysegid]],XMphase[tmpXE[xsegid+1]]-1, YMphase[tmpYE[ysegid+1]]-1)

        if(DEBUG_BIRTH_EXT == TRUE) {  cat("[prodPhi] xsegid: ", xsegid, ", ysegid: ", ysegid, ",  segcoord ")
                                       print.table(segcoord) }
        ## get the predictor data
        x = extractNodes(X, segcoord, xlocs,F)
      
        ## get the target data
        y = extractNodes(Y, segcoord, xlocs,F)
      
        ## number of locations
        omega = length(y)
      
        ## calculate projection matrix
        Pr = computeProjection(as.matrix(x[,which(S2Dall == 1)]), delta2)

#        prodPhiPlus = prodPhiPlus * gamma((v0+omega)/2) * ((gamma0+ t(y) %*% Pr %*% y)/2)^(-(v0+omega)/2)

        if( (dim(Pr)[1] != length(y)) || (dim(Pr)[2] != length(y))) {
          #browser()
        }
        
        tryCatch({
          sumPhiPlus  = sumPhiPlus  + lgamma((v0+omega)/2) + (-(v0+omega)/2) * log( (gamma0+ t(y) %*% Pr %*% y)/2)
        }, error = function(e) {
          cat("Caught error \n ")
          print(e)
          #browser()
        })

      }
    }

    ## Computation of the proposal Ratio:
    ##
    ## This operation corresponds to (kW - e) / (kW - e*) at the end of Eq. 4.25, where 'e' is the number of impossible position changes
    ## e  = sum(((E[2:(s+2)]-E[1:(s+1)]) <= minPhase) * nbmove)
    ## e* = sum(((Estar[2:(s+2)]-Estar[1:(s+1)]) <= minPhase) * nbmove)
    ## W = 2
    ## k = nr. of segments
    ##
    ## In consequence: (2*s-sum(((E[2:(s+2)]-E[1:(s+1)]) <= minPhase) * nbmove)) = W*k - e -> nr. of possible shift
    ## ---> The ratio is larger than 1 if the shift decreases the nr. of possible shifts and smaller than one if the shift introduces more shift opportunities
    
    ## Vector of length the current number of phases= c(1,2,2,...,2,2,1) i.e. the number of CP that can potentially be shifted into each phase
    k = length(E) - 2
    nbmove = c(1,array(2,k-1),1)
    propRatio = (2*k-sum(((E[2:(k+2)]-E[1:(k+1)]) <= minPhase) * nbmove))/(2 * k - sum(((Eshift[2:(k+2)]-Eshift[1:(k+1)]) <= minPhase) * nbmove))
    
    
    ## calculate acceptance probability
    ##alpha1 = prodPhiPlus / prodPhi 
    logLR = sumPhiPlus - sumPhi
    
    if(DEBUG_BIRTH_EXT == TRUE) {  cat("logLR: ", logLR, ", propRatio: ", propRatio, "\n");
                                   cat("logLR + log(propRatio): ", logLR+log(propRatio), ", exp(..): ", exp(logLR+log(propRatio)), "\n");
                                   cat("exp(logLR)*propRatio: ", exp(logLR)*propRatio, "\n")
                                 }

    
    ## Computation of alpha
    if(!is.nan(logLR) & (logLR+log(propRatio))<0){ # in case adding up the log yields < 0
      alpha = min(c(1, exp(logLR)*propRatio)) # do none log calculation
    } else { # else alpha would be higher than 1 , because exp(0) = 1, so we can set it anyways
      alpha = 1
    }


    ## Sample u to decide whether the CP shift is accepted or not
    u = runif(1,0,1)
    
    ## Boolean for the acceptation of the CP death move (=1 if birth accepted, 0 otherwise)

    if(u <= alpha){
      ## Acceptation of the death of the selected CP
      ## Move acceptation boolean =1
      accept = 1

      ## save new CP vector
      if(ALTERX) { XE = Eshift } else { YE = Eshift }
      
      if(DEBUGLVL1) { cat("s") }
      
    }
  } 
  
  ##  Return all variables
  ## (+ variable move describing the move type  (1= CP birth, 2= CP death, 3= CP shift, 4= Update phases)
  return(list(XE=XE, YE=YE, S2Dall=S2Dall, B2Dall=B2Dall, Sig2_2Dall=Sig2_2Dall, accept=accept, move=3))
}   # end CP.SHIFT




###################################################################
# Update phases
###################################################################

phase.update <- function(XE, YE, S2Dall, B2Dall, Sig2_2Dall, X, Y, GLOBvar, HYPERvar,  DEBUGLVL1 = F, DEBUGLVL2 = F ) {

  
  
  ### assignement of global variables used here ###
  q = GLOBvar$q
  smax = GLOBvar$smax
  
  ### assignement of hyperparameters variables used here ###
  c = HYPERvar$c
  alphalbd = HYPERvar$alphalbd
  betalbd = HYPERvar$betalbd
  alphad2 = HYPERvar$alphad2
  betad2 = HYPERvar$betad2
  v0 = HYPERvar$v0
  gamma0 = HYPERvar$gamma0
  ### end assignement ###

  if(DEBUGLVL2) { cat("\n -- START PHASE.UPDATE --\n\n") }

  ## current number of edges
  s = sum(S2Dall) - 2

  ## expected nr of edges (mean)
  Lambda = rgamma(1, shape=s + alphalbd, rate=1 + betalbd)
  
  ## Compute acceptation probability vector rho
  rho3 = computeRho3(s, 0, smax, c, Lambda)
  
  ## Sample u
  u = runif(1, 0, 1)
  
  ## Compute the corresponding move (Edge birth, Edge death or Update the regression coefficient) 
  bduout = bdu.homogeneousStructure(u, rho3, X, Y, XE, YE, S2Dall, Sig2_2Dall, q, v0, gamma0, smax, GLOBvar, HYPERvar, DEBUGLVL1, DEBUGLVL2)


  ## (+ variable move describing the move type  (1= CP birth, 2= CP death, 3= CP shift, 4= Update segments)
  return( list( XE=XE, YE=YE, S2Dall=bduout$S2Dall, B2Dall=bduout$B2Dall, Sig2_2Dall=bduout$Sig2_2Dall, move=4, accept=0, move1=bduout$move, accept1=bduout$accept, delta2=bduout$delta2))

}


bdu.homogeneousStructure <- function(u, rho3, X, Y, XE, YE, S2Dall, Sig2_2Dall, q, v0, gamma0, smax, GLOBvar, HYPERvar, DEBUGLVL1 = F, DEBUGLVL2 = F){

  ### INPUT:u,rho,s=S[i,],sig2=Sig2[i],delta2.
  ###	x: the data in state i in columns 
  ###	ni: total nb of repeated measurements	
  ### OUTPUT: newS,newSig2,newB.
  ### depends on: 
  ### q the number of predictors
  ### constant v0, gamma0.
    
  XMphase = GLOBvar$XMphase
  YMphase = GLOBvar$YMphase
  xlocs = GLOBvar$xlocs
  ylocs = GLOBvar$ylocs
  delta2 = HYPERvar$delta2
	
  ## Variable move describing the move type  (1= Edge birth, 2= Edge death, 3= Update coefficient, default=3)
  move = 3

  ## Boolean indicating whether the move is accepted or not (=1 if accepted, 0 otherwise, default=0)
  accept = 0

  ## New edges vector, to be returned at the end of the function
  newS = S2Dall
  
  ## Current number of edges
  s = sum(S2Dall) - 2 
  
  ## Choose between flip move and other moves
  choice = runif(1, 0, 1)
  
  ## Flip Move
  if(s > 0 && choice > 0.75) {    
  
    ## flip move is 4
    move = 4

    tryCatch({
    
      ## Sample the original parent
      parent.orig = sample(c(which(S2Dall[1:q]==1), which(S2Dall[1:q]==1)),1) # needed when there is only one position  S[1:q]==1
      
    }, error = function(e) {
      cat("FIXME: Caught error in bdu.homogeneousStructure while doing parent.orig = sample(c..)\n ")
      print(e)
      stop("stopping")
    })
    
                                        # Sample the new parent
    parent.new = sample(c(which(S2Dall==0), which(S2Dall==0)), 1) # needed when there is only one position  S==0
        
    stmp = S2Dall
    stmp[parent.orig] = 0
    stmp[parent.new]  = 1

    rfliplog = 0    

    ## loop over each current segment
    for(xsegid in 1:(length(XE)-1)) {

        for(ysegid in 1:(length(YE)-1)) {
    
          ## get segment coordinates
          segcoord = c(XMphase[XE[xsegid]], YMphase[YE[ysegid]],XMphase[XE[xsegid+1]]-1, YMphase[YE[ysegid+1]]-1)
          
          if(DEBUGLVL2) {  cat("[prodPhiStar-flip] xsegid: ", xsegid, ", ysegid: ", ysegid, ",  segcoord ")
                           print.table(segcoord) }
          
          ## get the predictor and target data
          x = extractNodes(X, segcoord, xlocs,F)
          y = extractNodes(Y, segcoord, xlocs,F)
          
	  ## number of locations
          omega = length(y)
          
          ## calculate projection matrix
          Pr = computeProjection(as.matrix(x[,which(S2Dall == 1)]), delta2)       # current structure 
          Prstar = computeProjection(as.matrix(x[,which(stmp == 1)]), delta2) # proposed structure

#	  rflip = rflip * ((gamma0 + t(y) %*% Prstar %*% y)/(gamma0 + t(y) %*% Pr %*% y))^(-(length(y) + v0)/2)

          ## this try catch is only for debugging
          tryCatch({

            ldebug = log(gamma0 + t(y) %*% Prstar %*% y)
	
          }, warning = function(e) {
            cat("Caught warning while doing log(gamma0 + t(y) %*% Prstar %*% y): ")
            print(e)

            cat("gamma0: ") 
            print.table(gamma0)
            cat("saving all parameters to file debugsave: ")
            save(gamma0, y, Prstar, file="debugsave")
            
            stop("DEBUG-E2: failed because of NaN , see output")    
            
          })

	# log version
          
          tryCatch({
             rfliplog = rfliplog + (-(length(y) + v0)/2) * ( log(gamma0 + t(y) %*% Prstar %*% y) - log(gamma0 + t(y) %*% Pr %*% y))
           }, error = function(e) {
             cat("Caught error \n ")
             print(e)
             #browser()
           })
          
	}
    }   

    u = runif(1,0,1)
    
    if(u <= min(1,exp(rfliplog))) {
      accept = 1
      newS = stmp

      if(DEBUGLVL1) { cat("f") }

    }
    
  } else {

    ##########################
    ## Birth of an edge move
    ##########################
    if(u < rho3[1] && s < smax && (length(which(S2Dall == 0)) > 0 ) ){

      
      ## Variable move describing the move type  (1= Edge birth, 2= Edge death, 3= Update coefficient)
      move = 1
      
      tryCatch({

        ## Sample the additional edge
        sstar = sample(c(which(S2Dall==0), which(S2Dall==0)), 1) # needed when there is only one position  S2Dall==0
      }, error = function(e) {

      	cat("FIXME: Caught error in bdu.homogeneousStructure() - birth edge move, sample() for sstar:\n")
       	print(e)
	stop("stopping")
      })
      

      ## Proposed edges vector (with an additional edge)
      stmp = S2Dall
      stmp[sstar] = 1
      
      ## product over posterior probs.       
      rbirth = 1
      
      ## loop over each current segment
      for(xsegid in 1:(length(XE)-1)) {

        for(ysegid in 1:(length(YE)-1)) {

          ## get segment coordinates
          segcoord = c(XMphase[XE[xsegid]], YMphase[YE[ysegid]],XMphase[XE[xsegid+1]]-1, YMphase[YE[ysegid+1]]-1)
          
          if(DEBUGLVL2) {  cat("[prodPhiPlus] xsegid: ", xsegid, ", ysegid: ", ysegid, ",  segcoord ")
                           print.table(segcoord) }
          
          ## get the predictor and target data
          x = extractNodes(X, segcoord, xlocs,F)
          y = extractNodes(Y, segcoord, xlocs,F)
          
          ## calculate projection matrix
          Pr =     computeProjection(as.matrix(x[,which(S2Dall == 1)]), delta2)       # current structure 
          Prplus = computeProjection(as.matrix(x[,which(stmp == 1)]), delta2)     # + 1 edge

          ## number of locations
          omega = length(y)

          ## Compute birth ratio, no log needed because no critical gamma() calculation
          ## orig 1D homog.:
          ##rbirth =    rbirth * ((gamma0 + t(y) %*% Pxlp1 %*% y) /(gamma0 + t(y) %*% Pxl %*% y))^(-(length(y) + v0)/2)/sqrt(1 + delta2)

          tryCatch({
            rbirth =rbirth * ((gamma0 + t(y) %*% Prplus %*% y) / (gamma0 + t(y) %*% Pr %*% y))^(-(v0 + omega)/2)/sqrt(1 + delta2)
          }, error = function(e) {
            cat("Caught error \n ")
            print(e)
            #browser()
          })
        
        }
      }

      ## Sample u 
      u = runif(1,0,1)
        
      if(u <= min(1,rbirth)){
        accept = 1
        newS = stmp
        if(DEBUGLVL1) { cat("e") }
      }
            	  
     } else {
       
      if(u < rho3[2] & s > 0){  # makes sure at least one edge exists (despite the bias and SAC edge)

        ##########################
        ## Death of an edge move
        ##########################

        ## Variable describing the move type  (1 for Edge birth, 2 for Edge death, 3 for Update coefficient)
        move=2

        tryCatch({

          ## Sample the edge to remove
          sstar = sample(c(which(S2Dall[1:q]==1), which(S2Dall[1:q]==1)),1) # needed when there is only one position  S[1:q]==1

	}, error = function(e){
          cat("FIXME: Caught error in bdu.homogeneousStructure() - edge death move while doing sample() for sstar:\n")
          print(e)
          stop("stopping")
        })

        ## Proposed edges vector (after taking away one edge)
        stmp = S2Dall
        stmp[sstar] = 0
        
        rdeath =1
                        
        ## loop over each current segment
        for(xsegid in 1:(length(XE)-1)) {

          for(ysegid in 1:(length(YE)-1)) {

            ## get segment coordinates
            segcoord = c(XMphase[XE[xsegid]], YMphase[YE[ysegid]],XMphase[XE[xsegid+1]]-1, YMphase[YE[ysegid+1]]-1)
          
            if(DEBUGLVL2) {  cat("[prodPhiPlus] xsegid: ", xsegid, ", ysegid: ", ysegid, ",  segcoord ")
                             print.table(segcoord) }
          
            ## get the predictor and target data
            x = extractNodes(X, segcoord, xlocs,F)
            y = extractNodes(Y, segcoord, xlocs,F)
          
            ## calculate projection matrix
            Pr =      computeProjection(as.matrix(x[,which(S2Dall == 1)]), delta2)       # current structure 
            Prminus = computeProjection(as.matrix(x[,which(stmp == 1)]), delta2)   # + 1 edge
          
            ## complies to TVDBN_SH1D
            tryCatch({
              rdeath = rdeath * ( (gamma0 + t(y) %*% Pr %*% y) / (gamma0 + t(y) %*% Prminus %*% y))^((length(y) + v0)/2)*(sqrt(1 + delta2))
            }, error = function(e) {
              cat("Caught error \n ")
              print(e)
              #browser()
            })
        
          }
        }

        ## Sample u 
        u<-runif(1,0,1)
	  
        if(u <= min(1,rdeath)){
          ## Boolean for the acceptation of the CP death move (=1 if birth accepted, 0 otherwise)
          accept = 1
          newS = stmp
          if(DEBUGLVL1) { cat("n") }
        }
      } 
    } 
  }
 
    
  ##
  ## Updating coefficients of each segment (regardless if structure changed)
  ##

  ## create regression parameters from scratch, this will replace the original one fully
  B2Dall = matrix(0, 0, q+4) # +4 because xcoord,ycoord,bias regr.,SAC regr.

  # update the variance sigma globally
  Sig2_2Dall = updateSigGlobal(xlocs, XMphase, YMphase, XE, YE, X,Y, S2Dall, delta2, HYPERvar$v0, HYPERvar$gamma0)
  
  ## loop over each current segment
  for(xsegid in 1:(length(XE)-1)) {

    for(ysegid in 1:(length(YE)-1)) {
      
      ## get segment coordinates
      segcoord = c(XMphase[XE[xsegid]], YMphase[YE[ysegid]],XMphase[XE[xsegid+1]]-1, YMphase[YE[ysegid+1]]-1)
      
      if(DEBUGLVL2) {  cat("[prodPhiPlus] xsegid: ", xsegid, ", ysegid: ", ysegid, ",  segcoord ")
                       print.table(segcoord) }
      
      ## get the predictor and target data
      x = extractNodes(X, segcoord, xlocs,F)
      y = extractNodes(Y, segcoord, xlocs,F)
                              
      ## sample edge weights
      ## AA22.02.2011, +2 because of additional bias and SAC edge
      newB = array(0, q+2)
      newB[which(newS == 1)] = sampleBxy(x[, which(newS == 1)], y, Sig2_2Dall, delta2)
         
      ## set B
      B2Dall = rbind(B2Dall, c(xsegid, ysegid, newB))
      
    }
  }
      
  
  ##  Return all variables
  return(list( u=u, move=move, accept=accept, S2Dall=newS, B2Dall=B2Dall, Sig2_2Dall=Sig2_2Dall)) 
}












#####################################################################################
## BIRTH OF A CHANGEPOINT
#####################################################################################

cp.birth <- function(E, Sall, Ball, Sig2all, X, Y, D, GLOBvar, HYPERvar, post_probs){
  # INPUT: E, Sall, Ball ,Sig2all, X, Y, D, GLOBvar, HYPPERvar
  # OUTPUT: 
  # depends on: .

  # current number of changepoints
  s = length(E) - 2
  
  ### assignement of global variables used here ###
  q = GLOBvar$q
  qmax = GLOBvar$qmax
  Mphase = GLOBvar$Mphase
  minPhase = GLOBvar$minPhase
  nbVarMax = GLOBvar$nbVarMax
  smax = GLOBvar$smax
  dyn = GLOBvar$dyn
  birth_proposal = GLOBvar$birth_proposals
  qmax = GLOBvar$qmax
  ### end assignement ###

  ### assignement of hyperparameters variables used here ###
  alphalbd = HYPERvar$alphalbd
  betalbd = HYPERvar$betalbd
  alphad2 = HYPERvar$alphad2
  betad2 = HYPERvar$betad2
  v0 = HYPERvar$v0
  gamma0 = HYPERvar$gamma0
  ### end assignement ###
  
  # When using the Hamming distance - draw distance, then choose changed edges
  # randomly
  ham_dist = 0;
  
  # Proposal 5 is a mixture of 2 and 3, so choose one of them. For the moment
  # assume equal probability.
  if(birth_proposal == 5) {
	choice = runif(1,0,1) < 1/2;
	
	birth_proposal = 2*choice + 3*(1-choice);
  }
  
  
  if(birth_proposal == 3) {
	  ham_dist = rpois(1, 0.5);
  }

  ## search for possible CP, not in E and not close to E if minPhase (length of phase) is > than 1
  toremove = E
  if(minPhase>1) for(i in 1:(minPhase-1)) toremove = c(toremove, E-i, E+i)
  # possible CPs are those not in 'toremove'
  possibleCP = setdiff((1+dyn):E[length(E)], toremove)
  #possibleCP = which(!((1+dyn):E[length(E)] %in% c(E,E+minPhase-1,E-minPhase+1)))
  ## Sample the new CP "estar"
  estar = sample(c(possibleCP, possibleCP),1)
  
  ## Position of the phase containing the new CP
  poskl = sum(E < estar)
  
  ## Current edges vector S in the phase containing the new CP
  Sold = Sall[poskl,]
  
  ## Current number of edges k in the phase containing the new CP
  k = sum(Sold) - 1

  ## Sample lambda
  lambda = rgamma(1, shape=alphalbd, rate=betalbd)

  ## Sample a new edges vector newS
    ## updated by Sophie 04/08/10:  update for the case homogeneousStructure=TRUE
  if(GLOBvar$homogeneousStructure){
    newS = Sold
	newRight = 0
    sL = Sold
    sR = Sold
  }else{
    newS = array(1, q+1)
	prop_prob = 1;
	
  if(birth_proposal == 1 || birth_proposal == 2) {
    
    newS[1:q] = 1:q %in% sample(1:q, sampleK(0, qmax, lambda, 1), replace=FALSE)
  } else if(birth_proposal == 3) {
  	changed = 1:q %in% sample(1:q, ham_dist, replace=FALSE)
	tempS = Sold[1:q];
	tempS[changed] = !tempS[changed];
	newS[1:q] = tempS
  } else if(birth_proposal == 4) {
	u = runif(1, 0, 1)
	cumul = post_probs$cumul;
	cum_prob = which(cumul[,dim(cumul)[2]] >= u);
    sampled_edges = cumul[cum_prob[1], 1:qmax]
    newS[1:q] = 1:q %in% sampled_edges;
		
    index = paste(sampled_edges, collapse=' ');
	probs = post_probs$probs;
	prop_prob = probs[[index]];
  }
 	
  ## Sample Right or Left (the position for the new edges vector: to the right or to the left of the new CP estar)
  if(runif(1,0,1) < 1/2){
    ##  New edges vector to the left of the new CP
    ## Boolean (= 1 if  the new edges vector is to the right of the new CP, 0 otherwise) 
    newRight = 0
    sL = newS
    sR = Sold
  } else {
    ## New edges vector to the right of the new CP
    ## Boolean (= 1 if  the new edges vector is to the right of the new CP, 0 otherwise)
    newRight = 1
    sR = newS
    sL = Sold
  }

 } # end if(GLOBvar$homogeneousStructure)

  
  ## Compute the matrices required for the computation of the acceptation probability alpha
  yL = Y[(Mphase[E[poskl]]:(Mphase[estar]-1))]
  xL = X[(Mphase[E[poskl]]:(Mphase[estar]-1)),]
  yR = Y[(Mphase[estar]:(Mphase[E[poskl+1]]-1))]
  xR = X[(Mphase[estar]:(Mphase[E[poskl+1]]-1)),]
  y2 = array(c(yL, yR))
  x2 = rbind(xL, xR)
  
  ## Updating parameters
  
  if(nbVarMax > 1){
    Sig2 = Sig2all[poskl]
  } else {
    Sig2 = Sig2all
  }

  delta2 = sampleDelta2(poskl, x2, q, Ball, Sall, Sig2, alphad2, betad2)
  
  ## Compute projection of the matrices required for the computation of the acceptation probability alpha
  PxL = computePx(length(yL), as.matrix(xL[,which(sL == 1)]), delta2)
  PxR = computePx(length(yR), as.matrix(xR[,which(sR == 1)]), delta2)
  Px2 = computePx(length(y2), as.matrix(x2[,which(Sold == 1)]), delta2)
  
  
  
	
  ## Compute the acceptation probability alpha
	# birth_proposal = 1: Old acceptance probability
	#                 = 2: Corrected acceptance probability
	#                 = 3: Hamming Distance-based acceptance probability
	#                 = 4: Changepoint-less posterior network probability-based 
	#                      acceptance probability
	alpha = switch(birth_proposal, 
                 bp.computeAlpha(1, sum(newS)-1, s, Mphase[E[poskl]]
				 , Mphase[estar], Mphase[E[poskl+1]], yL, PxL, yR, PxR, y2, Px2, D, delta2, q, smax, v0, gamma0),
				 bp.computeAlpha_updated(1, sum(newS)-1, s, Mphase[E[poskl]]
				 , Mphase[estar], Mphase[E[poskl+1]], yL, PxL, yR, PxR, y2, Px2, D, delta2, q, smax, v0, gamma0),
				 bp.computeHammingAlpha(1, sum(newS)-1, s, Mphase[E[poskl]]
				 , Mphase[estar], Mphase[E[poskl+1]], yL, PxL, yR, PxR, y2, Px2, D, delta2, q, smax, v0, gamma0, ham_dist),
				 bp.computePrecompAlpha(1, sum(newS)-1, s, Mphase[E[poskl]]
				 , Mphase[estar], Mphase[E[poskl+1]], yL, PxL, yR, PxR, y2, Px2, D, delta2, q, smax, v0, gamma0, prop_prob))
				 
  #alpha
  
  ## Sample u to conclude either to  acceptation or to rejection
  u = runif(1,0,1)
  
  ## Boolean for the acceptation of the CP birth move initially set to 0 (=1 if birth accepted, 0 otherwise)
  accept = 0
  
  if(!is.nan(alpha) & u <= alpha){
    ## Acceptation of the birth of the new CP
    ## Move acceptation boolean =1
    accept=1
    
    ## Compute new Sig2
    if(nbVarMax>1){
      newSig2 = array(0, s+2)
      newSig2[(1:(s+2))[-c(poskl,poskl+1)]] = Sig2all[(1:(s+1))[-c(poskl)]]
    }
    ## Compute new regression parameters newB
    newB = matrix(0, s+2, q+1)
    newB[(1:(s+2))[-c(poskl,poskl+1)],] = Ball[(1:(s+1))[-c(poskl)],]
    
    ## Update newB 
    if(newRight == 0){
      ## Update the phase to the left of the new CP (in newB)
      if(nbVarMax>1){
        newSig2[poskl] = sampleSig2(yL,PxL,v0,gamma0)
        newSig2[poskl+1] = Sig2all[poskl]
        Sig2all = newSig2
        Sig2 = newSig2[poskl]
      }
      
      newB[poskl+1,] = Ball[poskl,]
      newB[poskl, which(newS == 1)] = sampleBxy(xL[,which(newS==1)], yL, Sig2, delta2)
    } else {
      ## Update the phase to the right of the new CP (in newB)
      if(nbVarMax>1){
        newSig2[poskl+1] = sampleSig2(yR, PxR, v0, gamma0)
        newSig2[poskl] = Sig2all[poskl]
        Sig2all = newSig2
        Sig2 = newSig2[poskl+1]
      }
      
      newB[poskl,] = Ball[poskl,]
      newB[poskl+1, which(newS == 1)] = sampleBxy(xR[, which(newS == 1)], yR, Sig2, delta2)
    }

    ##if(nbVarMax==1){
    ##  Sig2all=updateSigSolo(X,Y,E,Sall,Ball,Sig2all)
    ##}
    
    ## Update current model and parameters
    Ball = newB
    Sall = (abs(Ball)>0)*1
    E = sort(c(E,estar))
  }

  ##  Return all variables
  ## (+ variable move describing the move type  (1= CP birth, 2= CP death, 3= CP shift, 4= Update phases)
  return(list(E=E, Sall=Sall, Ball=Ball, Sig2all=Sig2all, accept=accept, move=1, alpha=alpha, estar=estar, lambda=lambda))
}


#####################################################################################
#DEATH OF A CHANGEPOINT
#####################################################################################

cp.death <- function(E, Sall, Ball, Sig2all, X, Y, D, GLOBvar, HYPERvar, post_probs){
  ### INPUT:  E, Sall, Ball, Sig2all, X, Y, D, GLOBvar, HYPERvar
  ### OUTPUT: 
  ### depends on: .

  # current number of changepoints
  s = length(E) - 2
  
  ### assignement of global variables used here ###
  q = GLOBvar$q
  Mphase = GLOBvar$Mphase
  nbVarMax = GLOBvar$nbVarMax
  smax = GLOBvar$smax
  birth_proposal = GLOBvar$birth_proposals
  qmax = GLOBvar$qmax
  ### end assignement ###

  ### assignement of hyperparameters variables used here ###
  alphad2 = HYPERvar$alphad2
  betad2 = HYPERvar$betad2
  v0 = HYPERvar$v0
  gamma0 = HYPERvar$gamma0
  ### end assignement ###
  
  ## Sample the CP to be removed
  estar = sample(c(E[2:(length(E)-1)], E[2:(length(E)-1)]), 1)
  ##  Position of the phase starting at the selected CP
  poskstar = sum(E <= estar)

  # Proposal 5 is a mixture of 2 and 3, so choose one of them. For the moment
  # assume equal probability.
  if(birth_proposal == 5) {
	choice = runif(1,0,1) < 1/2;
	
	birth_proposal = 2*choice + 3*(1-choice);
  }
  
  # Placeholder variables (only used if correct birth proposal method in use)
  ham_dist = -1;
  prop_prob = -1;
  
  ## Sample the edge vector to be maintained  for merge phase after removing the CP
  ## (newRight =1 if the edge vector to the Right of the CP is maintained, =0 otherwise)
  newRight = sample(0:1,1)
  
  ## Position of the line of the phase to be removed in the matrices Sall and Ball
  away = (1-newRight) * poskstar + newRight * (poskstar-1)
  
  if(birth_proposal == 3) {
	  before = poskstar - 1;
	  
	  B_1 = Ball[poskstar,]
    S_1 = (abs(B_1) > 0) * 1
    
    B_2 = Ball[before,]
    S_2 = (abs(B_2) > 0) * 1

  	ham_dist = sum(abs(S_2[1:q] - S_1[1:q]));
  } else if(birth_proposal == 4) {
  	B_away = Ball[away,];
  	S_away = (abs(B_away) > 0) * 1;
  	
    edges = which(S_away[1:q] == 1)
  	index = paste(edges, collapse=' ')
  	
  	missing = qmax - length(strsplit(index, ' ')[[1]]);
  	zeros = "";
  	
  	if(missing > 0) {
  	  zeros = c();
  		for (i in 1:missing) {
  		  zeros = c(zeros, 0);
  	  }
  	  
  	  if(nchar(index) == 0) {
  	    index = paste(zeros, collapse=' ');
  	  } else {
  	  	index = paste(paste(zeros, collapse=' '), index, sep=' ');
  	  }
  	}
  	
  	probs = post_probs$probs;

  	prop_prob = probs[[index]];

  }
  
  ## Compute the matrices required for the computation of the acceptation probability alpha
  yL = Y[(Mphase[E[poskstar-1]]:(Mphase[estar]-1))]
  xL = X[(Mphase[E[poskstar-1]]:(Mphase[estar]-1)),]
  yR = Y[(Mphase[estar]:(Mphase[E[poskstar+1]]-1))]
  xR = X[(Mphase[estar]:(Mphase[E[poskstar+1]]-1)),]
  y2 = array(c(yL,yR))
  x2 = rbind(xL,xR)
 
  if(nbVarMax>1){
    Sig2 = Sig2all[poskstar - 1 + newRight]
  } else {
    Sig2 = Sig2all
  }

  ## Update delta
  delta2 = sampleDelta2(poskstar-1+newRight, x2, q, Ball, Sall, Sig2, alphad2, betad2)

  ## Compute projection of the matrices required for the computation of the acceptation probability alpha
  #modif 17 avril 
  Px2 = computePx(length(y2), as.matrix(x2[,which(Sall[poskstar-1+newRight,] == 1)]), delta2) 
  PxL = computePx(length(yL), as.matrix(xL[,which(Sall[poskstar-1,] == 1)]), delta2)
  PxR = computePx(length(yR), as.matrix(xR[,which(Sall[poskstar,] == 1)]), delta2)

  ## Compute the acceptation probability alpha
  alpha = switch(birth_proposal, 
                 bp.computeAlpha(-1, sum(Sall[poskstar-1+newRight,])-1, s-1, 
                   Mphase[E[poskstar-1]], Mphase[estar], Mphase[E[poskstar+1]], 
                   yL, PxL, yR, PxR, y2, Px2, D, delta2, q, smax, v0, gamma0),
				 bp.computeAlpha_updated(-1, sum(Sall[poskstar-1+newRight,])-1, s-1, 
                   Mphase[E[poskstar-1]], Mphase[estar], Mphase[E[poskstar+1]], 
                   yL, PxL, yR, PxR, y2, Px2, D, delta2, q, smax, v0, gamma0),
				 bp.computeHammingAlpha(-1, sum(Sall[poskstar-1+newRight,])-1, s-1, 
                   Mphase[E[poskstar-1]], Mphase[estar], Mphase[E[poskstar+1]], 
                   yL, PxL, yR, PxR, y2, Px2, D, delta2, q, smax, v0, gamma0,
                   ham_dist),
				 bp.computePrecompAlpha(-1, sum(Sall[poskstar-1+newRight,])-1, s-1, 
                   Mphase[E[poskstar-1]], Mphase[estar], Mphase[E[poskstar+1]], 
                   yL, PxL, yR, PxR, y2, Px2, D, delta2, q, smax, v0, gamma0,
                   prop_prob))
  
  ## Sample u to conclude either to  acceptation or to rejection
  u = runif(1,0,1)

  ## Boolean for the acceptation of the CP death move initially set to 0 (=1 if birth accepted, 0 otherwise)
  accept = 0
  
  if(!is.nan(alpha) & u <= alpha){
    ## Acceptation of the death of the selected CP
    ## Move acceptation boolean =1
    accept=1
   
    ## Remove the CP in E and the phase in the matrices Sall and Ball
    if(nbVarMax>1){
      Sig2all = Sig2all[(1:(s+1))[-c(away)]]
    }
    E = E[(1:(s+2))[-c(poskstar)]]
    Sall = matrix(Sall[(1:(s+1))[-c(away)],], s, q+1)
    Ball = matrix(Ball[(1:(s+1))[-c(away)],], s, q+1)
  }

  ##  Return all variables
  ## (+ variable move describing the move type  (1= CP birth, 2= CP death, 3= CP shift, 4= Update phases)
  return( list( E=E, Sall=Sall, Ball=Ball, Sig2all=Sig2all, accept=accept, move=2, alpha=alpha, estar=estar))
}



#####################################################################################
#MOVE OF A CHANGEPOINT
#####################################################################################

cp.shift <- function(E, Sall, Ball, Sig2all, X, Y, GLOBvar, HYPERvar){

  ### INPUT: u,rho,E,S,B,Sig2,X,Y
  ### OUTPUT: 
  ### depends on: .

  # current number of changepoints
  s = length(E) - 2
  
  ### assignement of global variables used here ###
  q = GLOBvar$q
  Mphase = GLOBvar$Mphase
  minPhase = GLOBvar$minPhase
  nbVarMax = GLOBvar$nbVarMax
  smax = GLOBvar$smax
  ### end assignement ###

  ### assignement of hyperparameters variables used here ###
  alphad2 = HYPERvar$alphad2
  betad2 = HYPERvar$betad2
  v0 = HYPERvar$v0
  gamma0 = HYPERvar$gamma0
  ### end assignement ###
  
  ## Sample the CP to be shifted
  estar = sample(c(E[2:(length(E)-1)], E[2:(length(E)-1)]), 1)

  ## Position of the phase starting at the selected CP
  poskstar = sum(E <= estar)

  ## Possible new position for the selected CP (CP-1, CP+1)
  newCPs = c(E[poskstar]-1,E[poskstar]+1)
  # remove positions that create too short phases (minPhase)
  newCPs = newCPs[which( !( newCPs %in% c(E[poskstar-1],E[poskstar-1]+minPhase-1, E[poskstar+1],E[poskstar+1]-minPhase+1,E)))]
  
  ## Boolean for the acceptation of the CP shift move initially set to 0 (=1 if birth accepted, 0 otherwise)
  accept = 0
  
  ## If there is at least one option to shift the selected CP 
  if(length(newCPs) > 0){
    ## Sample new CP position
    newCP = sample( c(newCPs, newCPs), 1)
    
    ## Compute the matrices required for the computation of the acceptation probability alpha
    ## Matrices for the model before shifting the CP
    yL = Y[ Mphase[E[poskstar-1]]:(Mphase[estar]-1) ]
    xL = X[ Mphase[E[poskstar-1]]:(Mphase[estar]-1), ]
    yR = Y[ Mphase[estar]:(Mphase[E[poskstar+1]]-1) ]
    xR = X[ Mphase[estar]:(Mphase[E[poskstar+1]]-1), ]

	
    ## Matrices for the model after shifting the CP
    yLStar = Y[ Mphase[E[poskstar-1]]:(Mphase[newCP]-1) ]
    xLStar = X[ Mphase[E[poskstar-1]]:(Mphase[newCP]-1), ]
    yRStar = Y[ Mphase[newCP]:(Mphase[E[poskstar+1]]-1) ]
    xRStar = X[ Mphase[newCP]:(Mphase[E[poskstar+1]]-1), ]
    
    if(nbVarMax>1){
      ## Update delta for the left phase
      delta2 = sampleDelta2(poskstar-1, xL, q, Ball, Sall, Sig2all[poskstar-1], alphad2, betad2)
      ## Compute projection of the matrices for the left phase
      PxL = computePx(length(yL), as.matrix(xL[, which(Sall[poskstar-1,] == 1)]), delta2)
      PxLStar = computePx(length(yLStar), as.matrix(xLStar[, which(Sall[poskstar-1,] == 1)]), delta2)
      ## Update delta for the right phase     
      delta2 = sampleDelta2(poskstar, xR, q, Ball, Sall, Sig2all[poskstar], alphad2, betad2)
      ## Compute projection of the matrices for the right phase
      PxR = computePx(length(yR), as.matrix(xR[, which(Sall[poskstar,] == 1)]), delta2)
      PxRStar = computePx(length(yRStar), as.matrix(xRStar[, which(Sall[poskstar,] == 1)]), delta2)
    } else {
      ## Update delta for the left phase
      delta2 = sampleDelta2(poskstar-1, xL, q, Ball, Sall, Sig2all, alphad2, betad2)
      ## Compute projection of the matrices for the left phase
      PxL = computePx(length(yL), as.matrix(xL[, which(Sall[poskstar-1,] == 1)]), delta2)
      PxLStar = computePx(length(yLStar), as.matrix(xLStar[, which(Sall[poskstar-1,] == 1)]), delta2)
      
      ## Update delta for the right phase     
      delta2 = sampleDelta2(poskstar, xR, q, Ball, Sall, Sig2all, alphad2, betad2)
      ## Compute projection of the matrices for the right phase
      PxR = computePx(length(yR), as.matrix(xR[, which(Sall[poskstar,] == 1)]), delta2)
      PxRStar = computePx(length(yRStar), as.matrix(xRStar[, which(Sall[poskstar,] == 1)]), delta2)
  
    }

    ## Compute the logarithm of the Likelihood Ratio (LR)
    ## logLR = log(gamma(((Mphase[newCP]-Mphase[E[poskstar-1]])+v0)/2)*gamma(((Mphase[E[poskstar+1]]-Mphase[newCP])+v0)/2))-log((gamma(((Mphase[E[poskstar]]-Mphase[E[poskstar-1]])+v0)/2)*gamma(((Mphase[E[poskstar+1]]-Mphase[E[poskstar]])+v0)/2)))+((Mphase[E[poskstar]]-Mphase[E[poskstar-1]])+v0)/2*log((gamma0+t(yL)%*%PxL%*%yL)/2)+((Mphase[E[poskstar+1]]-Mphase[E[poskstar]])+v0)/2*log((gamma0+t(yR)%*%PxR%*% yR)/2)-(((Mphase[newCP]-Mphase[E[poskstar-1]])+v0)/2)*log((gamma0+t(yLStar)%*%PxLStar %*%yLStar)/2)-((Mphase[E[poskstar+1]]-Mphase[newCP])+v0)/2*log((gamma0+t(yRStar) %*% PxRStar%*%yRStar)/2)
    
    ## Modifie par Sophie 01/03/08
    ## remplacement : log(gamma()) -> lgamma()
    logLR=lgamma(((Mphase[newCP]-Mphase[E[poskstar-1]])+v0)/2)+lgamma(((Mphase[E[poskstar+1]]-Mphase[newCP])+v0)/2)-lgamma(((Mphase[E[poskstar]]-Mphase[E[poskstar-1]])+v0)/2)- lgamma(((Mphase[E[poskstar+1]]-Mphase[E[poskstar]])+v0)/2)+ ((Mphase[E[poskstar]]-Mphase[E[poskstar-1]])+v0)/2* log((gamma0+t(yL)%*%PxL%*%yL)/2)+((Mphase[E[poskstar+1]]-Mphase[E[poskstar]])+v0)/2*log((gamma0+t(yR)%*%PxR%*% yR)/2)-(((Mphase[newCP]-Mphase[E[poskstar-1]])+v0)/2)*log((gamma0+t(yLStar)%*%PxLStar %*%yLStar)/2)-((Mphase[E[poskstar+1]]-Mphase[newCP])+v0)/2*log((gamma0+t(yRStar) %*% PxRStar%*%yRStar)/2)
  

    ## New CP vector Estar
    Estar = E
    Estar[poskstar] = newCP

    ## Computation of the proposal Ratio
    ## Vector of length the current number of phases= c(1,2,2,...,2,2,1) i.e. the number of CP that can potentially be shifted into each phase
    nbmove = c(1,array(2,s-1),1)
    propRatio = (2*s-sum(((E[2:(s+2)]-E[1:(s+1)]) <= minPhase) * nbmove))/(2 * s - sum(((Estar[2:(s+2)]-Estar[1:(s+1)]) <= minPhase) * nbmove))

    ## Computation of alpha
    if(!is.nan(logLR) & (logLR+log(propRatio))<0){
      alpha = min(c(1, exp(logLR)*propRatio))
    } else {
      alpha = 1
    }

    ## Sample u to decide whether the CP shift is accepted or not
    u = runif(1,0,1)

    ## Boolean for the acceptation of the CP death move (=1 if birth accepted, 0 otherwise)

    if(u <= alpha){
      ## Acceptation of the death of the selected CP
      ## Move acceptation boolean =1
      accept = 1
      E[poskstar] = newCP
    }
  }

  ##  Return all variables
  ## (+ variable move describing the move type  (1= CP birth, 2= CP death, 3= CP shift, 4= Update phases)
  return(list(E=E, Sall=Sall, Ball=Ball, Sig2all=Sig2all, accept=accept, move=3))
}




###################################################################
# Update phases
###################################################################

phase.update <- function(E, Sall, Ball, Sig2all, X, Y, GLOBvar, HYPERvar){
  
  # current number of changepoints
  s = length(E) - 2
  
  ### assignement of global variables used here ###
  q = GLOBvar$q
  qmax = GLOBvar$qmax
  Mphase = GLOBvar$Mphase
  nbVarMax = GLOBvar$nbVarMax
  smax = GLOBvar$smax
  qmax = GLOBvar$qmax
  ### end assignement ###

  ### assignement of hyperparameters variables used here ###
  c = HYPERvar$c
  alphalbd = HYPERvar$alphalbd
  betalbd = HYPERvar$betalbd
  alphad2 = HYPERvar$alphad2
  betad2 = HYPERvar$betad2
  v0 = HYPERvar$v0
  gamma0 = HYPERvar$gamma0
  ### end assignement ###
  
  ## For each current phase
  for (phase in E[1:(s+1)]){
  
    ## Position of the phase
    posPhase = which(E==phase)

    ## Parameters
    B = Ball[posPhase,]
    S = (abs(B) > 0) * 1
    k = sum(S)-1
    if(nbVarMax >1){
      Sig2 = Sig2all[posPhase]
    } else {
      Sig2 = Sig2all
    }
    
    ## Observations in the chosen phase
    y = Y[ Mphase[phase]:(Mphase[E[posPhase+1]]-1) ]
    x = X[ Mphase[phase]:(Mphase[E[posPhase+1]]-1), ]
    
    ## Updating hyperparameters
    delta2 = rinvgamma(1, shape=k + alphad2, scale=betad2 + B[which(S==1)] %*% t(x[,which(S==1)]) %*% x[,which(S==1)] %*% B[which(S==1)] / (2*Sig2) ) 
    lambda = rgamma(1, shape=k + alphalbd, rate=1 + betalbd)

    ## Compute acceptation probability vector rho
    rho3 = computeRho3(k, 0, qmax, c, lambda)
   
    ## Sample u
    u = runif(1, 0, 1)

    ## Compute the corresponding move (Edge birth, Edge death or Update the regression coefficient) 
    bduout = bdu(u, rho3, x, y, S, Sig2, delta2, q, v0, gamma0, qmax)
    
    Sall[posPhase,] = bduout$newS
    Ball[posPhase,] = bduout$newB

    ## Update Sig2: case MultiVar
    if(nbVarMax >1){
      Sig2all[posPhase] = updateSigMulti(phase, X, Y, E, Sall, Ball, Sig2, Mphase, alphad2, betad2, v0, gamma0)
    } #end update Sig2
  } # end update each phase

  ## Update Sig2: case UniVar
  if(nbVarMax == 1){
    Sig2all = updateSigSolo(X, Y, E, Sall, Ball, Sig2all, Mphase, alphad2, betad2, v0, gamma0)
  }

  

  ##  Return all variables
  ## (+ variable move describing the move type  (1= CP birth, 2= CP death, 3= CP shift, 4= Update phases)
  return( list( E=E, Sall=Sall, Ball=Ball, Sig2all=Sig2all, move=4, accept=0, move1=bduout$move, accept1=bduout$accept, lambda=lambda))

}


bdu <- function(u, rho3, x, y, S, Sig2, delta2, q, v0, gamma0, qmax){

  ### INPUT:u,rho,s=S[i,],sig2=Sig2[i],delta2.
  ###	x: the data in state i in columns 
  ###	ni: total nb of repeated measurements	
  ### OUTPUT: newS,newSig2,newB.
  ### depends on: 
  ### q the number of predictors
  ### constant v0, gamma0.

  ## Variable move describing the move type  (1= Edge birth, 2= Edge death, 3= Update coefficient, default=3)
  move = 3

  ## Boolean indicating whether the move is accepted or not (=1 if accepted, 0 otherwise, default=0)
  accept = 0

  ## New edges vector, to be returned at the end of the function
  newS = S


  
  ## Current number of edges
  l = sum(S) - 1 # = L[i]

  ## To update Sig2: matPx= Pxl, Pxlp1, or Pxlm1 depending the computed move (birth, death or update).
  ## Compute the projection matrix with the current edge ("Pxl")
  Pxl = computePx(length(y), x[,which(S == 1)], delta2)
 
  if(u < rho3[1] && l < qmax){
    ## Variable move describing the move type  (1= Edge birth, 2= Edge death, 3= Update coefficient)
    move = 1

    ## Sample the additional edge
    sstar = sample(c(which(S==0), which(S==0)), 1) # needed when there is only one position  S==0

    ## Proposed edges vector (with an additional edge)
    stmp = S
    stmp[sstar] = 1

    ## Compute the projection matrix with an additional edge ("Pxl plus 1")
    Pxlp1 = computePx(length(y), x[,which(stmp == 1)], delta2)

    ## Compute birth ratio
    rbirth = ((gamma0 + t(y) %*% Pxlp1 %*% y)/(gamma0 + t(y) %*% Pxl %*% y))^(-(length(y) + v0)/2)/sqrt(1 + delta2)

    ## Sample u 
    u = runif(1,0,1)
                                     
    if(u <= min(1,rbirth)){
      accept = 1
      newS = stmp
    }
    
    ## at this stage : 
    ## newS=S if birth rejected
    ##    =stmp if birth accepted
	
   } else {
    if(u < rho3[2]){
      ## Variable describing the move type  (1 for Edge birth, 2 for Edge death, 3 for Update coefficient)
      move=2
      
      ## Sample the added predictor
      sstar = sample(c(which(S[1:q]==1), which(S[1:q]==1)),1) # needed when there is only one position  S[1:q]==1
      ## Proposed edges vector (after taking away one edge)
      stmp = S
      stmp[sstar] = 0
      
      ## Compute the projection matrix after removimg one edge ("Pxl minus 1")
      Pxlm1 = computePx(length(y), x[,which(stmp==1)], delta2)

      ## Compute death ratio
      rdeath=((gamma0 + t(y) %*% Pxl %*% y)/(gamma0 + t(y) %*% Pxlm1 %*% y))^((length(y) + v0)/2)*(sqrt(1 + delta2))

      ## Sample u 
      u<-runif(1,0,1)
      
      if(u <= min(1,rdeath)){
        ## Boolean for the acceptation of the CP death move (=1 if birth accepted, 0 otherwise)
        accept = 1
        newS = stmp
       }
    }
  }
    
  ## Updating coefficients 
  newB = array(0, q+1)
  if(sum(newS) > 0){
    newB[which(newS == 1)] = sampleBxy(x[, which(newS==1)], y, Sig2, delta2)
  }
 
  ##  Return all variables
  return(list( newS=newS, newB=newB, u=u, move=move, accept=accept))
}







###################################################################
# Update phases for HOMOGENEOUS STRUCTURE
###################################################################

phase.update.homogeneousStructure <- function(E, Sall, Ball, Sig2all, X, Y, GLOBvar, HYPERvar){
  
  # current number of changepoints
  s = length(E) - 2
  
  ### assignement of global variables used here ###
  q = GLOBvar$q
  qmax = GLOBvar$qmax
  Mphase = GLOBvar$Mphase
  nbVarMax = GLOBvar$nbVarMax
  smax = GLOBvar$smax
  qmax = GLOBvar$qmax
  ### end assignement ###

  ### assignement of hyperparameters variables used here ###
  c = HYPERvar$c
  alphalbd = HYPERvar$alphalbd
  betalbd = HYPERvar$betalbd
  alphad2 = HYPERvar$alphad2
  betad2 = HYPERvar$betad2
  v0 = HYPERvar$v0
  gamma0 = HYPERvar$gamma0
  ### end assignement ###
  
  B = Ball[1,]
  S = (abs(B) > 0) * 1
  k = sum(S)-1
  
  scaleForHomogeneousStruct = betad2
  
  ## For each current phase
  for (phase in E[1:(s+1)]){
  
    ## Position of the phase
    posPhase = which(E==phase)

    ## Parameters
    B = Ball[posPhase,]
    #S = (abs(B) > 0) * 1
    #k = sum(S)-1
    if(nbVarMax >1){
      Sig2 = Sig2all[posPhase]
    } else {
      Sig2 = Sig2all
    }
    
    ## Observations in the chosen phase
    y = Y[ Mphase[phase]:(Mphase[E[posPhase+1]]-1) ]
    x = X[ Mphase[phase]:(Mphase[E[posPhase+1]]-1), ]
    
    scaleForHomogeneousStruct= scaleForHomogeneousStruct + B[which(S==1)] %*% t(x[,which(S==1)]) %*% x[,which(S==1)] %*% B[which(S==1)] / (2*Sig2) 
	
  } #end  for (phase in E[1:(s+1)]){
  
    ## Updating hyperparameters
    delta2 = rinvgamma(1, shape=k + alphad2, scale=scaleForHomogeneousStruct)
    lambda = rgamma(1, shape=k + alphalbd, rate=1 + betalbd)

    ## Compute acceptation probability vector rho
    rho3 = computeRho3(k, 0, qmax, c, lambda)
   
    ## Sample u
    u = runif(1, 0, 1)

    ## Compute the corresponding move (Edge birth, Edge death or Update the regression coefficient) 
    bduout = bdu.homogeneousStructure(u, rho3, X, Y, E, Sall, Sig2all, delta2, q, v0, gamma0, qmax, GLOBvar, HYPERvar)
  
  Sall = bduout$Sall
  Ball = bduout$Ball
  Sig2all = bduout$Sig2all

## Removed by Sophie 03/09/2010  
#    ## Update Sig2: case MultiVar
#    if(nbVarMax >1){
#      for (phase in E[1:(s+1)]){
#        
#        posPhase = which(E==phase)
#
#        Sig2all[posPhase] = updateSigMulti(phase, X, Y, E, Sall, Ball, Sig2, Mphase, alphad2, betad2, v0, gamma0)
#      }
#      
#    } #end if(nbVarMax >1){
#
#  ## Update Sig2: case UniVar
#  if(nbVarMax == 1){
#    Sig2all = updateSigSolo(X, Y, E, Sall, Ball, Sig2all, Mphase, alphad2, betad2, v0, gamma0)
#  }

  
  ##  Return all variables
  ## (+ variable move describing the move type  (1= CP birth, 2= CP death, 3= CP shift, 4= Update phases)
  return( list( E=E, Sall=Sall, Ball=Ball, Sig2all=Sig2all, move=4, accept=0, move1=bduout$move, accept1=bduout$accept, lambda=lambda))

}


bdu.homogeneousStructure <- function(u, rho3,X, Y, E, Sall, Sig2all, delta2, q, v0, gamma0, qmax, GLOBvar, HYPERvar){

  ### INPUT:u,rho,s=S[i,],sig2=Sig2[i],delta2.
  ###	x: the data in state i in columns 
  ###	ni: total nb of repeated measurements	
  ### OUTPUT: newS,newSig2,newB.
  ### depends on: 
  ### q the number of predictors
  ### constant v0, gamma0.
  
  s =length(E)-2
  Mphase = GLOBvar$Mphase
  nbVarMax = GLOBvar$nbVarMax
	
  ## Variable move describing the move type  (1= Edge birth, 2= Edge death, 3= Update coefficient, default=3)
  move = 3

  ## Boolean indicating whether the move is accepted or not (=1 if accepted, 0 otherwise, default=0)
  accept = 0

  ## New edges vector, to be returned at the end of the function
  S    = Sall[1,]
  newS = S   
  ## Current number of edges
  l = sum(S) - 1 # = L[i]
  
  # Choose between flip move and other moves
  choice = runif(1, 0, 1)
  
  # Flip Move
  if(l > 0 && choice > 0.75) {    
    # flip move is 4
    move = 4
    
    ## Sample the original parent
    parent.orig = sample(c(which(S[1:q]==1), which(S[1:q]==1)),1) # needed when there is only one position  S[1:q]==1
    # Sample the new parent
    parent.new = sample(c(which(S==0), which(S==0)), 1) # needed when there is only one position  S==0
    
    stmp = S
    stmp[parent.orig] = 0
    stmp[parent.new]  = 1
    
    rflip = 1
    
    for(phase in E[1:(s+1)]) {
      ## Position of the phase
      posPhase = which(E==phase)
    
      ## Observations in the chosen phase
      y = Y[ Mphase[phase]:(Mphase[E[posPhase+1]]-1) ]
      x = X[ Mphase[phase]:(Mphase[E[posPhase+1]]-1), ]
    
      ## To update Sig2: matPx= Pxl, Pxlp1, or Pxlm1 depending the computed move (birth, death or update).
      ## Compute the projection matrix with the current edge ("Pxl")
      Pxl = computePx(length(y), x[,which(S == 1)], delta2)
    
      ## Compute the projection matrix with an additional edge ("Pxl plus 1")
      Pxlp1 = computePx(length(y), x[,which(stmp == 1)], delta2)
      
      rflip = rflip * ((gamma0 + t(y) %*% Pxlp1 %*% y)/(gamma0 + t(y) %*% Pxl %*% y))^(-(length(y) + v0)/2)
    }
    
    u = runif(1,0,1)
    
    if(u <= min(1,rflip)) {
      accept = 1
      newS = stmp
      Sall= t(matrix(newS,q+1, s+1)) 
    }
    
  } else {  
    ##########################
    ## Birth of an edge move
    ##########################
    if(u < rho3[1] && l < qmax){
 
    ## Variable move describing the move type  (1= Edge birth, 2= Edge death, 3= Update coefficient)
    move = 1
    
    ## Sample the additional edge
    sstar = sample(c(which(S==0), which(S==0)), 1) # needed when there is only one position  S==0

    ## Proposed edges vector (with an additional edge)
    stmp = S
    stmp[sstar] = 1
 
    # updated by Sophie for Homogeneous Structure : sum over the posterior ratios
	  #### HERE : rbirth =1 (not 0)
	  rbirth =1
	  ## For each current phase
	  for (phase in E[1:(s+1)]){
    
      ## Position of the phase
      posPhase = which(E==phase)
    
      ## Observations in the chosen phase
      y = Y[ Mphase[phase]:(Mphase[E[posPhase+1]]-1) ]
      x = X[ Mphase[phase]:(Mphase[E[posPhase+1]]-1), ]
    
      ## To update Sig2: matPx= Pxl, Pxlp1, or Pxlm1 depending the computed move (birth, death or update).
      ## Compute the projection matrix with the current edge ("Pxl")
      Pxl = computePx(length(y), x[,which(S == 1)], delta2)
    
      ## Compute the projection matrix with an additional edge ("Pxl plus 1")
      Pxlp1 = computePx(length(y), x[,which(stmp == 1)], delta2)
    
      ## Compute birth ratio
	  ## HERE : * (not +)
      rbirth = rbirth * ((gamma0 + t(y) %*% Pxlp1 %*% y)/(gamma0 + t(y) %*% Pxl %*% y))^(-(length(y) + v0)/2)/sqrt(1 + delta2)
    
	  } # end for (phase in E[1:(s+1)]){
      
	  ## Sample u 
      u = runif(1,0,1)
    
    
      if(u <= min(1,rbirth)){
        accept = 1
        newS = stmp
	    Sall= t(matrix(newS,q+1, s+1)) 
	   # 		print(newS)
	  #	print(Sall)
      }
      
      ## at this stage : 
      ## newS=S if birth rejected
      ##    =stmp if birth accepted
	  
     } else {
      if(u < rho3[2]){
	  #print("removing an edge??")
	    ##########################
	    ## Death of an edge move
        ##########################
        ## Variable describing the move type  (1 for Edge birth, 2 for Edge death, 3 for Update coefficient)
        move=2
        
        ## Sample the edge to remove
        sstar = sample(c(which(S[1:q]==1), which(S[1:q]==1)),1) # needed when there is only one position  S[1:q]==1
        ## Proposed edges vector (after taking away one edge)
        stmp = S
        stmp[sstar] = 0
        
	    # updated by Sophie for Homogeneous Structure : sum over the posterior ratios
	    rdeath =1
	    ## For each current phase
	    for (phase in E[1:(s+1)]){
    
        ## Position of the phase
        posPhase = which(E==phase)
    
        ## Observations in the chosen phase
        y = Y[ Mphase[phase]:(Mphase[E[posPhase+1]]-1) ]
        x = X[ Mphase[phase]:(Mphase[E[posPhase+1]]-1), ]
    
        ## To update Sig2: matPx= Pxl, Pxlp1, or Pxlm1 depending the computed move (birth, death or update).
        ## Compute the projection matrix with the current edge ("Pxl")
        Pxl = computePx(length(y), x[,which(S == 1)], delta2)
    
        ## Compute the projection matrix after removimg one edge ("Pxl minus 1")
        Pxlm1 = computePx(length(y), x[,which(stmp==1)], delta2)
    
        ## Compute death ratio
        rdeath= rdeath*((gamma0 + t(y) %*% Pxl %*% y)/(gamma0 + t(y) %*% Pxlm1 %*% y))^((length(y) + v0)/2)*(sqrt(1 + delta2))
	    } # end for (phase in E[1:(s+1)])
	    
        ## Sample u 
        u<-runif(1,0,1)
         #print(u)
	     #print(rdeath)  
	  
        if(u <= min(1,rdeath)){
          ## Boolean for the acceptation of the CP death move (=1 if birth accepted, 0 otherwise)
          accept = 1
          newS = stmp
	  	Sall= t(matrix(newS,q+1, s+1)) 
	  	#print(newS)
	  	#print(Sall)
         }
      } # end if  if(u < rho3[2]){
    } # end else if(u < rho3[1] && l < qmax){
    }
    	 
  ## Updating coefficients 
  Ball=matrix(0, s+1, q+1)
  newB = array(0, q+1)
  
  ## Moved below by Sophie 03/09/2010
  ## if(sum(newS) > 0){
    
  ## For each current phase
  for (phase in E[1:(s+1)]){
  
      ## Position of the phase
      posPhase = which(E==phase)

      ## Observations in the chosen phase
      y = Y[ Mphase[phase]:(Mphase[E[posPhase+1]]-1) ]
      x = X[ Mphase[phase]:(Mphase[E[posPhase+1]]-1), ]
	
      if(nbVarMax >1){
        Sig2 = Sig2all[posPhase]
      } else {
        Sig2 = Sig2all
      }
	
      ## 'if(sum(newS) > 0){' Added by Sophie 03/09/2010 
      if(sum(newS) > 0){
         # sample newB
         newB[which(newS == 1)] = sampleBxy(x[, which(newS==1)], y, Sig2, delta2)
         Ball[posPhase, ] = as.matrix(newB)
      } # end   if(sum(newS) > 0)
	  
      ## Added by Sophie 03/09/2010
      Sig2all[posPhase] = updateSigMulti(phase, X, Y, E, Sall, Ball, Sig2, Mphase, HYPERvar$alphad2, HYPERvar$betad2, HYPERvar$v0, HYPERvar$gamma0)
  
   } # end for (phase in E[1:(s+1)]){
 
  ##} # end   if(sum(newS) > 0){
  
  ##  Return all variables
  return(list( newS=newS, newB=newB, u=u, move=move, accept=accept, Ball=Ball, Sall=Sall, Sig2all=Sig2all)) 
}



bdu.homogeneousStructure_noflip<- function(u, rho3,X, Y, E, Sall, Sig2all, delta2, q, v0, gamma0, qmax, GLOBvar, HYPERvar){

  ### INPUT:u,rho,s=S[i,],sig2=Sig2[i],delta2.
  ###	x: the data in state i in columns 
  ###	ni: total nb of repeated measurements	
  ### OUTPUT: newS,newSig2,newB.
  ### depends on: 
  ### q the number of predictors
  ### constant v0, gamma0.
  
  s =length(E)-2
  Mphase = GLOBvar$Mphase
  nbVarMax = GLOBvar$nbVarMax
	
  ## Variable move describing the move type  (1= Edge birth, 2= Edge death, 3= Update coefficient, default=3)
  move = 3

  ## Boolean indicating whether the move is accepted or not (=1 if accepted, 0 otherwise, default=0)
  accept = 0

  ## New edges vector, to be returned at the end of the function
  S = Sall[1,]
  newS =S   
  ## Current number of edges
  l = sum(S) - 1 # = L[i]
 ##########################
 ## Birth of an edge move
 ##########################
  if(u < rho3[1] && l < qmax){
 # print("adding an edge??")
    ## Variable move describing the move type  (1= Edge birth, 2= Edge death, 3= Update coefficient)
    move = 1

    ## Sample the additional edge
    sstar = sample(c(which(S==0), which(S==0)), 1) # needed when there is only one position  S==0

    ## Proposed edges vector (with an additional edge)
    stmp = S
    stmp[sstar] = 1
 
    # updated by Sophie for Homogeneous Structure : sum over the posterior ratios
	#### HERE : rbirth =1 (not 0)
	rbirth =1
	## For each current phase
	for (phase in E[1:(s+1)]){
  
    ## Position of the phase
    posPhase = which(E==phase)

    ## Observations in the chosen phase
    y = Y[ Mphase[phase]:(Mphase[E[posPhase+1]]-1) ]
    x = X[ Mphase[phase]:(Mphase[E[posPhase+1]]-1), ]
  
    ## To update Sig2: matPx= Pxl, Pxlp1, or Pxlm1 depending the computed move (birth, death or update).
    ## Compute the projection matrix with the current edge ("Pxl")
    Pxl = computePx(length(y), x[,which(S == 1)], delta2)

    ## Compute the projection matrix with an additional edge ("Pxl plus 1")
    Pxlp1 = computePx(length(y), x[,which(stmp == 1)], delta2)

    ## Compute birth ratio
	## HERE : * (not +)
    rbirth = rbirth * ((gamma0 + t(y) %*% Pxlp1 %*% y)/(gamma0 + t(y) %*% Pxl %*% y))^(-(length(y) + v0)/2)/sqrt(1 + delta2)
 
	} # end for (phase in E[1:(s+1)]){
    
	## Sample u 
    u = runif(1,0,1)
 

    if(u <= min(1,rbirth)){
      accept = 1
      newS = stmp
	  Sall= t(matrix(newS,q+1, s+1)) 
	 # 		print(newS)
	#	print(Sall)
    }
    
    ## at this stage : 
    ## newS=S if birth rejected
    ##    =stmp if birth accepted
	
   } else {
    if(u < rho3[2]){
	#print("removing an edge??")
	  ##########################
	  ## Death of an edge move
      ##########################
      ## Variable describing the move type  (1 for Edge birth, 2 for Edge death, 3 for Update coefficient)
      move=2
      
      ## Sample the edge to remove
      sstar = sample(c(which(S[1:q]==1), which(S[1:q]==1)),1) # needed when there is only one position  S[1:q]==1
      ## Proposed edges vector (after taking away one edge)
      stmp = S
      stmp[sstar] = 0
      
	  # updated by Sophie for Homogeneous Structure : sum over the posterior ratios
	  rdeath =1
	  ## For each current phase
	  for (phase in E[1:(s+1)]){
  
      ## Position of the phase
      posPhase = which(E==phase)

      ## Observations in the chosen phase
      y = Y[ Mphase[phase]:(Mphase[E[posPhase+1]]-1) ]
      x = X[ Mphase[phase]:(Mphase[E[posPhase+1]]-1), ]
  
      ## To update Sig2: matPx= Pxl, Pxlp1, or Pxlm1 depending the computed move (birth, death or update).
      ## Compute the projection matrix with the current edge ("Pxl")
      Pxl = computePx(length(y), x[,which(S == 1)], delta2)

      ## Compute the projection matrix after removimg one edge ("Pxl minus 1")
      Pxlm1 = computePx(length(y), x[,which(stmp==1)], delta2)

      ## Compute death ratio
      rdeath= rdeath*((gamma0 + t(y) %*% Pxl %*% y)/(gamma0 + t(y) %*% Pxlm1 %*% y))^((length(y) + v0)/2)*(sqrt(1 + delta2))
	  } # end for (phase in E[1:(s+1)])
	  
      ## Sample u 
      u<-runif(1,0,1)
       #print(u)
	   #print(rdeath)  
	
      if(u <= min(1,rdeath)){
        ## Boolean for the acceptation of the CP death move (=1 if birth accepted, 0 otherwise)
        accept = 1
        newS = stmp
		Sall= t(matrix(newS,q+1, s+1)) 
		#print(newS)
		#print(Sall)
       }
    } # end if  if(u < rho3[2]){
  } # end else if(u < rho3[1] && l < qmax){
    
    	 
  ## Updating coefficients 
  Ball=matrix(0, s+1, q+1)
  newB = array(0, q+1)
  
  ## Moved below by Sophie 03/09/2010
  ## if(sum(newS) > 0){
    
	## For each current phase
	for (phase in E[1:(s+1)]){
  
      ## Position of the phase
      posPhase = which(E==phase)

      ## Observations in the chosen phase
      y = Y[ Mphase[phase]:(Mphase[E[posPhase+1]]-1) ]
      x = X[ Mphase[phase]:(Mphase[E[posPhase+1]]-1), ]
	
	if(nbVarMax >1){
      Sig2 = Sig2all[posPhase]
    } else {
      Sig2 = Sig2all
    }
	
	  ## 'if(sum(newS) > 0){' Added by Sophie 03/09/2010 
	  if(sum(newS) > 0){
	    # sample newB
        newB[which(newS == 1)] = sampleBxy(x[, which(newS==1)], y, Sig2, delta2)
	    Ball[posPhase, ] = as.matrix(newB)
	  } # end   if(sum(newS) > 0)
	  
	  ## Added by Sophie 03/09/2010
	  Sig2all[posPhase] = updateSigMulti(phase, X, Y, E, Sall, Ball, Sig2, Mphase, HYPERvar$alphad2, HYPERvar$betad2, HYPERvar$v0, HYPERvar$gamma0)
  
    } # end for (phase in E[1:(s+1)]){
 
  ##} # end   if(sum(newS) > 0){
  
  ##  Return all variables
  return(list( newS=newS, newB=newB, u=u, move=move, accept=accept, Ball=Ball, Sall=Sall, Sig2all=Sig2all)) 
}














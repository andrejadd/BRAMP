###############################################
###        tvDBN  :  initialisation        ####
###############################################



## Sample k from a truncated Poisson distribution, (eq 4.5)
sampleK <- function(mini, maxi, lambda, nb){
  # AA: mini: minimal nr. of endges (set to 0)
  #     maxi: max nr. edges (=smax)
  #     lambda: nr. of expected edges
  #     nb: 1
  # it basically says, sample the probability of each k in [mini..maxi] following a poission distribution with given lambda
  # this is the part with prob, then select one k in [mini..maxi] dependent on the probability
  if( mini == maxi) { print("Error with sampling from a truncated Poisson: mini = maxi") }
  out = sample( mini:maxi, nb, replace=TRUE, prob=lambda^(mini:maxi)/apply(matrix(mini:maxi, 1, maxi-mini+1), 2, factorial))
  return(out)
  
}

##############################################
## Sample initial regression coefficient B 

sampleBinit <- function(Si, sig2, delta2, X, q){
  ### INPUT: s=S[i,],sig2=Sig2[i],delta2,
  ###        X the observed data for predictors.
  ###        q number of predictors
  ### OUTPUT: vector newB.
  ### depends on: q the number of predictors.
  
  ## AA22.02.2011, + bias and SAC edge
  newB <- array(0,q+2)
  
  for(l in which(Si == 1)){
    newB[l] <- rnorm(1, mean=0, sd=sqrt(delta2 * sig2 * t(X[,l]) %*% X[,l]))
  }
  return(newB)
}

sampledelta2Global <- function(X, Y, XE, YE, S2Dall, B2Dall, Sig2_2Dall, GLOBvar, HYPERvar, DEBUGLVL2 =F) {

  ## 
  XMphase = GLOBvar$XMphase
  YMphase = GLOBvar$YMphase
  xlocs = GLOBvar$xlocs
  ylocs = GLOBvar$ylocs
  
  scaleForHomogeneousStruct = HYPERvar$betad2

  ## loop over each current segment
  for(xsegid in 1:(length(XE)-1)) {
    
    for(ysegid in 1:(length(YE)-1)) {
      ## extract regress. coefficients
      B = B2Dall[which( B2Dall[,1] == xsegid &  B2Dall[,2] == ysegid),3:dim(B2Dall)[2]]
      
      ## get segment coordinates
      segcoord = c(XMphase[XE[xsegid]], YMphase[YE[ysegid]],XMphase[XE[xsegid+1]]-1, YMphase[YE[ysegid+1]]-1)
      
      if(DEBUGLVL2) {  cat("[extract delta] xsegid: ", xsegid, ", ysegid: ", ysegid, ",  segcoord ")
                       print.table(segcoord) }
      
      ## get the predictor and target data
      x = extractXPredictors(X, segcoord, xlocs,F)
      y = extractYTargets(Y, segcoord, xlocs,F)
      
      scaleForHomogeneousStruct= scaleForHomogeneousStruct + B[which(S2Dall==1)] %*% t(x[,which(S2Dall==1)]) %*% x[,which(S2Dall==1)] %*% B[which(S2Dall==1)] / (2*Sig2_2Dall) 
      	
    }
  }
  
  nrsegs = (length(XE)-1) * (length(YE)-1)

  ## Updating hyperparameters
  ## AA: do not 'sum(S2Dall) - 1' because the bias also has to be considered
  ## the shape is the count for the dimension of the projection matrix (nr. of counts), it must match in the following way
  ##     sum(S2Dall) * nrsegs = sum(dim((x[,which(S2Dall==1)]))[,2])  , where dim()[,2] gives the column count which is nr. parents + bias
  delta2 = rinvgamma(1, shape= sum(S2Dall) * nrsegs + HYPERvar$alphad2, scale=scaleForHomogeneousStruct)

  return(delta2)



}

##############################################
## Sample parameters

sampleParms <- function(X, GLOBvar, HYPERvar, DEBUGLVL1=F){
  ### assignement of global variables used here ###
  smax = GLOBvar$smax
  q = GLOBvar$q
  smax = GLOBvar$smax
  kmax = GLOBvar$kmax
  n = GLOBvar$n
  xlocs = GLOBvar$xlocs
  ylocs = GLOBvar$ylocs
  XMphase = GLOBvar$XMphase
  YMphase = GLOBvar$YMphase
  nbVarMax = GLOBvar$nbVarMax
  dyn = GLOBvar$dyn
  minPhase = GLOBvar$minPhase
  ### end assignement ###
  
  ### assignement of hyperparameters variables used here ###
  alphaD = HYPERvar$alphaD
  betaD = HYPERvar$betaD
  alphalbd = HYPERvar$alphalbd
  betalbd = HYPERvar$betalbd
  v0 = HYPERvar$v0
  gamma0 = HYPERvar$gamma0
  alphad2 = HYPERvar$alphad2
  betad2 = HYPERvar$betad2
  delta2 = HYPERvar$delta2
### end assignement ###

  ## sample the number of expected (mean) changepoints with rgamma and then the number of changepoints for the x and y axis
  k_x <- sampleK(0,kmax, rgamma(1, shape=alphaD, rate = betaD),1) # scale s= 1/rate => f(x)= 1/(s^a Gamma(a)) x^(a-1) e^-(x/s)
  k_y <- sampleK(0,kmax, rgamma(1, shape=alphaD, rate = betaD),1)
    
  ## AA: init CP vectors for x and y axis
  XE = c(1+dyn, xlocs+1)
  YE = c(1+dyn, ylocs+1)

  ## sample the position of the changepoints k
  cpt = k_x
  while(cpt > 0){
    ## search for possible CP, not in E and not close to E if minPhase (length of phase) is > than 1
    toremove = XE

    if(minPhase>1) for(i in 1:(minPhase-1)) toremove = c(toremove, XE-i, XE+i)

    ## possibles CPs are those not in 'toremove'
    possibleCP = setdiff((1+dyn):XE[length(XE)], toremove)

    ## make sure we dont fail here in case for any reason we have an empty possibleCP list
    if(length(possibleCP)== 0) {
      break;
    }

    ## sample one CP in possibleCP (the vector is double for sake of function sample when size is = to 1)
    cp = sample( c(possibleCP, possibleCP), 1)
    
    XE=sort(c(XE, cp))
    cpt = cpt-1
  }

  ## do the same for the y axis
  cpt = k_y
  
  while(cpt > 0){
    ## search for possible CP, not in E and not close to E if minPhase (length of phase) is > than 1
    toremove = YE

    if(minPhase>1) for(i in 1:(minPhase-1)) toremove = c(toremove, YE-i, YE+i)

    ## possibles CPs are those not in 'toremove'
    possibleCP = setdiff((1+dyn):YE[length(YE)], toremove)

    ## make sure we dont fail here in case for any reason we have an empty possibleCP list
    if(length(possibleCP)== 0) {
      break;
    }

    ## sample one CP in possibleCP (the vector is double for sake of function sample when size is = to 1)
    cp = sample( c(possibleCP, possibleCP), 1)
    
    YE=sort(c(YE, cp))
    cpt = cpt-1
  }


  if(DEBUGLVL1) { cat("XE: "); print.table(XE);
                cat("YE: "); print.table(YE);
                }


  
  ## sample model structure, only one structure because this is homogeneous
  S = matrix(0, 1, q+2)

  ## sample mean nr. edges (Lambda) with rgamma and sampleK (is the same inverse gamma for CPs and edges) 
  sPred = sampleK(0, smax, rgamma(1, shape=alphalbd, rate = betalbd), 1)    # scale s= 1/rate => f(x)= 1/(s^a Gamma(a)) x^(a-1) e^-(x/s)

  # if there are any edges..
  if(sPred>0){
    # set random position in Structure S to edge (1)
    S[1, sample(1:q, sPred, replace=FALSE)] = array(1, sPred) # structure du model (1 si pred in the model)
  }

  ## we assume that there is a constant in each model, this is the last element that corresponds to the bias
  S[, q+1] = array(1,1)
  S[, q+2] = array(1,1)

  ## Here start the real interesting segment dependend data, first read out how many segments we have, then give each segment a id
  ## Finally each segment specific parameter (set) is assigned a matrix where the first column states the segment id and the rest are the parameter(s)

  ## create the variance sigma - one global value  : IG(v0/2,gamma0/2), note this is not the real sigma prior because we lack the knowledge of regre. at this point
  Sig2_2Dall = rinvgamma(n=1, shape=v0/2, scale=gamma0/2)
  
  ## create the regr. coefficient matrix: c(xseg_id, yseg_id, .., regr. coeff,.., bias, SAC )
  B = matrix(0,nrow=0,q+4)
  
  ## loop over each possible segment and sample parameters
  for(xsegid in 1:(length(XE)-1)) {

    for(ysegid in 1:(length(YE)-1)) {
    
      ## get segment coordinates
      segcoord = c(XMphase[XE[xsegid]], YMphase[YE[ysegid]],XMphase[XE[xsegid+1]]-1, YMphase[YE[ysegid+1]]-1)
        
      ## get the predictor data
      x = extractXPredictors(X, segcoord, xlocs,F)

      ## this works only if there is no CP inbetween, E is a helper here but everything should be considered using XE and YE
      newB = sampleBinit( S[1,], Sig2_2Dall, delta2, x , q)

      B = rbind(B, c(xsegid, ysegid, newB))
    }
    
  }

  
  return(list(XE=XE, YE=YE, S2Dall=S[1,], B2Dall=B, Sig2_2Dall=Sig2_2Dall))
}



sampleBxy <- function(xi, y, Sig2_2Dall, delta2){
 
    tryCatch({

      Ml = (delta2 / (delta2+1)) * pseudoinverse(t(xi) %*% xi)
      out = mvrnorm(1, mu=Ml %*% t(xi) %*% y, Sigma=Sig2_2Dall * Ml)
	
    }, error = function(e) {
	cat("Caught error while doing sampleBxy: ")
       	print(e)

        cat("delta2: ", delta2, "\n") 

	cat("saving all parameters to file debugsave1: ")
	save(delta2, xi, y, Sig2, Ml, file="debugsave1")

    	stop("DEBUG-E1: failed because of wrong Sigma , see output")    

    })

  return(out)
}


## AA, change arguments, instead of extracting from Sall, Ball, X, Y - directly pass from moves.R since they are already used there
updateSigGlobal <- function(xlocs,XMphase, YMphase, XE, YE, X, Y, S, delta2, v0, gamma0){

  sumP = 0

  ## loop over each current segment
  for(xsegid in 1:(length(XE)-1)) {

    for(ysegid in 1:(length(YE)-1)) {

      ## get segment coordinates
      segcoord = c(XMphase[XE[xsegid]], YMphase[YE[ysegid]],XMphase[XE[xsegid+1]]-1, YMphase[YE[ysegid+1]]-1)

      ## get the predictor and target data
      x = extractXPredictors(X, segcoord, xlocs,F)
      y = extractYTargets(Y, segcoord, xlocs,F)

      matPx = computePx(length(y), x[, which(S == 1)], delta2)
    
      sumP = sumP + (t(y) %*% matPx %*% y)

    }
  }

  ## Inverse Gamma
  out = rinvgamma(1, shape=(v0 + length(Y)) / 2, scale= (gamma0 + sumP) / 2)
  return(out)
  
}






###-------------------- OLD FUNCTIONS ----------------------------------------------------------------

# This is the signal to noise ratio
old_sampleDelta2 <- function(pos, x, q, B, S, sig2, alphad2, betad2){
  # INPUT: pos, the considered state
  #        xPos, the observations of X in state i
  #        B,S,Sig2
  #        alphad2,betad2.
  # OUTPUT: delta2
  # depends on: . 
  plus = 0

  # extract Structure and Regr. Coeff., might be also done in calling code
  currS = S[which(S[,1] == pos), 2:dim(S)[2]]
  currB = B[which(B[,1] == pos), 2:dim(B)[2]]

    
  if(sum(currS)>0){
    Bi = currB[which(currS == 1)]
    xi = x[, which(currS == 1)]
    plus = Bi %*% t(xi) %*% xi %*% Bi / (2 * sig2)
  }
  
  out = rinvgamma(1, shape=sum(currS[1:q]) + alphad2, scale=betad2 + plus)
  return(out)
}

## AA, change arguments, instead of extracting from Sall, Ball, X, Y - directly pass from moves.R since they are already used there
old_updateSigMulti <- function(x, y, S, B, Sig2, delta2, alphad2, betad2, v0, gamma0){

  s = sum(S)-1

  ## AA, make no delta2 update here , but the global update right before calling this function
  ## 1D
  ##  delta2 = rinvgamma(1, shape= s + alphad2, scale=betad2 + Ball[posPhase, which(S == 1)] %*% t(x[, which(S == 1)]) %*% x[,which(S == 1)] %*% Ball[posPhase,which(S == 1)] / (2 * Sig2) )

  ##  delta2 = rinvgamma(1, shape= s + alphad2, scale=betad2 + B[which(S == 1)] %*% t(x[, which(S == 1)]) %*% x[,which(S == 1)] %*% B[which(S == 1)] / (2 * Sig2) )
  
  matPx = computePx(length(y), x[, which(S == 1)], delta2)
    
  total = t(y) %*% matPx %*%y

  ## Inverse Gamma
  out = rinvgamma(1, shape=v0/2 + length(y)/2, scale=(gamma0 + total)/2)
  return(out)
  
}


old_updateSigSolo <- function(X, Y, E, Sall, Ball, Sig2, Mphase, alphad2, betad2, v0, gamma0){
 
  total = 0
  for (phase in E[1:(length(E)-1)]){
    posPhase = which(E == phase)
    S = Sall[posPhase,]
    k = sum(Sall[posPhase,])-1
    #new definition of Mphase: -dyn not required!
    y = Y[(Mphase[phase]:(Mphase[E[posPhase+1]]-1))]
    x = X[(Mphase[phase]:(Mphase[E[posPhase+1]]-1)),]
    delta2 = rinvgamma(1, shape=k + alphad2, scale=betad2 + Ball[posPhase, which(S == 1)] %*% t(x[, which(S == 1)]) %*% x[, which(S == 1)] %*% Ball[posPhase, which(S == 1)] / (2 * Sig2) )
    matPx = computePx(length(y), x[,which(S == 1)], delta2)
    total = total + t(y) %*% matPx %*%y
  }
  
  newSig2=rinvgamma(1, shape=v0/2+length(y)/2, scale = (gamma0+total)/2)
  ## newSig2=rinvgamma(1, shape=v0/2+length(y)/(length(E)-1)/2, scale = (gamma0+total)/(length(E)-1)/2)
  ## newSig2=rinvgamma(1, shape=v0/2+length(y)/2, scale = (gamma0+total/(length(E)-1))/2)
  ## newSig2=rinvgamma(1, shape=v0/2+ mean(E[2:length(E)]-E[1:(length(E)-1)])*m/2, scale = (gamma0+total/(length(E)-1))/2)

  ## mean( rinvgamma(100, shape=v0/2+length(y)/2, scale = (gamma0+total)/2))
  ## var( rinvgamma(100, shape=v0/2+length(y)/2, scale = (gamma0+total)/2))
  ## plot.density(density(rinvgamma(100, shape=v0/2, scale = gamma0/2)))          
  ## mean( rinvgamma(100, shape=v0/2+length(y)/2, scale = (gamma0+total/(length(E)-1))/2))
  ## var( rinvgamma(100, shape=v0/2+length(y)/2, scale = (gamma0+total/(length(E)-1))/2))
  
  return(newSig2)
}


old_updateSig <- function(u, rho, X, Y, Sall, Ball, Sig2all){
  ################ ATTENTION ##################
  ## VARIABLES GLOBALES : target, n, Mphase, alphad2, betad2, v0, gamma0
  
  newSig2 = Sig2all
  #for (target in 1:p){
  vect = ((target-1) * n+1):(target * n)
  total = 0 
  for (phase in Eall[which(Eall %in% vect)]) {
    posPhase = which(Eall == phase)
    S = Sall[posPhase,]
    k = sum(Sall[posPhase,])-1
    y = Y[Mphase[phase]:(Mphase[Eall[posPhase+1]]-1)]
    x = X[Mphase[phase]:(Mphase[Eall[posPhase+1]]-1),]
    delta2 = rinvgamma(1, shape=k + alphad2, scale=betad2 + Ball[posPhase, which(S==1)] %*% t(x[, which(S == 1)]) %*% x[,which(S == 1)] %*% Ball[posPhase,which(S == 1)] / (2 * Sig2all[target]))
    matPx = computePx(length(y), x[,which(S == 1)], delta2)
    
    total = total + t(y) %*% matPx %*% y
  }
  newSig2[target] = rinvgamma(1, shape=v0/2 + length(y)/2, scale=(gamma0 + total / sum(Eall %in% vect))/2)
#}
#
  return(newSig2)
}

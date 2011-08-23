# function that read input data, either as matrix or as file name
readInput <- function(data){
  # if data is character, i.e. a file path,
  if(is.character(data)){
    # if file exists
    if(file.exists(data)){
       DATA = as.matrix(read.table(data))
       return(DATA)
     } else {
       stop(paste("input file path",data,"is not valid\n"))
     }
  }
  # if data is already a matrix
  if(is.matrix(data)){
    return(as.matrix(data))
  }

  ### AA:ADD_START
  # in case it is a data.frame, as it will be when read.table is used before
  if(is.data.frame(data)) {
    return(as.matrix(data))
  }
  ### AA:ADD_END
  
  # otherwise input data is incorrect :
  stop(paste("input data ",data," is not valid\n"))

}


# function which show teh available prior distributions for kmax given as a parameter
choosePriors<-function(kmax,priorsPath=paste(priorsPath,"k_priors.txt",sep="")){
	priors=read.table(priorsPath,header=T)
	index=which(priors[,1]==kmax)
	
	if(length(index)>1){
		par(mfrow=(c(ceiling(length(index)/2),2)),cex=1)
		for(i in index[2:1]){
			plot(0:kmax,priors[i,4:(kmax+4)],type="h",lwd=5,col=2,main=paste("alpha=",priors[i,2],", beta= ", priors[i,3]),ylab="Prior probability",xlab="Number of changepoints or TF")
		}
		for(i in index[3:length(index)]){
			plot(0:kmax,priors[i,4:(kmax+4)],type="h",lwd=5,col=4,main=paste("alpha=",priors[i,2],", beta= ", priors[i,3]),ylab="Prior probability",xlab="Number of changepoints or TF")
		}
	}else{
		par(mfrow=(c(1,1)),cex=1)
		plot(0:kmax,priors[index,4:(kmax+4)],type="h",lwd=5,col=2,main=paste("alpha=",priors[index,2],", beta= ", priors[index,3]),ylab="Prior probability",xlab="Number of changepoints or TF")
	}
	
}


## AA: called from main.R
computeRho4 = function(k, kmin, kmax, c, lambda){
  # INPUT:l= k the number of hidden states
  #          kmax the maximal number of hidden states
  #          c constante
  #          lambda simulation parameter for the number of hidden states
  # OUTPUT: rho
  # depends on: (/)

  
  rho=array(1,4)
  if(k == kmax) { rho[1] = 0 } else { rho[1] = c*min(1,lambda/(k+1)) }
  if(k == kmin) { rho[2] = rho[1] } else { rho[2] = rho[1]+c*min(1,k/lambda) }
  if(k > 0) { rho[3] = rho[2]+(1-rho[2])/3 } else{ rho[3] = rho[2] }

  return(rho)
}

computeRho3 = function(k, kmin, kmax, c, lambda){
  # INPUT:l= k the number of hidden states
  #          kmax the maximal number of hidden states
  #          c constante
  #          lambda simulation parameter for the numbr of hidden states
  # OUTPUT: rho
  # depends on: (/)

  rho=array(1,3)
  if(k == kmax) { rho[1]=0 } else { rho[1] = c*min(1,lambda/(k+1)) }
  if(k == kmin) { rho[2] = rho[1] } else { rho[2] = rho[1]+c*min(1,k/lambda) }

  return(rho)
}


bp.computeHammingAlpha<-function(birth,lNew,kminus,Ekl,Estar,Ekr,yL,PxL,yR,PxR,y2,Px2,D,delta2, q, smax, v0, gamma0, ham_dist){
  # birth = 1 for birth, -1 for death.
  # lNew = number of edges in the new phase
  # kminus = minimal number of changepoints between the 2 compared models (=s for birth, s-1 for death) -> INUTILE !!!
  # Ekl = 
  # Estar =
  # yL,yR, y2 : response data  (left, right, both)
  # PxL, PxR, Px2 : projection matrix (left, right, both)
  # D : hyperparms for the number of edges in each phase. (Number of edges s ~ truncated Poisson P(D). )
  # delta2 : hyperparms for empirical covariance (can be seen as the expected  signal-to-noise ratio)  ~IG(alphad2,betad2)

  # modifie par Sophie 01/03/09 
  # Motif : la function gamma() prend en parms des valeurs < 172 -> remplacement par lgamma().
  # LR=factorial(q-lNew)/factorial(q)/(sqrt(delta2+1))^(lNew+1)*gamma(((Estar-Ekl)+v0)/2)*gamma(((Ekr-Estar)+v0)/2)/gamma(((Ekr-Ekl)+v0)/2)*(((gamma0+t(y2) %*% Px2 %*% y2)/2)/((gamma0+t(yL) %*% PxL %*%yL)/2))^(((Estar-Ekl)+v0)/2)* (((gamma0+t(y2) %*% Px2 %*% y2)/2)/((gamma0+t(yR) %*% PxR %*% yR)/2))^(((Ekr-Estar))/2)/((gamma0+t(yR) %*% PxR %*% yR)/2)^(v0/2)*(gamma0/2)^(v0/2)/gamma(v0/2)/(sum(D^(0:smax)/factorial(0:smax)))
  # res = (LR)^birth
  # return(min(1,res))

  # modifie par Sophie 22/05/09 
  # logR=log(factorial(q-lNew)/factorial(q)/(sqrt(delta2+1))^(lNew+1))+lgamma(((Estar-Ekl)+v0)/2)+lgamma(((Ekr-Estar)+v0)/2)-lgamma(((Ekr-Ekl)+v0)/2)+(((Estar-Ekl)+v0)/2)*log(((gamma0+t(y2) %*% Px2 %*% y2)/2)/((gamma0+t(yL) %*% PxL %*%yL)/2))+(((Ekr-Estar))/2)*log (((gamma0+t(y2) %*% Px2 %*% y2)/2)/((gamma0+t(yR) %*% PxR %*% yR)/2))-(v0/2)*log((gamma0+t(yR) %*% PxR %*% yR)/2)+(v0/2)*log(gamma0/2)-lgamma(v0/2)-log(sum(D^(0:smax)/factorial(0:smax)))+lNew+log(D)
  logR=log(factorial(q-lNew)/factorial(q)/(sqrt(delta2+1))^(lNew+1)*D^lNew)+
    lgamma(((Estar-Ekl)+v0)/2)+lgamma(((Ekr-Estar)+v0)/2)-
    lgamma(((Ekr-Ekl)+v0)/2)+(((Estar-Ekl)+v0)/2)*log(((gamma0+t(y2) %*% 
    Px2 %*% y2)/2)/((gamma0+t(yL) %*% PxL %*%yL)/2))+(((Ekr-Estar))/2)*
    log(((gamma0+t(y2) %*% Px2 %*% y2)/2)/((gamma0+t(yR) %*% PxR %*% yR)/2))-
    (v0/2)*log((gamma0+t(yR) %*% PxR %*% yR)/2)+(v0/2)*log(gamma0/2)-
    lgamma(v0/2)-log(sum(D^(0:smax)/factorial(0:smax)))

    #test = - dpois(ham_dist, 0.5, log=TRUE) + log(choose(q, ham_dist));

  logR = logR - dpois(ham_dist, 0.5, log=TRUE) + log(choose(q, ham_dist))
  
  logR=birth*logR


  if(logR>0){
    res=1
  } else {
    res=exp(logR)
  }
  
#### *D^lNew ?
  return(res)
}

cp.computeAlpha <- function(birth, X, Y, xlocs, ylocs, ALTERX, XMphase, YMphase, E, Eplus, OtherE,  poskl, HYPERvar, S2Dall, D, DEBUG_BIRTH_EXT=F) {

  gamma0 = HYPERvar$gamma0
  v0 = HYPERvar$v0
  delta2 = HYPERvar$delta2
  

  ###### product of posterior prob of current state : Phi_h
  ##
  #prodPhi = 1     # product over current state
  sumPhi = 0
  
  ## 1. extract parameters that are in the current state segments
  if(ALTERX) { xsegid = poskl; tmpYE = OtherE; tmpXE = E } else { ysegid = poskl; tmpXE = OtherE; tmpYE = E  }
  
  for(i in 1:(length(OtherE)-1)) {

    ## Assign proper segment id
    if(ALTERX) { ysegid = i } else { xsegid = i }
    
     ## get segment coordinates
    segcoord = c(XMphase[tmpXE[xsegid]], YMphase[tmpYE[ysegid]],XMphase[tmpXE[xsegid+1]]-1, YMphase[tmpYE[ysegid+1]]-1)


    if(DEBUG_BIRTH_EXT == TRUE) {  cat("[prodPhi] xsegid: ", xsegid, ", ysegid: ", ysegid, ",  segcoord ")
                                  print.table(segcoord) }
    
    ## get the predictor data
    x = extractXPredictors(X, segcoord, xlocs,F)

    ## get the target data
    y = extractYTargets(Y, segcoord, xlocs,F)
   
    ## number of locations
    omega = length(y)

    ## calculate projection matrix
    Pr = computePx(omega, as.matrix(x[,which(S2Dall == 1)]), delta2)

    #prodPhi = prodPhi * gamma((v0+omega)/2) * ((gamma0+ t(y) %*% Pr %*% y)/2)^(-(v0+omega)/2)
    sumPhi  = sumPhi  + lgamma((v0+omega)/2) + (-(v0+omega)/2) * log( (gamma0+ t(y) %*% Pr %*% y)/2)

  
  }


  ###### product of posterior prob of next state : Phi^+_h
  ##
  #prodPhiPlus = 1 # over next state
  sumPhiPlus = 0

  ## Assign proper CP vector to temporary type (used in segcoord extraction)
  if(ALTERX) {tmpYE = OtherE; tmpXE = Eplus } else { tmpXE = OtherE; tmpYE = Eplus }
  
  for(j in poskl:(poskl+1)) {

    for(i in 1:(length(OtherE)-1)) {

      ## Assign proper segment id
      if(ALTERX) { xsegid = j; ysegid = i } else { xsegid = i; ysegid = j}
      
      ## get segment coordinates
      segcoord = c(XMphase[tmpXE[xsegid]], YMphase[tmpYE[ysegid]],XMphase[tmpXE[xsegid+1]]-1, YMphase[tmpYE[ysegid+1]]-1)

      if(DEBUG_BIRTH_EXT == TRUE) {  cat("[prodPhiPlus] xsegid: ", xsegid, ", ysegid: ", ysegid, ",  segcoord ")
                                     print.table(segcoord) }
      ## get the predictor data
      x = extractXPredictors(X, segcoord, xlocs,F)
      
      ## get the target data
      y = extractYTargets(Y, segcoord, xlocs,F)
      
      ## number of locations
      omega = length(y)
      
      ## calculate projection matrix
      Pr = computePx(omega, as.matrix(x[,which(S2Dall == 1)]), delta2)

    #  prodPhiPlus = prodPhiPlus * gamma((v0+omega)/2) * ((gamma0+ t(y) %*% Pr %*% y)/2)^(-(v0+omega)/2)
      sumPhiPlus  = sumPhiPlus  + lgamma((v0+omega)/2) + (-(v0+omega)/2) * log( (gamma0+ t(y) %*% Pr %*% y)/2)


    }
  }

  ## nr edges s
  s = sum(S2Dall) - 1

  ## nr CPs k
  k = length(E) - 2

  ## this is the exponent |Phi+| - |Phi| = |Phi| = k of other axis, in eq. (15), since |Phi+| = 2|Phi| and |Phi_h| = |CPs other axis|   
  dSegmentNr = length(OtherE) - 1

#  alpha1a =  D / (c - 1 - k) * ( (gamma0/2)^(v0/2) / (gamma(v0/2) * (delta2+1)^((s+1)/2) ))
#  alpha1b = prodPhiPlus / prodPhi
  
#  alpha2a =  log(D / (c - 1 - k))+(v0/2)*log(gamma0/2) - lgamma(v0/2) - log((delta2+1)^((s+1)/2)) 
#  alpha2b = sumPhiPlus - sumPhi

#  alpha1 = D / (c - 1 - k) * ( (gamma0/2)^(v0/2) / (gamma(v0/2) * (delta2+1)^((s+1)/2) ))^dSegmentNr *  prodPhiPlus / prodPhi

  if(ALTERX) { c = xlocs } else { c = ylocs}

  alpha = birth * ( log(D / (c - 1 - k)) + dSegmentNr *( (v0/2)*log(gamma0/2) - lgamma(v0/2) - log((delta2+1)^((s+1)/2))) + sumPhiPlus - sumPhi )

  return(exp(alpha))


}

bp.computeAlpha_updated<-function(birth,lNew,Ekl,Estar,Ekr,yL,PxL,yR,PxR,y2,Px2,D,delta2, q, smax, v0, gamma0){
  # birth = 1 for birth, -1 for death.
  # lNew = number of edges in the new phase
  # Ekl = 
  # Estar =
  # yL,yR, y2 : response data  (left, right, both)
  # PxL, PxR, Px2 : projection matrix (left, right, both)
  # D : hyperparms for the number of edges in each phase. (Number of edges s ~ truncated Poisson P(D). )
  # delta2 : hyperparms for empirical covariance (can be seen as the expected  signal-to-noise ratio)  ~IG(alphad2,betad2)

  # modifie par Sophie 16/09/09 
  # -log(sum(D^(0:smax)/factorial(0:smax)))+log(factorial(q-lNew)/factorial(q)*D^lNew)


  logR=  +(v0/2)*log(gamma0/2)-lgamma(v0/2)  -log( (sqrt(delta2+1))^(lNew+1) ) +lgamma(((Estar-Ekl)+v0)/2)+lgamma(((Ekr-Estar)+v0)/2)-lgamma(((Ekr-Ekl)+v0)/2) +(((Estar-Ekl)+v0)/2)*log(((gamma0+t(y2) %*% Px2 %*% y2)/2)/((gamma0+t(yL) %*% PxL %*%yL)/2))  +(((Ekr-Estar))/2)*log (((gamma0+t(y2) %*% Px2 %*% y2)/2)/((gamma0+t(yR) %*% PxR %*% yR)/2))  -(v0/2)*log((gamma0+t(yR) %*% PxR %*% yR)/2)

  logR=birth*logR

  if(logR>0){
    res=1
  } else {
    res=exp(logR)
  }
  
  return(res)
}

bp.computePrecompAlpha<-function(birth,lNew,kminus,Ekl,Estar,Ekr,yL,PxL,yR,PxR,y2,Px2,D,delta2, q, smax, v0, gamma0, prop_prob){
  # birth = 1 for birth, -1 for death.
  # lNew = number of edges in the new phase
  # kminus = minimal number of changepoints between the 2 compared models (=s for birth, s-1 for death) -> INUTILE !!!
  # Ekl = 
  # Estar =
  # yL,yR, y2 : response data  (left, right, both)
  # PxL, PxR, Px2 : projection matrix (left, right, both)
  # D : hyperparms for the number of edges in each phase. (Number of edges s ~ truncated Poisson P(D). )
  # delta2 : hyperparms for empirical covariance (can be seen as the expected  signal-to-noise ratio)  ~IG(alphad2,betad2)

  # modifie par Sophie 01/03/09 
  # Motif : la function gamma() prend en parms des valeurs < 172 -> remplacement par lgamma().
  # LR=factorial(q-lNew)/factorial(q)/(sqrt(delta2+1))^(lNew+1)*gamma(((Estar-Ekl)+v0)/2)*gamma(((Ekr-Estar)+v0)/2)/gamma(((Ekr-Ekl)+v0)/2)*(((gamma0+t(y2) %*% Px2 %*% y2)/2)/((gamma0+t(yL) %*% PxL %*%yL)/2))^(((Estar-Ekl)+v0)/2)* (((gamma0+t(y2) %*% Px2 %*% y2)/2)/((gamma0+t(yR) %*% PxR %*% yR)/2))^(((Ekr-Estar))/2)/((gamma0+t(yR) %*% PxR %*% yR)/2)^(v0/2)*(gamma0/2)^(v0/2)/gamma(v0/2)/(sum(D^(0:smax)/factorial(0:smax)))
  # res = (LR)^birth
  # return(min(1,res))

  # modifie par Sophie 22/05/09 
  # logR=log(factorial(q-lNew)/factorial(q)/(sqrt(delta2+1))^(lNew+1))+lgamma(((Estar-Ekl)+v0)/2)+lgamma(((Ekr-Estar)+v0)/2)-lgamma(((Ekr-Ekl)+v0)/2)+(((Estar-Ekl)+v0)/2)*log(((gamma0+t(y2) %*% Px2 %*% y2)/2)/((gamma0+t(yL) %*% PxL %*%yL)/2))+(((Ekr-Estar))/2)*log (((gamma0+t(y2) %*% Px2 %*% y2)/2)/((gamma0+t(yR) %*% PxR %*% yR)/2))-(v0/2)*log((gamma0+t(yR) %*% PxR %*% yR)/2)+(v0/2)*log(gamma0/2)-lgamma(v0/2)-log(sum(D^(0:smax)/factorial(0:smax)))+lNew+log(D)
  logR=log(factorial(q-lNew)/factorial(q)/(sqrt(delta2+1))^(lNew+1)*D^lNew)+
    lgamma(((Estar-Ekl)+v0)/2)+lgamma(((Ekr-Estar)+v0)/2)-
    lgamma(((Ekr-Ekl)+v0)/2)+(((Estar-Ekl)+v0)/2)*log(((gamma0+t(y2) %*% 
    Px2 %*% y2)/2)/((gamma0+t(yL) %*% PxL %*%yL)/2))+(((Ekr-Estar))/2)*
    log(((gamma0+t(y2) %*% Px2 %*% y2)/2)/((gamma0+t(yR) %*% PxR %*% yR)/2))-
    (v0/2)*log((gamma0+t(yR) %*% PxR %*% yR)/2)+(v0/2)*log(gamma0/2)-
    lgamma(v0/2)-log(sum(D^(0:smax)/factorial(0:smax)))

    #test1 = log(factorial(q-lNew)/factorial(q)) + log(D^lNew);
    #test2 = -log(sum(D^(0:smax)/factorial(0:smax)))
   
    #browser()
	logR = logR - log(prop_prob);
	
  logR=birth*logR


  if(logR>0){
    res=1
  } else {
    res=exp(logR)
  }
  
#### *D^lNew ?
  return(res)
}

bp.computeAlpha<-function(birth,lNew,kminus,Ekl,Estar,Ekr,yL,PxL,yR,PxR,y2,Px2,D,delta2, q, smax, v0, gamma0){
  # birth = 1 for birth, -1 for death.
  # lNew = number of edges in the new phase
  # AA: not used: kminus = minimal number of changepoints between the 2 compared models (=s for birth, s-1 for death) -> INUTILE !!!
  # Ekl = 
  # Estar =
  # Ekr = 
  # yL,yR, y2 : response data  (left, right, both)
  # PxL, PxR, Px2 : projection matrix (left, right, both)
  # D : hyperparms for the number of edges in each phase. (Number of edges s ~ truncated Poisson P(D). )
  # delta2 : hyperparms for empirical covariance (can be seen as the expected  signal-to-noise ratio)  ~IG(alphad2,betad2)

  # modifie par Sophie 01/03/09 
  # Motif : la function gamma() prend en parms des valeurs < 172 -> remplacement par lgamma().
  # LR=factorial(q-lNew)/factorial(q)/(sqrt(delta2+1))^(lNew+1)*gamma(((Estar-Ekl)+v0)/2)*gamma(((Ekr-Estar)+v0)/2)/gamma(((Ekr-Ekl)+v0)/2)*(((gamma0+t(y2) %*% Px2 %*% y2)/2)/((gamma0+t(yL) %*% PxL %*%yL)/2))^(((Estar-Ekl)+v0)/2)* (((gamma0+t(y2) %*% Px2 %*% y2)/2)/((gamma0+t(yR) %*% PxR %*% yR)/2))^(((Ekr-Estar))/2)/((gamma0+t(yR) %*% PxR %*% yR)/2)^(v0/2)*(gamma0/2)^(v0/2)/gamma(v0/2)/(sum(D^(0:smax)/factorial(0:smax)))
  # res = (LR)^birth
  # return(min(1,res))

  # modifie par Sophie 22/05/09 
  # logR=log(factorial(q-lNew)/factorial(q)/(sqrt(delta2+1))^(lNew+1))+lgamma(((Estar-Ekl)+v0)/2)+lgamma(((Ekr-Estar)+v0)/2)-lgamma(((Ekr-Ekl)+v0)/2)+(((Estar-Ekl)+v0)/2)*log(((gamma0+t(y2) %*% Px2 %*% y2)/2)/((gamma0+t(yL) %*% PxL %*%yL)/2))+(((Ekr-Estar))/2)*log (((gamma0+t(y2) %*% Px2 %*% y2)/2)/((gamma0+t(yR) %*% PxR %*% yR)/2))-(v0/2)*log((gamma0+t(yR) %*% PxR %*% yR)/2)+(v0/2)*log(gamma0/2)-lgamma(v0/2)-log(sum(D^(0:smax)/factorial(0:smax)))+lNew+log(D)

  logR = log( factorial (q-lNew) /factorial(q) / (sqrt(delta2+1))^(lNew+1)*D^lNew)+lgamma(((Estar-Ekl)+v0)/2)+lgamma(((Ekr-Estar)+v0)/2)-lgamma(((Ekr-Ekl)+v0)/2)+(((Estar-Ekl)+v0)/2)*log(((gamma0+t(y2) %*% Px2 %*% y2)/2)/((gamma0+t(yL) %*% PxL %*%yL)/2))+(((Ekr-Estar))/2)*log (((gamma0+t(y2) %*% Px2 %*% y2)/2)/((gamma0+t(yR) %*% PxR %*% yR)/2))-(v0/2)*log((gamma0+t(yR) %*% PxR %*% yR)/2)+(v0/2)*log(gamma0/2)-lgamma(v0/2)-log(sum(D^(0:smax)/factorial(0:smax)))

  logR=birth*logR


  if(logR>0){
    res=1
  } else {
    res=exp(logR)
  }
  
#### *D^lNew ?
  return(res)
}


computePostProbNoChangePoints <- function(n, v0, delta2, D, q, l, y, Px) {  
	# Terms that are commented out are constant under static changepoints
	#gamma_term = log(gamma0/2) - lgamma(vo/2) - 
	#  log(sum(D^(0:smax)/factorial(0:smax)));

	#pi_term = - (n/2) * log(2*pi);
	
	q_term = log(factorial(q-l)/factorial(q)/(sqrt(delta2+1))^(l+1)*D^l)
	
	#epsilon_term = lgamma((v0+n)/2)
	
	Px_term = - ((v0 + n)/2) * log((v0 + t(y) %*% Px %*% y)/2)
	
	logProb= q_term + Px_term;
  
  res=exp(logProb)
  return(res)
}

computePx = function(len, x, delta2){
  # INPUT: len, delimiting breakpoints.
  #        x, the observations of X in the corresponding state
  #        delta2.
  # OUTPUT: the projection matrix Px.
  # depends on: . 
  moins = matrix(0,len,len)
  if(prod(dim(x))>0){
    moins = (delta2/(delta2+1))* x%*%pseudoinverse(t(x)%*%x)%*%t(x)
  }
  Px=diag(1,len)-moins
  return(Px)
}


convertModel <- function(binary){
  return(binary %*% 2^(0:(length(binary)-1)))
}


# converts model into binary vector
tobinaire <- function(Model, size){
  bin = array(0, size)
  for(i in (size-1):0){
    bin[i+1] = Model %/% 2^i
    Model = Model %% 2^i
  }
  return(bin)
}


## Convert L in (1..2^q) to S[1:q] (decimal to binary)
LtoS = function(L){
  ###### ATTENTION ######
  # VARIABLES GLOBALES : q
  s = array(1:q)
  for (i in (q-1):0){
    s[q-i] = floor(L / 2^i)
    L = L %% 2^(i)
  }
  return(s[1:q])
}







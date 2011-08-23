## Compute the prior odd ratio for CP postion of TF
priorOdd<-function(pk,total){
  cpt=0
  for (k in 0:(length(pk)-1)){
    cpt=cpt+k/total*pk[k+1]
  }
  return(cpt/(1-cpt))
}
  
## Compute the pior odd ratio for k, the number of changepoint position or the number of TF
kPriorDistribution<-function(kpriorsfile, kmax,alpha,beta){ 
  priors=read.table(kpriorsfile,header=T)
  index1=which(priors[,1]==kmax)
  index2=which(priors[,2]==alpha)
  index3=which(priors[,3]==beta)
  
  pos=which(index1 %in% index2 %in% index3)
  
  if(length(pos)==0){
    stop("The prior probability for these parameters not available here. Please look at the available parameters with function showHyperParms(kmax)")
    pk=NULL
  }else{	 	
    pk=priors[index1[pos],4:(kmax+4)]
  }
  return(pk)
}
   
# kpriorsfile = "k_priors.txt"   
# pkCP=kPriorDistribution(kpriorsfile, kmax=4,alpha=1,beta=0.5)
# pkTF=kPriorDistribution(kpriorsfile, kmax=4,alpha=1,beta=0.5)

HyperParms <- function(alphaCP, betaCP, alphaTF, betaTF, pkCP, pkTF, kpriorsfile, n, q, kmax, smax, dyn, BFOut){

  #####################################################################
  ### rjMCMC hyperparameters 
  #####################################################################

  ### level 1 (4 moves: CP birth, CP death, CP update or phase update)

  ## for birth/death/move/update phase acceptation
  cD = 0.1
  
  
  ### level 2 (4 moves)(for each current phase: Pred birth, Pred death or Regression Coefficient update)
  ## for Pred birth/death acceptation
  c = 0.5

  ## for each hidden state (model selection)
  # for sigma sampling (v0=0,gamma0=0 => p(sig2) ~ 1/sig2 ???????
  # sig2 ~ IG (v0/2,gamma0/2)
  v0 = 1 
  gamma0 = 0.1

  # for the signal-to-noise ratio
  # delta2 ~ IG(alphad2,betad2)
  alphad2 = 2
  betad2 = 0.2 
  
  #######################################################################
  #######################################################################
  
 
  ## for the number of Change Points (CP)
  # for D sampling (D ~ Ga(alphaD,betaD)
  alphaD = alphaCP
  betaD = betaCP

  ## sample delta2, this is not the real one , but for the real one (calculated in sampledelta2Global()) we need data that comes later
  delta2 = rinvgamma(1, shape=alphad2, scale=betad2)

  ## for the number of Transcription Factor (TF)
  ## for lambda ~ Ga(alphalbd,betalbd)
  alphalbd = alphaTF
  betalbd = betaTF
  priorOddCP=NULL
  priorOddTF=NULL
  
  if(is.null(pkCP) && BFOut){
	pkCP=kPriorDistribution(kpriorsfile, kmax,alphaD,betaD)
	if(is.null(pkCP)){
		stop("Please either specify a prior Distribution for CPs number or choose alphaCP and betaCP among the parameters shown by the function showHyperParms(kmax)")
	}
  }
  
  priorOddCP=priorOdd(pkCP,n-1-dyn)
  
  if(is.null(pkTF)&& BFOut){
    pkTF=kPriorDistribution(kpriorsfile, smax,alphalbd,betalbd)
    if(is.null(pkTF)){
	  stop("Please either specify a prior Distrdibution for TFs number or choose alphaCP and betaCP among the parameters shown by the function showHyperParms(kmax)")
	}	
  }
  
  priorOddTF=priorOdd(pkTF,q)
  
  print("")
  print(paste("You chose alphaCP=",alphaCP,", betaCP= ", betaCP, " and alphaTF=",alphaTF,", betaTF= ", betaTF, " which gives the prior dictributions plotted beside") )
  print("")

  HYPERvar = list(cD=cD, alphaD=alphaD, betaD=betaD, c=c, v0=v0, gamma0=gamma0, alphad2=alphad2, betad2=betad2, alphalbd=alphalbd, betalbd=betalbd, delta2=delta2, pkCP=pkCP,pkTF=pkTF,priorOddCP=priorOddCP,priorOddTF=priorOddTF)
  return(HYPERvar)
}

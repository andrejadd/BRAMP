

HyperParms <- function(alphaTF, betaTF){

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
    
  ## sample delta2, this is not the real one , but for the real one (calculated in sampledelta2Global()) we need data that comes later
  delta2 = rinvgamma(1, shape=alphad2, scale=betad2)

  ## for the number of Transcription Factor (TF)
  ## for lambda ~ Ga(alphalbd,betalbd)
  alphalbd = alphaTF
  betalbd = betaTF
    
  cat("[INIT] Priors: aTF: ",alphaTF,", bTF: ", betaTF, "\n") 
  
  HYPERvar = list(cD=cD, c=c, v0=v0, gamma0=gamma0, alphad2=alphad2, betad2=betad2, alphalbd=alphalbd, betalbd=betalbd, delta2=delta2)

  return(HYPERvar)
}

####
####
####
#### main function

main <- function(X, Y, initiation, GLOBvar, HYPERvar){
  
  ### assignement of global variables used here ###
  niter = GLOBvar$niter
  smax = GLOBvar$smax
  
  q = GLOBvar$q
  birth_proposals = GLOBvar$birth_proposals
  ### end assignement ###

  ### assignement of hyperparameters variables used here ###
  cD = HYPERvar$cD
  alphaD = HYPERvar$alphaD
  betaD = HYPERvar$betaD
  ### end assignement ###
  
  ### assignement of initiation variables used here ###
  # initial state
  XE = initiation$initState$XE
  YE = initiation$initState$YE
  S2Dall = initiation$initState$S2Dall
  B2Dall = initiation$initState$B2Dall
  Sig2_2Dall = initiation$initState$Sig2_2Dall
  
  ## counting CP moves (CP Birth, CP death, CP move, Updating phases)
  cptMove2 = array(0,4)
  acceptMove2 = array(0,4)

  ## counting "Updating phases" moves (Edge Birth, Edge death, Udating regression coefficient)
  cptMove = array(0,3)
  acceptMove = array(0,3)
  
  ## everything important stored here and saved at the end to a file
  Structsamples = list(struct = list(), XE = list(), YE = list(), iter=list(), regression.coeff=list())
  counters = list()

  ## when to start saving the regression coefficients (takes up memory)
  start.save.regr = niter - floor(niter * 1/5)
  
  # do main iteration
  r = 1
  while(r < niter) {

    r = r + 1

    ## decide by same chance which axis to alter: x or y 
    ALTERX = T
    k = length(XE) - 2  # need this for next rgamma calculation
    kmax = GLOBvar$kmax.x       # need this for rho calculation (move prob)
    
    if(runif(1,0,1) < 0.5) {
      ALTERX = F
      k = length(YE) - 2
      kmax = GLOBvar$kmax.y
    }
    
    ## mean nr. of changepoints (lambda)    
    D = rgamma(1, shape= k + alphaD, rate = 1+betaD)

    ## the move probabilities
    rho = computeRho4(k, 0, kmax, cD, D)

    ## Sample u to choose one of the 4 moves : CP birth, CP death, CP shift, Update phases.
    u1 = runif(1, 0, 1)
    
    ## Run 1 out of the 4 moves (depending on the value of u)
    if (u1 < rho[1]){
      ## CP birth move: return the new model if the move is accepted, the previous model otherwise.
      out = cp.birth(ALTERX, XE, YE, S2Dall, B2Dall, Sig2_2Dall, X, Y, D, GLOBvar, HYPERvar, F,F)

    } else {
      if(u1 < rho[2]){
     
        ## CP death move: return the new model if the move is accepted, the previous model otherwise.
           out = cp.death(ALTERX, XE, YE, S2Dall, B2Dall, Sig2_2Dall, X, Y, D, GLOBvar, HYPERvar, F,F)

      } else {
        if(u1 < rho[3]){

          ## CP shift move: return the new model if the move is accepted, the previous model otherwise.
          out =  cp.shift(ALTERX, XE, YE, S2Dall, B2Dall, Sig2_2Dall, X, Y, GLOBvar, HYPERvar, F,F)

        } else {
          ## Update phases: return the new model if the move is accepted, the previous model otherwise.
          out = phase.update(XE, YE, S2Dall, B2Dall, Sig2_2Dall, X, Y, GLOBvar, HYPERvar, F, F)
       
          cptMove[out$move1] = cptMove[out$move1]+1
          acceptMove[out$move1] = acceptMove[out$move1]+out$accept1

        }
      }
    }
    
    ## 
    ## Apply changes to the current model
    ##
    XE = out$XE
    YE = out$YE
    B2Dall = out$B2Dall
    S2Dall = out$S2Dall
    Sig2_2Dall = out$Sig2_2Dall

    cptMove2[out$move] = cptMove2[out$move]+1
    acceptMove2[out$move] = acceptMove2[out$move]+out$accept

    ## In case I decide to make a delta2 update in phase.update (just before updating the regression parameters) then I should check if this update
    ## happened and dont do it here.
    ## However, right now I do the delta2 update only here. So updates to the regr. parameters will profit from it in the next call to phase.update
    #cat("FIXME: update delta2 from phase.update (if global delta2 in phase.update will be made)\n")
    
    ##
    ## Update delta2, global
    ##
    HYPERvar$delta2 = sampledelta2Global( X, Y, XE, YE, S2Dall, B2Dall, Sig2_2Dall, GLOBvar, HYPERvar, F)
        
    ##
    ## Save Data
    ## every 10th iteration (to save memory and MCMC thin-out)

    if((r %% 10) == 0) {

      # save data, this is written to the major outfile        
      Structsamples$struct[[length(Structsamples$struct) + 1]] = S2Dall
      Structsamples$XE[[length(Structsamples$XE) + 1]] = XE
      Structsamples$YE[[length(Structsamples$YE) + 1]] = YE
      Structsamples$iter[[length(Structsamples$iter) + 1]] = r

      ## only save in later stage, because needs memory
      if(r > start.save.regr) {
        Structsamples$regression.coeff[[length(Structsamples$regression.coeff) + 1]] = B2Dall
      } else {
        Structsamples$regression.coeff[[length(Structsamples$regression.coeff) + 1]] = 1
      }
      
      counters[[length(counters) + 1]] = list(cptMove2=cptMove2, acceptMove2=acceptMove2, cptMove=cptMove, acceptMove=acceptMove)
 
    }

    #  print status every X and save data to disk iteration 
    if((r %% 500) == 0) {
       cat("\nr: ", r, "\t")
       cat("mem(Structsamples): ", object.size(Structsamples)/1048600)
       
     }
        
  } # end iteration

  writeDataOut(GLOBvar, Structsamples, counters)
  
 
  return(1)
}


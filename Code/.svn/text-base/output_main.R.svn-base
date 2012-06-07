
output <- function(counters, listStock, GLOBvar, HYPERvar, OUTvar){

  ### assignement of global variables used here ###
  target = GLOBvar$target
  niter = GLOBvar$niter
  smax = GLOBvar$smax
  qmax = GLOBvar$qmax
  n = GLOBvar$n
  ### end assignement ###

  ### assignement of results variables used here ###
  # counters
  cptMove = counters$cptMove
  acceptMove = counters$acceptMove
  cptMove2 = counters$cptMove2
  acceptMove2 = counters$acceptMove2
  # listStock
  Estock = listStock$Estock
  Sstock = listStock$Sstock
  Bstock = listStock$Bstock
  Sig2stock = listStock$Sig2stock
  ### end assignement ###

  ### assignement of output variables used here ###
  WriteStock = OUTvar$WriteStock
  outputStockPath = OUTvar$outputStockPath
  outputResPath = OUTvar$outputResPath
  outputFile=OUTvar$outputFile
  Picture = OUTvar$Picture
  Format = OUTvar$Format
  ndeb = OUTvar$ndeb
  nfin = OUTvar$nfin
  ### end assignement ###
  
  ## 
  cexValue=1
  colCP=3
  colTF=4



  if(WriteStock){
    ##proportion of mouvements, all procedure
    tab=round(as.matrix(rbind(cptMove2/sum(cptMove2),acceptMove2/cptMove2)),3)
    colnames(tab)=c("Birth","Death","CP update", "Phases update")
    rownames(tab)=c("Prop Moves", "Prop Accept")
    tab[2,4]=1
    write.table(as.matrix(tab), file = paste(outputStockPath, outputFile,"_cptMoveAll_",target,".txt",sep=""))

    ##prop de mouvement, model upadting
    tab=round(as.matrix(rbind(cptMove/sum(cptMove),acceptMove/cptMove)),3)
    colnames(tab)=c("Birth","Death","Change")
    rownames(tab)=c("Prop Moves", "Prop Accept")
    tab[2,3]=1
    write.table(as.matrix(tab), file = paste(outputStockPath,outputFile,"_cptMovePhases_",target,".txt",sep=""))

    write.table(Estock,paste(outputStockPath,outputFile,"_stockE_",target,".txt",sep=""))
    write.table(Sstock,paste(outputStockPath,outputFile,"_stockS_",target,".txt",sep=""))
	write.table(Bstock[seq(1, dim(Bstock)[1], 10),], paste(outputStockPath,outputFile,"_stockB_",target,".txt",sep=""))
    write.table(Sig2stock,paste(outputStockPath,outputFile,"_stockSig2_",target,".txt",sep=""))
  }
  
  # Reduce matrices if analyses not on all iterations
  Estocki = Estock[ndeb:nfin,]
  Sstocki = Sstock[ndeb:nfin,]
  Bstocki = Bstock[ndeb:nfin,]
  Sig2stocki = Sig2stock[ndeb:nfin,]

  # count # CP per iteration
  NbOfStatesAll = apply(Estock[1:niter,] > 0, 1, sum) - 2
  if(WriteStock) write.table(as.matrix(NbOfStatesAll),paste(outputStockPath,outputFile,"_NBCPiter_",target,".txt",sep=""))
  NbOfStates = apply(as.matrix(Estocki[1:(nfin-ndeb+1),] > 0), 1, sum) - 2


  ## posterior distribution of # of CP
  countNbStates = array(0,smax+1)
  for (i in 1:(smax+1)){
    countNbStates[i] = sum(NbOfStates==(i-1))
  }
 
  ## ajoute par Sophie 01/03/09
  ## Output computed with Bayes Factor
  if(OUTvar$BF){
    dataout =	outputBF(NbOfStatesAll, countNbStates, Estocki, Sstocki, Bstocki, GLOBvar, HYPERvar, OUTvar)
  }
	
  ## ajoute par Sophie 01/03/09
  ## Simple Output 
  if(OUTvar$simple){
    # most represented number of ChangePoints
    NbCP = which(countNbStates==max(countNbStates))-1
    dataout = outputSimple(NbOfStatesAll, countNbStates, NbCP, Estocki, Sstocki, Bstocki, GLOBvar, OUTvar)
    
  } 

 
  # positions of ChangePoints
  bestE = dataout$E
  # matrix (# of phases) x (q): models for each phases 
  bestModels = dataout$bestModels

  print("results exported")
  return(list(bestE=bestE, bestModels=bestModels))
}

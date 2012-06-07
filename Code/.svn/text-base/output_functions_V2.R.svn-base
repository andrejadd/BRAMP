#utilities
isequal <- function(vect, starvect){
  # function that returns 
  out = sum(abs(vect - starvect))==0
  return(out)
}

# function that returns TRUE if phase start->end is in vector E
hasPhase <- function(E, start, end){
  if(!(start %in% E & end %in% E)) return(FALSE)
  else {
    if((end != start+1) & (sum((start+1):(end-1) %in% E)!=0)) return(FALSE)
    else return(TRUE)
  }
}

# function that returns start position in E
PhasePos <- function(E, start){
  posStart = which(E==start)
  return(posStart)
}

#function that displays density of coeffs for a given edge
CoeffDensity <- function(S, edge, Bstock, PhaseRow, PhaseCol, ltype, ph, q){
  rows = PhaseRow[which(S==edge)]
  cols = PhaseCol[which(S==edge)]
  B = NULL
  # extraction of the q+1 coefficients associated to model "edge" in current phase
  for(i in 1:length(rows)) { B = rbind(B, Bstock[rows[i],((cols[i]-1)*(q+1)+1):(cols[i]*(q+1))]) }
  for(TF in 1:(ncol(B)-1)) { lines(density(B[,TF]), col=TF, lty=ltype) }
  return(B)
}

# Simple output (greatest posterior distribution)

outputSimple <- function(NbOfStatesAll,countNbStates,NbCP, Estocki, Sstocki, Bstocki, GLOBvar, OUTvar){
  # function for simple output
  cexValue=1
  colCP=3
  colTF=4
  ### assignement of global variables used here ###
  target = GLOBvar$target
  smax = GLOBvar$smax
  q = GLOBvar$q
  n = GLOBvar$n
  #ajoute Sophie 17/08/09
  dyn = GLOBvar$dyn
  niter = GLOBvar$niter
  bestPosMat = GLOBvar$bestPosMat
  ##ajoute par Sophie 2 juillet09
  predNames=GLOBvar$predNames
  targetNames= GLOBvar$targetNames 
  ### end assignement ###

  ### assignement of output variables used here ###
  outputResPath = OUTvar$outputResPath
  outputFile = OUTvar$outputFile
  Picture = OUTvar$Picture
  Format=OUTvar$Format
  #secondModel retabli par sophie 17/08/09
  secondModel = OUTvar$secondModel
  WriteStock = OUTvar$WriteStock
  ### end assignement ###


  # 1. CHANGEPOINTS
		
    # most represented number of ChangePoints
    NbCP = which(countNbStates==max(countNbStates))-1
	
 if(Picture){
    # name of picture file
    pictFileName = paste(outputResPath, outputFile,"_AllFigures_OutSimple_", target, sep="")
    # create correct graphic
    if(Format != 0){
      jpeg(paste(pictFileName,".jpg",sep=""))
    } else {
      postscript(paste(pictFileName,".eps",sep=""))
    }
    par(mfrow=c(2,2),cex=cexValue, mar=c(3,2,2,1))
    plot(1:niter, NbOfStatesAll, type="l", xlab="Iteration", ylim=c(0,max(NbOfStatesAll)), ylab="Nb of ChangePoints", main=paste("Number of Changepoints accross iterations - Target", target), col=4, lwd=1)
    end = max(which(countNbStates > 0)) + 1
    barplot(countNbStates[1:end], names=0:(end-1), lwd=3,col=colCP, main=paste("Number of CP posterior distribution - Target", target), ylab="Density", xlab="Number of ChangePoints")
  }
  
  configurations=unique(Estocki)
	cpt=array(0,dim(configurations)[1])
	
	for(i in 1:dim(configurations)[1]){
		test=apply(t(t(Estocki)==configurations[i,]), 1, sum)
		cpt[i]=sum(test==(dim(configurations)[2]))
	}
	sum(cpt)
	cpt=cpt/sum(cpt)
	
	configurations=cbind(cpt, configurations)
	configurations[order(cpt,decreasing=T)[1:min(10,dim(configurations)[1])],]
		




  
  # count # of CP at each timepoint
  # as timepoint 1 is always a CP, no need to take it in count, so CPatTP[1]=0
  countCPatTP = array(0, n)
  for(tp in (Estocki[1,1]+1):n) countCPatTP[tp]= sum(Estocki==tp)
  names(countCPatTP)=1:n
  
  if(Picture){
    # posterior proba of CP at each timepoint :
    PPCPatTP = countCPatTP / nrow(Estocki)
    barplot(PPCPatTP, lwd=3, col=colCP, main = paste("CP Position posterior probability - Target",target), ylab="Posterior Probability",xlab="CP position", ylim=c(0,max(PPCPatTP)))
  } # end picture

  #### choose position(s) of CP / define phases of the model
  if(NbCP>0){
    CPpos = as.numeric(names(sort(countCPatTP, decreasing=TRUE)[1:NbCP]))
  } else {
    CPpos = NULL
  }
  Estar = sort(c(1+dyn, CPpos, n+1))

  # Plot Best CP Vector chosen
  if(Picture){
    CPvect = rep(0,smax+2)
    CPvect[Estar]=1
    plot(0, type="n", main="", xlab="", ylab="", bty="n", axes=FALSE, xlim=c(0,10), ylim=c(0,10))
    text(5,7, paste("Identified CP positions: \t", paste(Estar, collapse=" ")), adj=0.5, cex=1)
    #text(5,5, paste("CP vector\t", paste(CPvect, collapse=" ")), adj=0.5, cex=1)
  } # end Picture

  # 2. TRANSCRIPTION FACTORS
  # output models vector (1 model per phase)
  ModelsVector = NULL
  bestModels = NULL

  # for each phase:
  for(ph in 1:(length(Estar)-1)){
    start = Estar[ph]
    end = Estar[ph+1]
    if(Picture){
      plot(0, type="n", main="", xlab="", ylab="", bty="n", axes=FALSE, xlim=c(0,10), ylim=c(0,10))
      text(5,7, paste("Phase\t", ph), adj=0.5, cex=1)
      text(5,5, paste("time",start,"to",end,sep="\t"), adj=0.5, cex=1)
    }
    
    # search for index of iterations where phase is encountered
    PhaseRow = which(apply(Estocki, 1, hasPhase, start, end))
	
	#ajoute Sophie 17/08/09
	if(length(PhaseRow)==0){
		print(paste(" !!!! Problem segment" , ph, " !!!!"))
		bestModels = rbind(bestModels, array(NA, dim(bestModels)[2]))
	}else{
		# positions of start at each iteration
		PhaseCol = apply(Estocki[PhaseRow,], 1, PhasePos, start)
    
		## show density of edges
		# find models of corresponding phase (use of unique index of matrices - by rows)
		S = Sstocki[(PhaseCol-1)*nrow(Sstocki)+PhaseRow]
  
		# biggest code of model :
		if(q < 10){
			EdgeMax = rep(1,q) %*% 2^(0:(q-1))
			# count density of all models
			CountEdges = NULL
			for(edge in 0:EdgeMax) CountEdges = c(CountEdges, sum(S==edge))
			names(CountEdges)=0:EdgeMax
		} else {
			CountEdges = table(S)
		}
		if(Picture){
			# posterior density of models for this phase :
			PDedges = CountEdges / length(PhaseRow)
		barplot(PDedges, lwd=3, col=colTF, main = paste("Models density - Phase",ph , "- Target",target), ylab="Posterior Probability",xlab="models")
		}

		#### choose 1 or 2 best model(s)
		BestEdges = as.numeric(names(sort(CountEdges, decreasing=TRUE)[1:2]))
		#retire par Sophie 17/08/09
		#Qt = sort(CountEdges, decreasing=TRUE)[1:2]
		#if(Qt[2]>=(0.75*Qt[1])){
		#	secondModel=TRUE
		#} else { secondModel=FALSE }
		#secondModel=FALSE
    
		## show density of coeff for best models
		if(Picture){
			plot(0, type="n", main=paste("Coefficients density in phase", ph),  xlim=c(min(Bstocki),max(Bstocki)), ylim=c(0,10), xlab="coeff value", ylab="density")
			B1 = CoeffDensity(S=S, edge=BestEdges[1], Bstock=Bstocki, PhaseRow=PhaseRow, PhaseCol=PhaseCol, ltype=1, ph=ph, q=q)
			lgd = paste(1:q,"model",BestEdges[1])
			lgd.col = 1:q
			lgd.lty = rep(1,q)
			cpt = 1

			# do the same for second best model
			if(secondModel){
				B2 = CoeffDensity(S=S, edge=BestEdges[2], Bstock=Bstocki, PhaseRow=PhaseRow, PhaseCol=PhaseCol, ltype=2, ph=ph, q=q)
				lgd = c(lgd, paste(1:q,"model",BestEdges[2]))
				lgd.col = c(lgd.col, 1:q)
				lgd.lty = c(lgd.lty, rep(2,q))
				cpt = 2
			}
			legend("topright", legend=lgd, col=lgd.col, lty=lgd.lty, bty="n", cex=1)

		## show mean coeff for 2 best models for each TF
		B1mean = apply(abs(B1[,1:q]),2,mean)
		if(secondModel){
			B2mean = apply(abs(B2[,1:q]),2,mean)
			bar.col = NULL
			for(i in 1:q) { bar.col = c(bar.col, rep(i,2)) }
			barplot(rbind(B1mean,B2mean), beside=TRUE, width=1, names.arg=1:q, col=bar.col, density=rep(c(-1,20),2), angle=90, main=paste("TF mean coeff - Target",target, "Phase", ph, "models", paste(BestEdges,collapse=" and ")), ylim=c(0,max(B1,B2)))
		} else {
			barplot(B1mean, width=1, space=5, names.arg=1:q, col=1:q, main=paste("TF mean coeff - Target",target, "Phase", ph, "model", BestEdges[1]), ylim=c(0,max(B1)))
		}
		} # end of Coeff pictures

		# add BestModel to output models
		ModelsVector = c(ModelsVector, BestEdges[1])
		bestModels = rbind(bestModels, tobinaire(BestEdges[1],q))
		
	} # end 'if(length(PhaseRow)==0)'
		
	} # end for phases

    #Output Summary    
    outfile = paste(outputResPath, outputFile, "_OutSimple_Summary_",target,".txt",sep="")
    cat(paste("##### target",target,"#####\n"), file=outfile, append=TRUE)
    cat("##positions of change points: \n", file=outfile, append=TRUE)
    write(Estar, ncolumns=length(Estar), file=outfile, append=TRUE)
    cat("## TF models for each phase: \n", file=outfile, append=TRUE)
    write.table(bestModels, file=outfile, quote=FALSE, col.names=FALSE, row.names=paste("phase", 1:nrow(bestModels),sep="_"),append=TRUE)
	
	#Output Edges list
	outfile = paste(outputResPath, outputFile, "_0_OutSimple_EdgesList.txt",sep="")
	for(i in 1:(length(Estar)-1)){
		for(j in which(bestModels[i,1:q]==1)){
			write.table(matrix(c(j,predNames[j],target,targetNames[target],Estar[i],Estar[i+1]),1,6),file=outfile,append=TRUE,col.names=F,row.names=F)
		}
	}

	# Output CP list
	if(NbCP>0){	
		outfile = paste(outputResPath, outputFile, "_0_OutSimple_CPlist.txt",sep="")
		write.table(as.matrix(cbind(target,targetNames[target],Estar[2:(NbCP+1)])),file=outfile,append=TRUE,col.names=F,row.names=F)
	}
	
 
 
	# close graphic device
	if(Picture) dev.off()
  
  return(list(E=Estar, bestModels=bestModels))
}
  
 
# BF output (model with the greatest BF is selected) 
outputBF<-function(NbOfStatesAll, countNbStates, Estocki, Sstocki, Bstocki, GLOBvar, HYPERvar, OUTvar){
	# function for output with bayes factor
	cexValue=1
	colCP=3
	colTF=4
	### assignement of global variables used here ###
	target = GLOBvar$target
	smax = GLOBvar$smax
	qmax = GLOBvar$qmax
	q = GLOBvar$q
	n = GLOBvar$n
	  #ajoute Sophie 17/08/09
	dyn = GLOBvar$dyn
	niter = GLOBvar$niter
	bestPosMat = GLOBvar$bestPosMat
	pkTF=HYPERvar$pkTF
	##ajoute par Sophie 2 juillet09
	predNames=GLOBvar$predNames
	targetNames= GLOBvar$targetNames 
	### end assignement ###

	### assignement of output variables used here ###
	outputResPath = OUTvar$outputResPath
	outputFile=OUTvar$outputFile
	Picture = OUTvar$Picture
	Format=OUTvar$Format
	ndeb = OUTvar$ndeb
	nfin = OUTvar$nfin
	WriteStock = OUTvar$WriteStock

	

 if(Picture){
    # name of picture file
    pictFileName = paste(outputResPath, outputFile,"_AllFigures_OutBF_", target, sep="")
    # create correct graphic
    if(Format != 0){
      jpeg(paste(pictFileName,".jpg",sep=""))
    } else {
      postscript(paste(pictFileName,".eps",sep=""))
    }
    par(mfrow=c(2,2),cex=cexValue, mar=c(3,2,2,1))
    plot(1:niter, NbOfStatesAll, type="l", ylim=c(0,max(NbOfStatesAll)), xlab="Iteration", ylab="Nb of ChangePoints", main=paste("Number of Changepoints accross iterations - Target", target), col=4, lwd=1)

    end = max(which(countNbStates > 0)) + 1
    barplot(countNbStates[1:end], names=0:(end-1), lwd=3,col=colCP, main=paste("Number of CP posterior distribution - Target", target), ylab="Density", xlab="Number of ChangePoints")
  }



        
        
	# Bayes Factor for the # of CP
	PPnbCP=countNbStates/(nfin-ndeb+1)
	PPnbCP=PPnbCP[1:length(HYPERvar$pkCP)]
	names(PPnbCP)=0:smax
	##### save DATA ####
	outfile = paste(outputResPath, outputFile, "_OutBFData_",target,".txt",sep="")
        cat(paste("\n### Posterior distribution for the number of CP for target",target, "####\n"), file=outfile, append=TRUE)
        write.table(matrix(PPnbCP, nrow=1), sep="\t", quote=FALSE, row.names=FALSE, col.names=0:smax, file=outfile, append=TRUE)

	BFnbCP=PPnbCP/(1-PPnbCP)/HYPERvar$pkCP/(1-HYPERvar$pkCP)
	##### save DATA ####
	cat(paste("\n### Bayes Factor for the number of CP for target",target, "####\n"), file=outfile, append=TRUE)
        write.table(matrix(BFnbCP, nrow=1), sep="\t", quote=FALSE, row.names=FALSE, col.names=0:smax, file=outfile, append=TRUE)

	if(Picture){
		barplot(as.numeric(HYPERvar$pkCP/(1-HYPERvar$pkCP)), names=0:smax, lwd=3,col=colCP, main=paste("Prior Odd for the number of CP - Target", target),sub =paste("(Prior odd for each CP position = ",round(as.numeric(HYPERvar$priorOddCP),3),")" ), ylab="Prior Odd", xlab="Number of ChangePoints")
		barplot(as.numeric(pkTF/(1-pkTF)), names=0:qmax, lwd=3,col=colTF, main=paste("Prior Odd for the number of TF - Target", target), sub =paste("(Prior odd for each TF = ",round(as.numeric(HYPERvar$priorOddTF),3),")" ), ylab="Prior Odd", xlab="Number of TFs")

	}



	if(Picture){
		#end = max(which(countNbStates > 0)) + 1
		barplot(as.numeric(BFnbCP), names=0:(length(BFnbCP)-1), lwd=3,col=colCP, main=paste("Bayes Factor for the number of CP - Target", target), ylab="Bayes Factor", xlab="Number of ChangePoints")
		#barplot((BFnbCP), names=0:(length(BFnbCP)-1), lwd=3,col=colCP, main=paste("Bayes Factor for the number of CP - Target", target), ylab="Bayes Factor", xlab="Number of ChangePoints")

	}

	## Number of Changepoints :
	NbCP=which(BFnbCP==max(BFnbCP))-1




	# 1. CHANGEPOINTS
	# Posterior Probability for a CP at each TimePoint
	# as timepoint 1 is always a CP, no need to take it in count, so CPatTP[1]=0
	PPCPatTP = array(0, n)
	for(tp in (Estocki[1,1]+1):n) PPCPatTP[tp]= sum(Estocki==tp)
	PPCPatTP=PPCPatTP/(dim(Estocki)[1])
	##### save DATA ####
	cat(paste("\n### Posterior distrbution for CP position for target",target, "####\n"), file=outfile, append=TRUE)
    write.table(matrix(PPCPatTP, nrow=1), sep="\t", quote=FALSE, row.names=FALSE, col.names=1:n, file=outfile, append=TRUE)
	

	
	CPBF= 	PPCPatTP/(1-PPCPatTP)/as.numeric(HYPERvar$priorOddCP)
	names(CPBF)=1:n
	##### save DATA ####
	cat(paste("\n### Bayes Factor  forCP position for target",target, "####\n"), file=outfile, append=TRUE)
    write.table(matrix(CPBF, nrow=1), sep="\t", quote=FALSE, row.names=FALSE, col.names=1:n, file=outfile, append=TRUE)
	
		
	if(Picture){
		barplot(CPBF, lwd=3, col=colCP, main = paste("CP Position Bayes Factor - Target",target), ylab="Bayes Factor",xlab="CP position", ylim=c(0,max(CPBF)))
		## ?## lines (y=3)
	} # end picture

	#### choose position(s) of CP / define phases of the model
	if(NbCP>0){
		CPpos = sort(as.numeric(names(sort(CPBF, decreasing=TRUE)[1:NbCP])))
	} else {
		CPpos = NULL
	}
	
	## modifie par Sophie 01/03/09
	## Estar commence a 1+dyn
	Estar = sort(c(1+dyn, CPpos, n+1))
	
	
	##################### 	## modifie par Sophie 01/04/09
	####  NEW : in case there is 2 or more changepoints next to each other or very close...
	#####################
	if(NbCP>1){
			pb=which(Estar[2:(NbCP+1)]-Estar[1:(NbCP)]<=5)
			while(length(pb)>0){
				l=pb[1]
				toRemove=l-1+order(CPBF[Estar[c(l,l+1)]])[1]
				Estar=Estar[-c(toRemove)]
				NbCP=NbCP-1
				pb=which(Estar[2:(NbCP+1)]-Estar[1:(NbCP)]<=5)
		
			}
			
		}

	#####################
	#####################
	
	
  # Plot Best CP Vector chosen
#  if(Picture){
#    CPvect = rep(0,smax+2)
#    CPvect[Estar]=1
#    plot(0, type="n", main="", xlab="", ylab="", bty="n", axes=FALSE, xlim=c(0,10), ylim=c(0,10))
#    text(5,7, paste("Estar\t", paste(Estar, collapse=" ")), adj=0.5, cex=1)
#    text(5,5, paste("CP vector\t", paste(CPvect, collapse=" ")), adj=0.5, cex=1)
#  } # end Picture


 ###### 2. TRANSCRIPTION FACTORS ######
  # for all model edges storage
  bestModels = NULL 
  
  # for each phase
  for (l in 1:(length(Estar)-1)){
    start = Estar[l]
    end = Estar[l+1]
    # search for index of iterations where phase is encountered
    PhaseRow = which(apply(Estocki, 1, hasPhase, start, end))
	
	if(length(PhaseRow)==0){
		print(paste(" !!!! Problem segment" ,l, " !!!!"))
		bestModels = rbind(bestModels, array(NA, q*2))
	}else{
		# positions of start at each iteration
		PhaseCol = apply(Estocki[PhaseRow,], 1, PhasePos, start)
		
		# find TF coefficients of the corresponding phase
		Bphase = matrix(0, length(PhaseRow),q+1)
		for(i in 1:length(PhaseRow)){
			row = PhaseRow[i]
			col = (q+1)*(PhaseCol[i]-1)+1
			Bphase[i,] = Bstocki[row, col:(col+q)]
		}
   
		# posterior probability of the number of TFs :
		NbOfTFs=apply(abs(Bphase[,1:q])>0,1,sum)
		# posterior distribution of # of CP
		PPnbTF = array(0,qmax+1)
		for (i in 1:(qmax+1)){
			PPnbTF[i] = sum(NbOfTFs==(i-1))
		}
		PPnbTF=PPnbTF/(dim(Bphase)[1]) ## 
		##### save DATA ####
		cat(paste("\n### Posterior distribution for the number of TF for target ",target," in phase ",l, " ####\n"), file=outfile, append=TRUE)
		write.table(matrix(PPnbTF, nrow=1), sep="\t", quote=FALSE, row.names=FALSE, col.names=0:qmax, file=outfile, append=TRUE)

		#if(Picture){
		#	barplot(PPnbTF, names=0:(length(PPnbTF)-1), lwd=3,col=colCP, main=paste("Posterior Distribution for the number of TF - Target", target), ylab="Bayes Factor", xlab="Number of TFs")
		#}

		BFnbTF=PPnbTF/(1-PPnbTF)/pkTF/(1-pkTF)
		##### save DATA ####
		## Number of TFs :
		NbTF=which(BFnbTF==max(BFnbTF))-1
		BFnbTF[which(BFnbTF==Inf)]=100000
		##### save DATA ####
		cat(paste("\n### Bayes Factor for the number of TF for target ",target," in phase ",l, " ####\n"), file=outfile, append=TRUE)
		write.table(matrix(BFnbTF, nrow=1), sep="\t", quote=FALSE, row.names=FALSE, col.names=0:qmax, file=outfile, append=TRUE)

		if(Picture){
			barplot(as.numeric(BFnbTF), names=0:(length(PPnbTF)-1), lwd=3,col=colCP, main=paste("Bayes Factor for the number of TF - Target", target), ylab="Bayes Factor", xlab="Number of TFs")
		}

		# posterior probability of each TF :
		PPTF = apply(abs(Bphase[,1:q]) > 0, 2, sum) / (length(PhaseRow))  
		names(PPTF)=1:q
		##### save DATA ####
		cat(paste("\n### posterior Probabilities of TF for target",target, "in phase",l,"####\n"), file=outfile, append=TRUE)
		write.table(matrix(PPTF,nrow=1),quote=FALSE, row.names=FALSE, file=outfile, col.names=1:q, append=TRUE)

	
		# calculate bayes factor 
		BFTF = PPTF / (1-PPTF) / as.numeric(HYPERvar$priorOddTF)
		BFTF[which(BFTF==Inf)]=100000
		##### save DATA ####
		cat(paste("\n### Bayes Factor of CP position for target",target, "####\n"), file=outfile, append=TRUE)
		write.table(matrix(BFTF, nrow=1), sep="\t", quote=FALSE, row.names=FALSE, col.names=1:q, file=outfile, append=TRUE)
	
		VectorTF=(BFTF > sort(BFTF,decreasing=T)[NbTF+1])*1
	
		# vector TF edges 1/0 (1 if coeff > 0)
		Sphase = (abs(Bphase[,1:q])>0)*1
		posVectorTF = which(apply(Sphase , 1, isequal, VectorTF))
		# mean coeff for each TF 
		#coef = apply(Bphase[posVectorTF,], 2, mean)
		##### save DATA ####
		#cat(paste("\n### Best Model: Mean Coeff of TF for target",target, "in phase",l,"####\n"), file=outfile, append=TRUE)
		#write.table(matrix(coef, nrow=1), row.names=FALSE, quote=FALSE, file=outfile, append=TRUE)
		#if (WriteStock){
		#  cat(paste("\n### Best Model: Coeff of TF for target",target, "in phase",l,"####\n"), file=paste(outputStockPath, outputFile,"_stockCoef",target,".txt",sep=""), append=TRUE)
		#  write.table(Bphase[posVectorTF,], file=paste(outputStockPath,outputFile,"_stockCoef",target,".txt",sep=""), quote=FALSE, row.names=FALSE, append=TRUE)
		#}

       
		if(Picture){
			# case Bayes Factor
			maintitle = paste("TF Bayes Factor - Target ", target, " Phase ", l, " (nb=", length(PhaseRow), sep="")
			maintitle = paste(maintitle, ")", sep="")
			barplot(BFTF, width=1, space=5, names.arg=bestPosMat[target,], main=maintitle, ylab = "Bayes Factor", col=1:q, ylim=c(0,max(BFTF, OUTvar$BFthres[2])))
			abline(h =OUTvar$BFthres[2], lty=3)  
			abline(h =OUTvar$BFthres[1], lty=3, col="grey")
		} # end Picture

		# store best model
		bestModels = rbind(bestModels, c(VectorTF,round(BFTF,1)))
	} # end if length(PhaseRow)=0
		
	} # end for each Phase

    #Output Summary
	# CP
	outfile = paste(outputResPath, outputFile, "_OutBF_SummaryCP_",target,".txt",sep="")
	write(as.array(CPBF), ncolumns=length(CPBF), file=outfile, append=F)
	CP=array(0,length(CPBF))
	CP[1:length(Estar)]=Estar
    write(CP, ncolumns=length(CPBF), file=outfile, append=TRUE)

	# TF
	outfile = paste(outputResPath, outputFile, "_OutBF_SummaryTF_",target,".txt",sep="")
	#names(bestModels)=c(paste("TF", 1:q, sep="_"),paste("BF", 1:q, sep="_"))
	write.table(bestModels, file=outfile, quote=FALSE, col.names=c(paste("TF", 1:q, sep="_"),paste("BF", 1:q, sep="_")), row.names=paste("phase", 1:nrow(bestModels),sep="_"),append=F)
 	
	
    #Output Edges list
	outfile = paste(outputResPath, outputFile, "_0_OutBF_EdgesList.txt",sep="")
	for(i in 1:(length(Estar)-1)){
		for(j in which(bestModels[i,1:q]==1)){
			write.table(matrix(c(j,predNames[j],target,targetNames[target],Estar[i],Estar[i+1],round(bestModels[i,j+q],1)),1,7),file=outfile,append=TRUE,col.names=F,row.names=F)
		}
	}

	# Output CP list
	if(NbCP>0){	
		outfile = paste(outputResPath, outputFile, "_0_OutBF_CPlist.txt",sep="")
		write.table(as.matrix(cbind(target,targetNames[target],Estar[2:(NbCP+1)],round(CPBF[Estar[2:(NbCP+1)]],1),round(max(BFnbCP),1))),file=outfile,append=TRUE,col.names=F,row.names=F)
	}
	
  
	# close graphic device
	if(Picture) dev.off()

   	if(Picture){
		postscript(paste(outputResPath, outputFile,"_CPBF_Figures_OutBF_", target,".eps", sep=""))
		barplot(CPBF, lwd=3, col=colCP, main = paste("CP Position Bayes Factor - Target",target), ylab="Bayes Factor",xlab="CP position", ylim=c(0,max(CPBF)))
		dev.off()
	} # end picture

  return(list(E=Estar, bestModels=bestModels))


  }


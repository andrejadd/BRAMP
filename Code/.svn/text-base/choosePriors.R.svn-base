choosePriors<-function(kmax,priorsPath=paste(codePath,"k_priors.txt",sep="")){
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

#Exemple
#choosePriors(4)
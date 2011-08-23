

# column 1: iteration , column 2: psrf value

iterma = matrix(nrow=0,ncol=2)

for(node in 1:30) {

  input = paste("psrf_n21_i", node, ".data", sep="")

  if(!file.exists(input)) {
    cat("file ", input, " does not exists, skipping..\n")
    next
  }
  
  psrfm = as.matrix(read.table(file=input, header=FALSE))

  iters = psrfm[,1]

  iters = sort(unique(iters))


  for(r in iters) {
		
	
    psrfv = psrfm[which(psrfm[,1] == r),]
	
	
    if(length(psrfv) >= 1) {

      
                                        # this checks if only one sample was returned
      if(is.null(dim(psrfv))) {
        if(psrfv[2] < 1.05) {
          cat("found below th 1\n")
          iterma = rbind(iterma, c(r,node))
        }	
        
        
      } else { # this is the case for many samples
        
        verysmall = psrfv[which(psrfv[,2] <= 1.05),1]
        
        if(length(verysmall) >= 1) {		
          
          cat("found below th - length(verysmall):", length(verysmall), "\n")
          
          itertmp = matrix(node, nrow=length(verysmall), ncol=2)
          
          itertmp[1:length(verysmall),1] = t(verysmall)
          
          print.table(itertmp)
          
          iterma = rbind(iterma, itertmp)
          
          
        }
        
      }
    }
    
  }
}

print.table(iterma)


boxplot(iterma[,1] ~ iterma[,2], xlab="Nodes", ylab="Iteration")


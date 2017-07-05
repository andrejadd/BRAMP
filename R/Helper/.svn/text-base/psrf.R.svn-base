# Calculate the potential scale reduction factor (by Rubin and Gelman)
#
# Input:
# parameters - A list of MCMC trajectories, where each trajectory is a matrix
# with nbPars rows and nbIterations columns, where nbPars is the number of
# parameters and nbIterations is the number of samples. 
psrf <- function(parameters) {
  
  # number of seqences
  nbSeq=length(parameters)
  #cat("[psrf()] nbSeq: ", nbSeq, "\n")

  nbIterations=dim(parameters[[1]])[2]
  #cat("[psrf()] nbIterations: ", nbIterations, "\n")

  nbPars = dim(parameters[[1]])[1]
  #cat("[psrf()] nbPars: ", nbPars, "\n")
  
  # compute B
  seq_means = matrix(0, nbPars, nbSeq)
  
  for(i in 1:nbSeq) {
   seq_means[, i] = apply(parameters[[i]], 1, mean)
  }

  B = nbIterations / (nbSeq-1) * apply((seq_means-matrix(apply(seq_means,1,mean) , nbPars, nbSeq) )^2 , 1 , sum)
  
  #par(mfrow=c(1,1),cex=1.5)
  #barplot(B, names.arg=1:nbPars,main="Between-sequence variance of the changepoint posterior probability", xlab="Timepoint", ylab="Potential reduction scale factor" )
  
 
  # compute W
  diffs = matrix(0, nbPars, nbSeq)
  
  for(i in 1:nbSeq) {
    # calculate the within segment variances s^2
    diffs[, i] = apply((parameters[[i]] - kronecker(seq_means[, i], matrix(1, 1 , nbIterations)))^2, 1, sum)
  }
  
  W = apply(diffs, 1, sum) / (nbSeq*(nbIterations-1))
  
  #par(mfrow=c(1,1),cex=1.5)
  #barplot(W, names.arg=1:nbPars,main="Within-sequence variance of the changepoint posterior probability", xlab="Timepoint", ylab="Potential reduction scale factor" )
  
  seq_overest = (nbSeq + 1) / nbSeq;
  #cat("[psrf()] seq_overest ", seq_overest, "\n")

  # this is probably not necessary, 
  #PSRF = array(0,nbPars)
  #if(B == 0) {
  #  PSRF = seq_overest * ((nbIterations-1)/nbIterations) - (nbIterations - 1)/(nbIterations*nbSeq)
  #} else {
  #  PSRF = seq_overest * ((nbIterations-1)/nbIterations+ B/(W * nbIterations)) -  (nbIterations - 1)/(nbIterations*nbSeq)
  #}


  PSRF = array(0,nbPars)
  
  # look if there is a B for one of the parameters
  for(i in 1:nbPars) {
    if(W[i] == 0.0) {
      PSRF[i] = seq_overest * ((nbIterations-1)/nbIterations) -  (nbIterations - 1)/(nbIterations*nbSeq)
    } else {
      PSRF[i] = seq_overest * ((nbIterations-1)/nbIterations+ B[i]/(W[i] * nbIterations)) -  (nbIterations - 1)/(nbIterations*nbSeq)
    }
  }

  #par(mfrow=c(1,1),cex=1.5)
  #barplot(PSRF2, names.arg=1:nbPars, main="PSRF of the changepoint posterior probability", xlab="Timepoint", ylab="Potential reduction scale factor" )
 
  return(PSRF)
}

# AA:ADDED FCT
# Call this function to calculate the Potential Scale Reduction Factors for a set of parameters in each iteration for several sequences
# Note, r must be passed because Sstock is none growing (iteration can not be derived)
# Input: Sstock
calcPSRFWrapper <- function(Sstock, nrwindows=5, fixedWndSize=NULL ) {

  r = length(Sstock)
  q = length(Sstock[[1]])
  
  start_it = 1
  fixedWndSize = as.integer(r/nrwindows)

#  cat("==start_it: ", start_it, ", r: ", r, ", fixedWndSize: " , fixedWndSize , "\n")

  
  # each element in this list is a matrix where the colums are iterations and the rows correspond to each node, each node shows
  # if there is an edge to our target (1) or not (0)
  matlist = list()
  
  # a matrix with the number of nodes rows (q) and no column yet
  edgematrix = matrix(0,q,0)

  seq_index = fixedWndSize

  for(iter in 1:(r-1)) {

    seq_index = seq_index - 1;
    
    # this vector holds the edges of the current iteration i
#    edgevector = integer.base.b(x=Sstock[i,1] , ndigits=q)

    edgevector = Sstock[[iter]]
            
    # add to matrix column-wise
    edgematrix = cbind(edgematrix, edgevector)

    # if this is true, it means the current sequence was all read in
    if(seq_index == 0) {
      # DEBUG
      #cat("seq finished on iteration i: ", i, "\n")
      #cat("edgematrix dims: n=" , dim(edgematrix)[1] , ", m=", dim(edgematrix)[2], "\n")
      
      # reset for the new sequence
      seq_index = fixedWndSize
            
      # copy current edgematrix to the list
      matlist[[length(matlist)+1]] = edgematrix
      
      # create new matrix
      edgematrix = matrix(0,q,0)
    }
    
  }


  
  for(w in 1:nrwindows) {
    
    cat("\nmatlist[[", w, "]] - dim:\n")
    print.table(dim(matlist[[w]]))
    
    for(s in 1:dim(matlist[[w]])[2] ) {

      print.table( matlist[[w]][,s] )

    }

  }




  
 # this first gives an array with the PSRF for each parameter (edge) - from psrf() and then  we are interested in the worsed value (max)
  return(max(psrf(matlist)))
}


calcPSRFWrapperSPEED <- function(Sstock, nrwindows=5, maxSamplelen=NULL, RUNINFOS=FALSE) {


  r = length(Sstock)
  q = length(Sstock[[1]])

   if(!is.null(maxSamplelen)) {
     maxIterations = nrwindows * maxSamplelen

     if( (r - maxIterations) > 0) {
       start_it = r - maxIterations
       fixedWndSize = maxSamplelen
     } else {
       start_it = 1
       fixedWndSize = as.integer(r/nrwindows)
     }
   } else {
     start_it = 1
     fixedWndSize = as.integer(r/nrwindows)
   }
  
#  cat("start_it: ", start_it, ", r: ", r, ", q: ", q, ", nrwindows: ", nrwindows, ", maxSamplelen: ", maxSamplelen, ", fixedWndSize: " , fixedWndSize , "\n")

  
  # each element in this list is a matrix where the colums are iterations and the rows correspond to each node, each node shows
  # if there is an edge to our target (1) or not (0)
  matlist = list()

  # a matrix with the number of nodes rows (q) and no column yet
  #cat("nr. possible parents q: ", q, "\n")
  #edgematrix = matrix(0,q,0)
  edgematrix = matrix(0,q,fixedWndSize)

  seq_index = 1
  
  # construct a matrix list, each matrix holds a sequence of networks
  for(i in start_it:r) {

   
    
    # extract edge vector of this iteration of phase 1
    edgevector = Sstock[[i]]#[1,]
    
    edgematrix[1:q,seq_index] = t(edgevector)

    #cat("\nedgematrix in iter ", i, " and seq_index ", seq_index, "\n")
    #print.table(edgematrix)
    
    seq_index = seq_index + 1
  

                                        # if this is true, it means the current sequence was all read in
    if(seq_index == fixedWndSize+1) {
      
      # DEBUG
      #cat("seq finished on iteration i: ", i, "\n")
      #cat("edgematrix dims: n=" , dim(edgematrix)[1] , ", m=", dim(edgematrix)[2], "\n")
      
      # reset for the new sequence
      seq_index = 1#fixedWndSize
            
      # copy current edgematrix to the list
      matlist[[length(matlist)+1]] = edgematrix
      
      # reset edge matrix
      edgematrix = matrix(0,q,fixedWndSize)
    }
    
  }

    
  startt = proc.time()[3]
  maxpsrf = max(psrf(matlist))

  if(RUNINFOS) {
   
    #cat("\tpsrf() took ", (proc.time()[3] - startt), "Sec.\n")
    
    
    memsum = sum(sapply(ls(), fcta <- function(i){
      objs = object.size(get(i))/1048600
                                        # report objects greater than 10MB
      if( objs > 1) { cat("\tpsrf() ", i, " : ", objs, " MB\n" ) }
      return(objs)
    } ))
    
    #cat("\tpsrf() TotalMem : ", memsum, " MB\n")

                                        # explicitely remove this list
  }
  
  rm(matlist)
  
  return(maxpsrf)

}

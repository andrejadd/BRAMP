

timeout <- function(expr, seconds = 60)  {   ## used to wait for a keystroke, use with ' z <- try(silent=TRUE, timeout(readline(prompt="Hit me: "), seconds=5))'
         # Set up a background process that will send a signal
         # to the current R process after 'seconds' seconds.
         # Evaluate expr with an interrupt handler installed
         # to catch the interrupt.
         # If expr finishes before that time it will kill the killer.
  
         killer.pid <- system(intern = TRUE, paste(" (sleep", seconds, " ; kill -INT", Sys.getpid(), ")>/dev/null&\n echo $!"))
         on.exit(system(paste("kill", killer.pid, "> /dev/null 2>&1")))
         withCallingHandlers(expr, interrupt=function(...)stop("Timed  out", call.=FALSE))
 }



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

## AA: not used, but maybe later for defining priors in text file
##
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



##
## The acceptance probability of doing a cut shift. This is easy since we do not change the nr. of segments or nr. of edges, only
## the membership of the locations to one of the two segments.
##
## The prior and proposal are also invariant and are terminated in the posterior ratio but check this again ! (FIXME)
##
##

shift.computeAlpha <- function(X, Y, old.set, proposed.set, seg.ids, HYPERvar, DEBUG_BIRTH_EXT=0) {

  edge.struct = old.set$edge.struct
  
  ##
  ## Posterior of old set (no shift) - Current State
  ##
  sumPhi = 0

  ## get first segment
  x.child1 = extractData(old.set, X, seg.ids[1])
  y.child1 = extractData(old.set, Y, seg.ids[1])

tryCatch( {
  Proj.child1 = computeProjection(as.matrix(x.child1[,which(edge.struct == 1)]), HYPERvar$delta2)
}, error=function(e) { 
   rnd.id = ceiling(runif(1,1,100000))
   save(file=paste("debug.out.", rnd.id, sep=""), x.child1, seg.ids, proposed.set, old.set, X, Y, edge.struct) 
   print(e)
   stop("encountered Error wrote to debug file id: ", rnd.id, "\n")
} 
)
  omega.child1 = length(y.child1)
  
  ## get second segment
  x.child2 = extractData(old.set, X, seg.ids[2])
  y.child2 = extractData(old.set, Y, seg.ids[2])
  Proj.child2 = computeProjection(as.matrix(x.child2[,which(edge.struct == 1)]), HYPERvar$delta2)
  omega.child2 = length(y.child2)

  tryCatch({
    sumPhi  = lgamma((HYPERvar$v0 + omega.child1) / 2) + (-(HYPERvar$v0 + omega.child1) / 2) * log( (HYPERvar$gamma0 + t(y.child1) %*% Proj.child1 %*% y.child1)/2) +  lgamma((HYPERvar$v0 + omega.child2) / 2) + (-(HYPERvar$v0 + omega.child2) / 2) * log( (HYPERvar$gamma0 + t(y.child2) %*% Proj.child2 %*% y.child2) / 2)
  }, error = function(e) {
    cat("Caught error \n ")
    print(e)
  })


  ##
  ## Posterior of proposed set (with shift) - Next State
  ##
  sumPhiPlus = 0

  ## get first segment of proposed set
  x.child1 = extractData(proposed.set, X, seg.ids[1])
  y.child1 = extractData(proposed.set, Y, seg.ids[1])
 
tryCatch( {
 Proj.child1 = computeProjection(as.matrix(x.child1[,which(edge.struct == 1)]), HYPERvar$delta2)
}, error=function(e) { 
   rnd.id = ceiling(runif(1,1,100000))
   save(file=paste("debug.out.", rnd.id, sep=""), x.child1, seg.ids, old.set, proposed.set, X, Y, edge.struct) 
   print(e)
   stop("encountered Error wrote to debug file id: ", rnd.id, "\n")
} 
)
 
  omega.child1 = length(y.child1)

  
  ## get second segment
  x.child2 = extractData(proposed.set, X, seg.ids[2])
  y.child2 = extractData(proposed.set, Y, seg.ids[2])
  Proj.child2 = computeProjection(as.matrix(x.child2[,which(edge.struct == 1)]), HYPERvar$delta2)
  omega.child2 = length(y.child2)

  tryCatch({
    sumPhiPlus  = lgamma((HYPERvar$v0 + omega.child1) / 2) + (-(HYPERvar$v0 + omega.child1) / 2) * log( (HYPERvar$gamma0 + t(y.child1) %*% Proj.child1 %*% y.child1)/2) +  lgamma((HYPERvar$v0 + omega.child2) / 2) + (-(HYPERvar$v0 + omega.child2) / 2) * log( (HYPERvar$gamma0 + t(y.child2) %*% Proj.child2 %*% y.child2) / 2)
  }, error = function(e) {
    cat("Caught error \n ")
    print(e)
  })


  ##
  ## compute acceptance chance
  ##

  alpha = min(1, exp( sumPhiPlus - sumPhi ))
 
  return(alpha)

}

###########################################################################################################################################
##
## Compute acceptance of segment split (cut), returns the log of the acceptance probability
##            Take inverse of this to get the probability of a merge (cut remove)
##
## Why I'm not passing Node pointers instead of segment IDs and Budgets? Because child node pointer, e.g., not known yet for cut move
##
############################################################################################################################################

cp.computeAlpha <- function(HYPERvar, X, Y, old.set, proposed.set, parent.id, child.ids, parent.budget, child.budget, DEBUG_BIRTH_EXT=0) {

  edge.struct = old.set$edge.struct

  ## does not matter from which set, stays the same
  nr.edges = sum(edge.struct) - old.set$additional.parents

  ###### Posterior Prob. of Current state - this is the state of a single segment (not split yet or merged)
  sumPhi = 0

#  if(DEBUG_BIRTH_EXT == 2) {  cat("  [computeAlpha] IDs: ", child.ids[1], " - ", child.ids[2]) }
      
  ## get the full data, i.e. all predictors of the parent
  x = extractData(old.set, X, parent.id)
  y = extractData(old.set, Y, parent.id)
  
  ## Parent: calculate projection matrix
  Pr = computeProjection(as.matrix(x[,which(edge.struct == 1)]), HYPERvar$delta2)
    
  ## Parent: number of locations 
  omega = length(y)

  ## original equation  without log transform (makes it necessary to multiply)
  ## prodPhi = prodPhi * gamma((v0+omega)/2) * ((gamma0+ t(y) %*% Pr %*% y)/2)^(-(v0+omega)/2)

  if( (dim(Pr)[1] != length(y)) || (dim(Pr)[2] != length(y))) {
    stop("dimension mismatch for projected target data in cp.computeAlpha")
  }

  tryCatch({
    sumPhi  = lgamma((HYPERvar$v0 + omega)/2) + (-(HYPERvar$v0 + omega)/2) * log( (HYPERvar$gamma0 + t(y) %*% Pr %*% y)/2)
  }, error = function(e) {
    stop("Caught error \n ")
    print(e)
  })


  ###### Posterior Prob two segments (next state for split, current state for merge)
  ##
  sumPhiPlus = 0


  ## get first child (which matches the parent id)
  x.child1 = extractData(proposed.set, X, child.ids[1])
  y.child1 = extractData(proposed.set, Y, child.ids[1])

  Proj.child1 = computeProjection(as.matrix(x.child1[,which(edge.struct == 1)]), HYPERvar$delta2)
  omega.child1 = length(y.child1)

  
  ## get second child
  x.child2 = extractData(proposed.set, X, child.ids[2])
  y.child2 = extractData(proposed.set, Y, child.ids[2])

  Proj.child2 = computeProjection(as.matrix(x.child2[,which(edge.struct == 1)]), HYPERvar$delta2)
  omega.child2 = length(y.child2)

  tryCatch({
    sumPhiPlus  = lgamma((HYPERvar$v0 + omega.child1) / 2) + (-(HYPERvar$v0 + omega.child1) / 2) * log( (HYPERvar$gamma0 + t(y.child1) %*% Proj.child1 %*% y.child1)/2) +  lgamma((HYPERvar$v0 + omega.child2) / 2) + (-(HYPERvar$v0 + omega.child2) / 2) * log( (HYPERvar$gamma0 + t(y.child2) %*% Proj.child2 %*% y.child2) / 2)
  }, error = function(e) {
    cat("Caught error \n ")
    print(e)
  })

  ## get the dimension for the halfperimeter, Ending Norm means it is scaled into the [0,1] interval
  parent.dim = getSegmentDimScale(old.set, parent.id)
  
  ## get the child dimensions
  child1.dim = getSegmentDimScale(proposed.set, child.ids[1])  
  child2.dim = getSegmentDimScale(proposed.set, child.ids[2])  

  ## log prior ratio
  log.prior.ratio = log(exp(-1 * sum(child1.dim) * child.budget)) + log(exp(-1 * sum(child2.dim) * child.budget)) - log(exp(-1 * sum(parent.dim) * parent.budget)) 

  #cat("prior ratio: " , exp(log.prior.ratio), "\n")
  
  ## Proposal ratio: Q(M_t+1, M_t) / Q(M_t+1; M_t )
  ##
  ##   where M_t+1 is the next state with an applied cut
  ##         M_t   is current state with not cut (equal to merged)
  ## Q is in both cases uniform discret distributed
  ## Q(M_t+1; M_t) : chance of choosing a segment for a cut ~ all leaf nodes
  ## Q(M_t, M_t+1) : chance of choosing a leaf pair for a merge
 
  ## the number of segments that can be cut corresponds to all the segments in our segment.map
  ## the leaf pairs that can be merged are in the proposed set (M_t+1) , in the case of a real merge move, this is inverted anyways later
  log.proposal.ratio = log( getNrSegments(old.set) / getNrLeafPairs(proposed.set$mondrian.tree) )  
 # cat("FIXME: check log.proposal.ratio\n")  

  ## take out segment prior
  alpha =  log.proposal.ratio + log.prior.ratio + ( (HYPERvar$v0 / 2)*log(HYPERvar$gamma0 / 2) - lgamma(HYPERvar$v0 / 2) - log((HYPERvar$delta2 + 1)^((nr.edges + 1) / 2)) ) + sumPhiPlus - sumPhi 
 
  
  return(alpha)

}


computeProjection = function(x, delta2){
  # INPUT: len, delimiting breakpoints.
  #        x, the observations of X in the corresponding state
  #        delta2.
  # OUTPUT: the projection matrix Px.
  # depends on: .

  ## the number of rows == number of locations == number of elements in the target vector == dim of projection matrix
  len = dim(x)[1]
  
  moins = matrix(0,len,len)

  if(prod(dim(x))>0){
    moins = (delta2/(delta2+1))* x%*%ginv(t(x)%*%x)%*%t(x)
  }

  Px=diag(1,len)-moins


  if(dim(Px)[1] != len) {
    cat("inside computeProjection\n")
    #browser()
  }
  if(len==1) {
    browser()
  }
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



sampleValidateCPs <-function (candidateCPs, min.seglocs, E, E.other, Y, ALTERX, xlocs, type, cp.pos) {

  ##
  ## This part is specific to 2D change-points and for the case that not all locations in the grid a valid (sampled) patches
  ##   With small segments it can happen that only very little actual sampled locations exist, if there is only 1 location for instance
  ##   the method will crash after calculating the projection matrix.
  ##   The code here makes sure each segment will have a minimum number of valid locations to work with, if not the change-point is rejected

  ## this is returned
  cp.new = NaN

  ## loop as long no CP was valid and elements in candidateCPs exist
  while(length(candidateCPs) > 0) {

    ## flag for invalid cp found
    invalid.cp = FALSE

    ## sample the new changepoint uniformly
    cp.tmp = sample(c(candidateCPs, candidateCPs),1)

    ## remove in the case cp.tmp was invalid and while continuous
    candidateCPs = candidateCPs[-which(candidateCPs == cp.tmp)]

    ## create candidate changepoint vector E
    if(type == "shift") {     ## replace the cp when shifted 
      E.tmp = E
      E.tmp[cp.pos] = cp.tmp
    } else {                  ## insert the cp 
      
      E.tmp = sort(c(E,cp.tmp))
    }

    ## get position of new cp
    E.new.pos = which(E.tmp == cp.tmp)
    
    ## extract cp of affected segments (left and right of new one )
    E.tmp = E.tmp[c(E.new.pos - 1, E.new.pos, E.new.pos + 1)]


    for(a.coord in 1:(length(E.tmp)-1)) {
      start.c1 = E.tmp[a.coord]
      end.c1 = E.tmp[a.coord+1]
      
      for (b.coord in 1:(length(E.other)-1) ) {
        
        start.c2 = E.other[b.coord]
        end.c2 = E.other[b.coord+1]

        ## construct actual coordinates, note the 3rd and 4th coord. are -1 because the segment only extends but not includes the right/bottom changepoint
        if(ALTERX) { segcoord = c(start.c1, start.c2, end.c1-1, end.c2-1) }
        else { segcoord = c( start.c2, start.c1, end.c2-1, end.c1-1) }

        ## extract valid locations
        y = extractNodes(Y, segcoord, xlocs, F)

        ## check if size is sufficient, then mark as invalid
        if(length(y) < min.seglocs) {
          #cat("nonvalid: ", length(y), " - ", type, ", ALTERX: ", ALTERX, "\n") 
          invalid.cp = TRUE
          break
        }
      }
    }

    ## check if cp creates valid segments
    if(!invalid.cp) {

      cp.new = cp.tmp
      break
    }
      
  }

  return(cp.new)
}

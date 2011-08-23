

integer.base.b <- function(x, b=2, ndigits=NULL) {
  
        xi <- as.integer(x) 

        if(any(is.na(xi) | ((x-xi)!=0))) 
                print(list(ERROR="x not integer", x=x)) 

        # check for 0, couldn't handle this before
        if(xi == 0) {
           if(is.null(ndigits)) {
             return(array(0, dim=1)) 
           } 
           return(array(0, dim=ndigits))  
        }
        
        N <- length(x) 

        xMax <- max(x)	

        # if not predefined, set number of needed binary digits
        if(is.null(ndigits)) {
          ndigits <- (floor(logb(xMax, base=2))+1) 
        }

        # create the binary vector with size ndigits
        Base.b <- array(NA, dim=c(N, ndigits)) 

        for(i in 1:ndigits){#i <- 1 
                #Base.b[, ndigits-i+1] <- (x %% b)
                Base.b[,i] <- (x%%b)
                x <- (x %/% b) 
        } 
        if(N ==1) Base.b[1, ] else Base.b 
} 

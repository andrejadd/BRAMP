
##
## Inverse Gamma
##

## evaluate the inverse gamma density
##
## Kevin Rompala 5/6/2003
## fixed KQ 3/8/2005

"dinvgamma" <-
  function(x, shape, scale = 1) {

    # error checking
    if(shape <= 0 | scale <=0) {
      stop("Shape or scale parameter negative in dinvgamma().\n")
    }
    
    alpha <- shape
    beta <- scale
   
    # done on log scale to allow for large alphas and betas
    log.density <- alpha * log(beta) - lgamma(alpha) -
       (alpha + 1) * log(x) - (beta/x)
    return(exp(log.density))
  }

## generate draws from the inverse gamma density (using
## the gamma simulator)
##
## Kevin Rompala 5/6/2003
## fixed KQ 3/8/2005

rinvgamma <- function(n, shape, scale) {

# AA , Debugging  
#  tryCatch({
#    aa = 1 / rgamma(n, shape=shape, scale=1/scale)
#  }, warning = function(w) {
#    print(w)
#    browser()
#  })
  
  return(1 / rgamma(n, shape=shape, scale=1/scale))
}

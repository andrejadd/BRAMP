

print.bramp <- function(x, ...) {
  
  cat("Call:\n")
  print(x$call)
  
  
  iters = x$MCMC.chain$Structsamples$iter
  cat("Chain iterations: ", iters[[length(iters)]], "\n")
  
      
  cat("\nEdge probabilities:\n")
  print(x$edge_probs)
  
  cat("\nMean coefficients:\n")
  print(x$mean_coefficients)
}

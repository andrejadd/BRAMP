source("runBRAM.R")

network = commandArgs()[4]
target  = commandArgs()[5]
run = commandArgs()[6]
bestPredictors = commandArgs()[7]
maxiter = commandArgs()[8]
start.budget = commandArgs()[9]



runBRAM(dataid=as.integer(network), target=as.integer(target), runid=as.integer(run), niter=as.integer(maxiter), start.budget=as.double(start.budget))


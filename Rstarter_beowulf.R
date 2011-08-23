source("TVDBNexample.R")

network = commandArgs()[4]
target  = commandArgs()[5]
run = commandArgs()[6]
bestPredictors = commandArgs()[7]
maxiter = commandArgs()[8]


TVDBNexample(modelid=as.integer(network), target=as.integer(target), runid=as.integer(run), niter=as.integer(maxiter))


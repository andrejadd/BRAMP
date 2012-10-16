source("runMethod.R")

network = commandArgs()[4]
target  = commandArgs()[5]
run = commandArgs()[6]
maxiter = commandArgs()[7]
data.prefix = commandArgs()[8]

runMethod(dataid=as.integer(network), target=as.integer(target), runid=as.integer(run), niter=as.integer(maxiter), data.prefix=data.prefix)


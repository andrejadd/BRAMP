source("TVDBNexample.R")

modelid = commandArgs()[4]
target  = commandArgs()[5]
run = commandArgs()[6]

TVDBNexample(modelid=as.integer(modelid), target=as.integer(target), runid=as.integer(run))



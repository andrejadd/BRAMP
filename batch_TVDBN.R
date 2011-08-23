

# networkid=commandArgs()[4]

source("TVDBNexample.R")

for(networkid in 1:30) {

  
  print("start with network: ", networkid)

  for(trials in 1:30) {
    TVDBNexample(networkid=networkid)
  }

}

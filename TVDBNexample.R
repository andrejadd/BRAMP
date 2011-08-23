

TVDBNexample <- function(modelid=NULL, target=NULL, runid=NULL, niter=NULL){

# remove all but arguments  
rm(list= ls()[ls()!="modelid" && ls()!= "target" && ls()!="runid" && ls()!="niter"])

# modeltype means static (0 would be dynamic, only used to read in appropriate ymatrix .Rdata files)
modeltype = 1

cat("dataid   : ", modelid, "\n")
cat("type : ", modeltype, "\n")
cat("target    : ", target, "\n")
cat("runid     : ", runid, "\n")


### specify here the TVDBN directory 
path="./"

# Unix time (seconds since 1970), should be fine for the cluster, might add/subtract milliseconds 
set.seed(as.numeric(Sys.time()))

## Import code
# code path
codePath=paste(path,"Code/",sep="")
source(paste(codePath,"runtvDBN.R",sep=""))

## function choosePriors:
## This function plots some examples of prior distribution for the number of breakpoints (changepoints, CP) and the number of incoming edges. 
## The prior distribution for the number of CP is a truncated poisson, with maximum=maxk and mean=lambda which sampled from a Gamma distribution:
##  k ~ Gamma(alphaCP,betaCP)  (default: alphaCP=0.5, betaCP =1)
## same thing for the number of incoming edge with maximum= maxTF (for Transcription Factor) and hyperparameters (alphaTF, betaTF)
## (Remark: truncated Poisson could be changed into a uniform ditribution)
#choosePriors(5,paste(codePath,"k_priors.txt",sep=""))

# read in the data
indata = paste(path,"../Data/Data_id", modelid, ".Rdata",sep="")
cat("read Rdata file : ", indata, "\n")

## the main sample data that is loaded is a target/predictor matrix [nodes x locations] and a SAC (spatial autocorrelation) matrix with
## dame size. The SAC matrix has for each target and location a SAC node which is included into the linear regression 
load(indata)

# total sample points
n=Model$xlocs*Model$ylocs

cat("total location points: " , n, " with x: ",  Model$xlocs, " , y: ",  Model$ylocs, "\n")

# number of parent nodes
q=dim(Model$Ymatrix)[1]-1

cat("nr. parent nodes: ", q, "\n")

# Maximum number of parent nodes (fan-in restriction). A low limit is needed 
# for birth proposals based on precomupted posterior distribution (method 4),
# otherwise you can set smax = q
smax = min(8,q);

# maximal number of CPs, tweak this, One option: make dependent of number of locations -> but then each axis should have its own kmax (FIXME)
kmax = max(floor((( Model$xlocs+ Model$ylocs)/2) / 10), 5)

# minima length of a segment (or a phase)
minPhase=2

#hyperparameters for the number of CP and incoming edges (TF) 
alphaCP=1
betaCP=0.5
alphaTF=1
betaTF=0.5

#1# run TVDBN procedure:
runtvDBN(fullData=Model$Ymatrix,
         sacData=Model$SAC.nodes,
         n=n,
         xlocs=Model$xlocs,
         ylocs= Model$ylocs,
         q=q,
         minPhase=minPhase, 
         kmax=kmax,
         smax=smax,
         alphaCP=alphaCP, betaCP=betaCP, alphaTF=alphaTF, betaTF=betaTF,
         niter=niter,  
	 outputFile=outputFile, 
	 modelid=modelid, target=target, runid=runid )

}

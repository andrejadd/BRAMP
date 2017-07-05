

save.file = TRUE;
gene.number = 10;

## AA
timeseries.length = 500;

nrsamples = 100;

cp_dist = list()

networkid = 6

filename = paste("CPsamples_n", networkid, "_i", sep="")
  

# CP counter
cp_dist = matrix(0, 1, timeseries.length+1);
  
for(gene in 1:gene.number) {
    
    # Load file containing changepoints (results$cp_samples) and networks
    # (results$edge_samples)
    load(file=paste(filename, gene, sep=""))
    
    cat("Gene: ", (gene - 1), " , nr. of CP samples: ", length(CPsamples), "\n")
    	
    # Calculate change-point posterior distribution
    cp_dist_tmp = matrix(0, 1, timeseries.length+1);
    
    # Sum changepoints for gene and increase CP counter for specific time position cp_inds.
    # -- What data is this

    if(length(CPsamples) >1) {

      for(sample in 1:length(CPsamples)) {
    
        cp_inds = CPsamples[[sample]]
        cp_dist_tmp[cp_inds] = cp_dist_tmp[cp_inds] + 1;
      
      }
    }
    
        
    # Add posterior distribution of changepoints for the current gene
    # Basically using P(AuB) = P(A) + P(B) - P(AnB) where I've assumed 
    # independence, i.e. P(AnB) = P(A)*P(B)
    
    cp_dist = cp_dist + (cp_dist_tmp/nrsamples) - ( (cp_dist_tmp/nrsamples) * cp_dist);

#    print.table(cp_dist)
      
  

}

# Plot changepoint posterior distribution

if(save.file) {
  postscript(file = paste("CPsamples_PPavg_n", networkid, ".eps", sep=""))
  par(lwd=4, ps=20)
}

plot(1, 1, type='n', ylim=c(0,1), xlim=c(1, timeseries.length),xlab="Timepoints", ylab="Posterior Probability");

#legend(1, 1, c("HetDBN-0", "HetDBN-Exp", "HetDBN-Bino1", "HetDBN-Bino2"), 
#    lty=c(1,1,1,1), pch=c(-1,1,4,0), y.intersp=1.5, lwd=2)


points(1:timeseries.length, cp_dist[,1:timeseries.length], type='l', xlab="Timepoints", ylab="Posterior Probability", ylim=c(0,1), lty=1);


#abline(v=30, lty=3)
#abline(v=75, lty=3)


#abline(v=130, lty=3)
#abline(v=390, lty=3)

#abline(v=13, lty=3)
#abline(v=26, lty=3)
#abline(v=39, lty=3)
#abline(v=52, lty=3)

#abline(v=12, lty=3)
#abline(v=35, lty=3)
#abline(v=47, lty=3)
#abline(v=61, lty=3)
#abline(v=75, lty=3)
#abline(v=90, lty=3)


if(save.file) {
dev.off()
}

#phases = list()
#browser()

#results = list(networks=networks, thresholded=thresholded);

#save('results', file=
# paste("../Results/drosophila_networks_exp.RData", sep=""))


#for(gene in 1:11) {
#  load(file=paste(filename, gene, sep="_"))
#  browser()
#  print(paste("Gene", gene))
#  print("Calculating Marginal Networks")
#  marginal_networks = marginal_comp(results$edge_samples, 9, 11, 
#    results$cp_samples)
#  
#  for(phase in 1:4) {
#    networks[[phase]][,gene] = marginal_networks$margin_net[[4]][,phase];
#  }
#  
#}
#
#for(phase in 1:4) {
#  networks[[phase]] = matrix(0, 9, 9)
#}
#
#print(networks)

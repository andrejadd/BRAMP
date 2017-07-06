\name{BRAMP}
\alias{BRAMP}
\title{Hierarchical Bayesian Regression and Mondrian Process for inference of species networks and spatial partitions.}
\description{
  Runs MCMC simulations for network inference in a two-dimensional grid.
}
\usage{
BRAMP(y, X, y_SAC_node, xlocs, ylocs, nr_iterations, chain_thinout, mcmc_obj, mcmc_rnd_seed, fixed_edges, edge_fanin)

}
\arguments{
  \item{y}{A m-length vector with the number of 'm' observations for the target node.}
  \item{X}{A n-by-m matrix with 'n' nodes and 'm' observations.}
  \item{y_SAC_node}{A m-length vector with the spatial autocorrelation data. Can be left NULL.}
  \item{xlocs}{An integer defining the number of observations along the x-axis.}
  \item{ylocs}{An integer defining the number of observations along the y-axis.}
  \item{nr_iterations}{The number of MCMC iterations to run.}
  \item{chain_thinout}{Save samples from the chain every 'chain_thinout' iteration.}
  \item{mcmc_obj}{A MCMC data structure of class 'bramp'. The same that
    is returned from BRAMP(). If this is passed as argument back to
    BRAMP _and_ 'nr_iterations' is greater than the iterations in 'mcmc_obj', the MCMC simulation is continued.}
  \item{mcmc_rnd_seed}{A integer value specifying the seed for the random generator.}
  \item{fixed_edges}{A n-length vector that defines fixed and active edges.}
  \item{edge_fanin}{The maximum number of incoming parents.}
}
\seealso{
  \code{\link{save}}.
}
\examples{
## save all data
save(list = ls(), file= "all.RData")

## restore the saved values to the current environment
load("all.RData")

## restore the saved values to the workspace
load("all.RData", .GlobalEnv)
}
\keyword{baysian regression}

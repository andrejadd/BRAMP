# BRAMP

This is the software implementation of the Bayesian regression and Mondrian process (BRAMP) model published here:

`Aderhold, Andrej, Dirk Husmeier, and V. Anne Smith. "Reconstructing ecological networks with hierarchical Bayesian regression and Mondrian processes." AISTATS. 2013.`

The main function is BRAMP(), which takes an Rdata input file, the number of iterations to run, and the target node for which to run the simulation as input. See run_Example.R on how to call BRAMP().


# Install

Inside R, execute the following commands:

```r
install.packages("devtools")
devtools::install_github("andrejadd/BRAMP")
``` 

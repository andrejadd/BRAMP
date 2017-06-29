# BRAMP

This is the software implementation of the Bayesian regression and Mondrian process (BRAMP) model published here:

`Aderhold, Andrej, Dirk Husmeier, and V. Anne Smith. "Reconstructing ecological networks with hierarchical Bayesian regression and Mondrian processes." AISTATS. 2013.`

The main function is BRAMP(), which takes an Rdata input file, the number of iterations to run, and the target node for which to run the simulation as input. See run_Example.R on how to call BRAMP().



## TODO

- fix move proposals in Code/main.R
- put delta2 where it belongs (its under Hyperparameters but should it?)
- bind edge weights to mondrian tree node instead of using a seperate lookup matrix - is it really feasable or makes things more complex?
- clean up more and simplify, are all fct. gone that I don't need?
- add the synthetic data generator? or better to leave it in one package? 
- put in basic evaluation script that is currently in another place.


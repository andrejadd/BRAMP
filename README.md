# BRAMP

This is the software implementation of the Bayesian regression and Mondrian process (BRAMP) model published here:

`Aderhold, Andrej, Dirk Husmeier, and V. Anne Smith. "Reconstructing ecological networks with hierarchical Bayesian regression and Mondrian processes." AISTATS. 2013.`

Needs cleaning up and proper generic entry function (replacing current run_Method.R) to make it easier to use for everyone else. 


## TODO

- fix move proposals in Code/main.R
- put delta2 where it belongs (its under Hyperparameters but should it?)
- bind edge weights to mondrian tree node instead of using a seperate lookup matrix - is it really feasable or makes things more complex?
- clean up more and simplify, are all fct. gone that I don't need?
- add the synthetic data generator? or better to leave it in one package? 
- put in basic evaluation script that is currently in another place.

Maybe not important for future use because its specific to an SGE/PBS cluster.:

- make the push.dsh.jobs.py script easier to read, maybe use a parameter vector such as job.params = [waittime, dataids, targets, run.ids, iterations, additional.param]. From this the JOB.queue is then created and I can have all important in one line

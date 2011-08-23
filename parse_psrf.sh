#!/bin/bash

# grep -l Network2 s_tvdbn*` 


#for a in `grep -li "Response variable  1" \`grep -li "Network1" s_tvdbn_n1.*\``
for a in `grep -li "Network2" s_tvdbn_n1.*`
do 
  #echo $a
	#if [  $a == "TVDBNexample.R" ]; then cat $a; fi; done
  gawk '{ if(/currpsrf/) {print $2,$4;} }' $a
done 



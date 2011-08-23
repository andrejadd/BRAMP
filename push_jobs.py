#!/usr/bin/python

import os
network = 2

# these nodes get started
#nodelist = [1,2,3]
nodes = 10

repeats = 5

for i in range(repeats):
    # loop through nodes and start
    for node in range(nodes):
# for node in nodelist:

        # create run file
        f = open('./s_dbn_v2', 'w')

        f.write('cd ~/TVDBN_static\n')
        f.write('cat Rtest_singleNode.R | R --vanilla --args ' + str(network) + " " +  str(node-1) + " " + str(i) + "\n")
        f.close()

        # exec qsub
        os.system('qsub -q idle.q ./sh_dbn_v2')



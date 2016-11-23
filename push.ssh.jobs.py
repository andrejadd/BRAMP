#!/usr/bin/python
#
#
#
# change log:
#    see 'git log'

# to change: maxiter, nodes.config


import os
import os.path
import time
import subprocess
import sys
import glob
import pdb   # to use debugging , pdb.set_trace()

######## PARAMETER ------------------------

SSH_PUBLICKEY_CONNECT = True   # set this to false if only to run it locally without dsh,ssh (when publickey ssh fails)
CPU_USAGE = 0.5

BEOWULFDIR = "../Helper/Beowulf_root/"  # here are the beowulf config settings 
WORKINGDIR = os.getcwd()         # get current dir, need to change to this on each node
STARTSCRIPTS = "Startscripts"
RESULTDIR = "Results/"  # use this to check if runs already exist
NODELOGDIR = "Nodelogs"

nodefile = BEOWULFDIR + '/nodes.list.euclids'
idfile = BEOWULFDIR + '/jobid.last'


interDispatchTime = 30           # time in seconds to wait until a full dispatch attempt is made again 
start_run_id = 1
nr_runs = 1                      # number of runs per unique job 
defaultnodes = range(1,11)      # the nodes to compute
maxiter = 400000                 # nr. of MCMC iteration steps


## Data ------------------------------------------

datasets = []

## Lotka-Volterra Simulated data
#datasets.append([range(201,231), "LOTKA.VOLTERRA"])
#datasets.append([range(801,831), "LOTKA.VOLTERRA", "LOTKA.VOLTERRA"])
datasets.append([range(401,431), "LOTKA.VOLTERRA", "LOTKA.VOLTERRA"])
#datasets.append([range(601,631), "LOTKA.VOLTERRA"])


## Information Sharing data with different epsilons
#datasets.append([range(1,21), "SYNTHETIC.INFSHARING.EPSILON.0", "SYNTHETIC.INFSHARING/MONDRIAN_CP/RndNodes.WithoutErrorBias"])
#datasets.append([range(1,21), "SYNTHETIC.INFSHARING.EPSILON.0.125", "SYNTHETIC.INFSHARING/MONDRIAN_CP/RndNodes.WithoutErrorBias"])
#datasets.append([range(1,21), "SYNTHETIC.INFSHARING.EPSILON.0.25", "SYNTHETIC.INFSHARING/MONDRIAN_CP/RndNodes.WithoutErrorBias"])
#datasets.append([range(1,21), "SYNTHETIC.INFSHARING.EPSILON.0.5", "SYNTHETIC.INFSHARING/MONDRIAN_CP/RndNodes.WithoutErrorBias"])
#datasets.append([range(1,21), "SYNTHETIC.INFSHARING.EPSILON.1", "SYNTHETIC.INFSHARING/MONDRIAN_CP/RndNodes.WithoutErrorBias"])



## add  exceptions here 
nodehash = {}

# read next available jobid from BEOWULF directory
try:
	f = open(idfile, 'r')
	jobstr = f.readline()
	jobid = int(jobstr.rstrip("\n"))
	f.close()

except IOError as (errno, strerror):
    print "I/O error({0}): {1}".format(errno, strerror) + " : " + idfile
    sys.exit()
except ValueError:
    print "Could not convert data to an integer."
    sys.exit()

print "Starting with job-id " + str(jobid)


# this list will include one entry per job and is constructed below
JOBs = []

# set vector of runs
runids = range(start_run_id,(start_run_id + nr_runs)) 
#pdb.set_trace()

for run in runids:
	for dset in datasets:
		for network in dset[0]:

			if nodehash.has_key(network):
				nodes = nodehash[network]
			else:
				nodes = defaultnodes
				
			for node in nodes:

				JOBs.append([jobid, network, node, run, maxiter, dset[1], dset[2]])
				print "added job: " + str(jobid) + ", model: " + str(network) + ", node: " + str(node) + ",run: " + str(run) + ", dataprefix: " + dset[1]
				jobid += 1
				





# write back jobid + 1
try:
	idfile = BEOWULFDIR + '/jobid.last'
	f = open(idfile, 'w')
	f.write(str(jobid + 1))
	f.close()

except IOError as (errno, strerror):
    print "I/O error({0}): {1}".format(errno, strerror) + " : " + idfile
    sys.exit()




# read in the available nodes from the config file
cnodelist = []

try:
	f = open(nodefile, 'r')

	for line in f:
		line = line.rstrip("\n")
		
		# exclude lines starting with # 
		if line.find("#",0,1) == -1:
			cnodelist.append(line)	
		
	f.close()

except IOError as (errno, strerror):
    print "I/O error({0}): {1}".format(errno, strerror) + " : " + nodefile
    sys.exit()
except ValueError:
    print "Could not extract nodes from config."
    sys.exit()

print cnodelist




# matrix, each row for one node, column 1 : node name, 2 : nr. CPUs, 3 : load avg 
# is filled on the CPU detect loop below
CPUs = []

if SSH_PUBLICKEY_CONNECT:
	# read out the nr. of cores per node
	for cnode in cnodelist:
		#tmppipe = subprocess.Popen('dsh -m ' + cnode + ' -- grep processor /proc/cpuinfo ', shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	
		tmppipe = subprocess.Popen('ssh ' + cnode + ' grep processor /proc/cpuinfo ', shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		# get output from pipe/
		procinfos = tmppipe.communicate()

		# check for error
		# check if there were errors (probably ssh failure)
		if len(procinfos[1]) > 0:   # some text here
			print("failure on cpuinfo " + cnode + ": " + procinfos[1] + "\n")
			continue
                                
		# get nr. of cpus on the node
		nrprocs = procinfos[0].count("processor") 
		print("cnode: " + cnode + " , nr. CPUs: " + str(nrprocs))
	
		# save ,this matrix is used to assign jobs
		CPUs.append([cnode, nrprocs, -1])
else:
	tmppipe = subprocess.Popen('grep processor /proc/cpuinfo ', shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	
	# get output from pipe
	procinfos = tmppipe.communicate()

	# check if there were errors (probably ssh failure)
	if len(procinfos[1]) > 0:   # some text here
		print("failure to get cpuinfo from /proc/cpuinfo\n")
	                                
	# get nr. of cpus on the node
	nrprocs = procinfos[0].count("processor") 
	print("  local nr. cores: " + str(nrprocs))
	
	# save ,this matrix is used to assign jobs
	CPUs.append(['localhost', nrprocs, -1])
	


# loop as long there are jobs present
while len(JOBs) > 0: 

	# get the load averages for the last minute for each node 	
	for nodeentry in CPUs:

		cnode = nodeentry[0]

		# get uptime info
		if SSH_PUBLICKEY_CONNECT:
			#tmppipe = subprocess.Popen("dsh -m " + cnode + " -- uptime", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

			tmppipe = subprocess.Popen("ssh " + cnode + " uptime", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		else:
			tmppipe = subprocess.Popen("uptime", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			
        	# returns a list, first element -> stdout, second -> stderr
		uptime = tmppipe.communicate()

        	# check if there were errors (probably ssh failure)
        	if len(uptime[1]) > 0:   # some text here
			print("failure getting uptime: \n\t  " + uptime[1] + "\t  .. set lavg to -1")
			nodeentry[2] = -1
			continue	

		# extract load avg 
		start = uptime[0].index("age:") + 4
	
		# here it is 
		minloadavg = float(uptime[0][start:start+4])

		print("load avg 1min: " + str(minloadavg))

		nodeentry[2] = minloadavg




	# check each processor and put a job on free slots (cpus)
	for nodeentry in CPUs:

		cnode = nodeentry[0]

		# in case getting lavg failed before, skip (will be tried next time)
		if nodeentry[2] == -1:
			continue
	
		# since load avg is not normalized in respect to the nr. of cpus it can be just subtracted
		cores_available = int(round(nodeentry[1] * CPU_USAGE))
		slotnr = max(0, int(round(cores_available - nodeentry[2])))

		print("cnode: " + cnode + ", cpus: " + str(nodeentry[1]) + ", lavg: " + str(nodeentry[2]) + ", slots: " + str(slotnr) + " (of " + str(cores_available) + " cores available)" )
		
		for slot in range(0,slotnr):

			# check if still jobs there
			if len(JOBs) == 0:
				print("no more jobs, finishing")
				sys.exit()
				
			# get and remove first job
			nextjob = JOBs.pop(0)
			
			print("\tpush job: " + str(nextjob))
			

			# define log file
			nodelogstd = NODELOGDIR + "/log_" + nextjob[5] + "_job" + str(nextjob[0]) + "_id" + str(nextjob[1]) + "_n" + str(nextjob[2]) + "_" + cnode
			nodelogerr = NODELOGDIR + "/err_" + nextjob[5] + "_job" + str(nextjob[0]) + "_id" + str(nextjob[1]) + "_n" + str(nextjob[2]) + "_" + cnode
			

			# create run file
			runfile = STARTSCRIPTS + "/job_ID" + str(nextjob[0]) + "_D" + str(nextjob[1]) + "_N" + str(nextjob[2]) + "_r" + str(nextjob[3]) + ".sh"

                        f = open(runfile, 'w')

                        # run the R code (niced)
                        f.write("nice -n 19 cat Rstarter_beowulf.R | nice -n 19 R --vanilla --args " + str(nextjob[1]) + " " +  str(nextjob[2]) + " " + str(nextjob[3]) + " " + str(nextjob[4]) + " " + str(nextjob[5]) + " " + str(nextjob[6]) + " >> " + nodelogstd + " 2>> " + nodelogerr + "\n")
                
                        f.close()

			# set executeable
			os.system("chmod u+x " + runfile)

                        # exec the shell script on the node
			if SSH_PUBLICKEY_CONNECT:
				cmttmp = 'ssh ' + cnode + ' \" cd ' + WORKINGDIR + ' ; ' + runfile + ' \"  &'
				
				#cmttmp = 'dsh -m ' + cnode + ' -- \" cd ' + WORKINGDIR + ' ; ' + runfile + ' \"  &'
			else:
			#	cmttmp = 'cd ' + WORKINGDIR + ' ; ' + runfile + ' &'
				cmttmp = runfile + ' &'
		
			##pdb.set_trace()
 			subprocess.Popen(cmttmp, shell=True)

			# no need to sleep but doing it anyways
			time.sleep(1)



	print("\nsleeping " + str(interDispatchTime) + " seconds\n")
	time.sleep(interDispatchTime)			
		

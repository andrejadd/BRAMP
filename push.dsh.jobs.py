#!/usr/bin/python
#
#
#
# change log:
#    15.02.2011   for each data, a range of nodes can be edited

# to change: maxiter, nodes.config


import os
import os.path
import time
import subprocess
import sys
import glob

######## PARAMETER ------------------------


BEOWULFDIR = "../Helper/Beowulf_root/"  # here are the beowulf config settings 
WORKINGDIR = os.getcwd()         # get current dir, need to change to this on each node
STARTSCRIPTS = "Startscripts"
RESULTDIR = "Results/"  # use this to check if runs already exist
NODELOGDIR = "Nodelogs"

nodefile = BEOWULFDIR + '/nodes.config'
idfile = BEOWULFDIR + '/jobid.last'

startbudget = 1

interDispatchTime = 300           # time in seconds to wait until a full dispatch attempt is made again 
start_run_id = 2
nr_runs = 1                       # number of runs per unique job 

defaultnodes = range(1,11)
maxiter = 5000                 # nr. of MCMC iteration steps

bestPredictors = 20 ## is obselete, keep for larger networks or remove


############# THE DATA ------------------------------------------

networks = range(40,43)

# weak predation:
#networks = range(201, 231) + range(401,431) + range(601,631) + range(801,831)
#networks = range(100801,100820)

#networks = range(161,171)

#networks = range(1000101,1000121) +  range(1000601,1000621) + range(1000801,1000821) 

#networks = range(601,630) + range(801,830)

#networks = range(100001,100031) + range(100101,100131) + range(100201,100231) + range(100301,100331) + range(100401,100431) + range(100601,100631) + range(100801,100831)

#networks = range(100101,100131) + range(100201,100231) + range(100301,100331) + range(100401,100431) + range(100601,100631) + range(100801,100831)
#networks = range(100201,100221) + range(100401,100421) + range(100601,100621) + range(100801,100821)



## Use this, if the different data sets have varying numbers of nodes
nodehash = {}

## e.g.:
#nodehash[3011] = range(1,15)




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

for network in networks:
	for run in runids:

		if nodehash.has_key(network):
			nodes = nodehash[network]
		else:
			nodes = defaultnodes

		for node in nodes:

			# check if log file already exists
			lfile = NODELOGDIR + "/log_oID*_D" + str(network) + "_N" + str(node) + "_*"
									
			# check if run already exists
			dfile = RESULTDIR + "SC2D_m" + str(network) + "_i" + str(node) + "_run" + str(run)
			
			logexists = "no"
			resexists = "no"

			if len(glob.glob(lfile)) == 1:
				logexists = "yes"
			if os.path.exists(dfile):
				resexists = "yes"
				

			# only create job when log file not exists, a result file not be overwritten by BRAM
#			if resexists == "no" and logexists == "no":
#			if logexists == "no":

			## FIXME: Need to rethink this> check for existing log file is not really working since it does not considers different runs
			##        So everything relayed on the existence of a result file. This checks now BRAM in order to be able to proceed with 
			##        a MCMC chain (open, check chain and iterations, continue or exit)

			JOBs.append([jobid, network, node, run, bestPredictors, maxiter, startbudget])
			print "adding job: " + str(jobid) + ", model: " + str(network) + ", node: " + str(node) + ",run: " + str(run)
				# increment, although, job not created, just to not conflict with running jobs ids and logs
			jobid += 1

#			else:

#				print "skip job, log exists: " + logexists + ", result exists: " + resexists



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

# read out the nr. of cores per node
for cnode in cnodelist:
	tmppipe = subprocess.Popen('dsh -m ' + cnode + ' -- grep processor /proc/cpuinfo ', shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	
	# get output from pipe
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



# loop as long there are jobs present
while len(JOBs) > 0: 

	# get the load averages for the last minute for each node 	
	for nodeentry in CPUs:

		cnode = nodeentry[0]

		# get uptime info
		tmppipe = subprocess.Popen("dsh -m " + cnode + " -- uptime", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
           
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
		slotnr = max(0, int(round(nodeentry[1] - nodeentry[2])))
		
		print("cnode: " + cnode + ", cpus: " + str(nodeentry[1]) + ", lavg: " + str(nodeentry[2]) + ", slots: " + str(slotnr))
		
		for slot in range(0,slotnr):

			# check if still jobs there
			if len(JOBs) == 0:
				print("no more jobs, finishing")
				sys.exit()
				
			# get and remove first job
			nextjob = JOBs.pop(0)
			
			print("\tpush job: " + str(nextjob))
			

			# define log file
			nodelogstd = NODELOGDIR + "/log_oID" + str(nextjob[0]) + "_D" + str(nextjob[1]) + "_N" + str(nextjob[2]) + "_" + cnode
			nodelogerr = NODELOGDIR + "/log_eID" + str(nextjob[0]) + "_D" + str(nextjob[1]) + "_N" + str(nextjob[2]) + "_" + cnode
			

			# create run file
			runfile = STARTSCRIPTS + "/job_ID" + str(nextjob[0]) + "_D" + str(nextjob[1]) + "_N" + str(nextjob[2]) + "_r" + str(nextjob[3]) + ".sh"

                        f = open(runfile, 'w')

                        # run the R code (niced)
                        f.write('nice -n 19 cat Rstarter_dshjob.R | nice -n 19 R --vanilla --args ' + str(nextjob[1]) + " " +  str(nextjob[2]) + " " + str(nextjob[3]) + " " + str(nextjob[4]) + " " + str(nextjob[5])+ " " + str(nextjob[6])+  " >> " + nodelogstd + " 2>> " + nodelogerr + "\n")
                
                        f.close()

			# set executeable
			os.system("chmod u+x " + runfile)

                        # exec the shell script on the node
                        cmttmp = 'dsh -m ' + cnode + ' -- \" cd ' + WORKINGDIR + ' ; ' + runfile + ' \"  &'
 			subprocess.Popen(cmttmp, shell=True)

			# no need to sleep but doing it anyways
			time.sleep(1)



	print("\nsleeping " + str(interDispatchTime) + " seconds\n")
	time.sleep(interDispatchTime)			
		

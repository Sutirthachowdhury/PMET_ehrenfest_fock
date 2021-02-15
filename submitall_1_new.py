#!/usr/bin/python
import os
import shutil
import time
lastNode = 0

offset= 00
nproc = 100


def sbatch(filename):
    """
    Submit a job and get the job id returned
    """
    submit = os.popen("sbatch  %s"%(filename)).read()
    subId = submit.split()[3].replace("\n","")
    return subId


def checkJob(jobId,user="schowdh4"):
    """
    If Job is running returns True
    """
    import os
    SQ = os.popen("squeue -u %s"%(user)).read().split("\n")[1:]
    jobfound = False
    for i in SQ:
	
       	try : 
		idn = int(i.split()[0])
       	except : 
		idn = False
	if (idn == int(jobId)):
		jobfound = True
		if i.split()[4] == "R":
			
			return i.split()[5] 
			
		else : 
		   print "Job status : \033[93m", i.split()[4] , "\033[0m"            
    if jobfound == False  : 
	raise
            
    return False

def getNodeId(jobId, user="schowdh4"):
    import os
    SQ = os.popen("squeue -u %s"%(user)).read().split("\n")[1:]
    for i in SQ:
        try:
            idn = int(i.split()[0])
	    if (idn == int(jobId)) and i.split()[4] == "R": 
		return i.split()[-1] 
	    if (idn == int(jobId)) and i.split()[4] != "R":
		return i.split()[4]  
	except:
	    return False
    return False

def wait(job):
	jobDone = False
	while jobDone == False:
		time.sleep(1)	
		jobDone = checkJob(job)
		time.sleep(3)	
	return jobDone


def nodeList():
	sinfo = os.popen("sinfo -p exciton")[1:]
	for i in sinfo:
		if i.replace("\n","").split()[4] == 'idle' :
			print "Idle List : ", i.replace("\n","").split()[5]
			lists = i.replace("\n","").split()[5].split(",")
			return lists

dirs = []

l = []

if (nproc > 9):
	for i in range(0,10):
		name = 'mapp-0' + str(i)
		dirs.append(name)
	for i in range(10,nproc):
		name = 'mapp-'  + str(i)
		dirs.append(name)
else:
	for i in range(0,nproc):
		name = 'mapp-'  + str(i)
		dirs.append(name)
info = open("submits.dat", "w+") 
currentpath = os.getcwd()
fake = open("fake.dat", "w+")
i = offset

while i < nproc:
        j = i
	workingd = currentpath + '/' + dirs[i]
	os.chdir(workingd)
	thisinfo = open("submit.dat", "w+")
	jobID = sbatch("submit.sbatch") 
	print "Job Submitted", jobID
	try : 
		Node = wait(jobID)
	except :
                print "Job Error! will try a rerun" 
		print "Fake Job to occupy node"
		os.chdir(currentpath)
		fakeID = sbatch("fake.sbatch")
		os.chdir(workingd)
		time.sleep(2) 
		fakeNode = getNodeId(fakeID) 
        	fake.write("%s\t%s\t%s\n"%(dirs[i], jobID, NODE))
        NODE = getNodeId(jobID) 
        if NODE == False : 
	        "Job run error at  \033[91m", dirs[i],  "\033[0m" 
		i = i - 1 
	else : 
		print "\033[94m" , dirs[i],"\033[0m submitted at \033[92m" ,NODE, "\033[0m with id = \033[91m", jobID, "\033[0m"
        	info.write("%s\t%s\t%s\n"%(dirs[i], jobID, NODE))
        	thisinfo.write("%s\t%s\t%s\n"%(dirs[i], jobID, NODE))
	i = i + 1
info.close()
fake.close()

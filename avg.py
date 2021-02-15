import numpy as np 
import sys

nproc = 100

dirs = []


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


filename = sys.argv[1]
dat = np.loadtxt(dirs[0]+"/"+filename) * 0.0 

for i in range(nproc):
    dat = dat + np.loadtxt(dirs[i]+"/"+filename) 
np.savetxt("avg-"+filename, dat/nproc)
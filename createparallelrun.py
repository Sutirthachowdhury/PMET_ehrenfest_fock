#!/usr/bin/python
import os
import shutil

nproc = 50

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

#generate folders
for i in range(0,nproc):
	if not(os.path.exists(dirs[i])):
			os.mkdir(dirs[i])



#copy executables
for i in range(0,nproc):
	shutil.copy2('rpmd.exe',dirs[i])
	shutil.copy2('steps.in',dirs[i])
	shutil.copy2('params.in',dirs[i])
	shutil.copy2('submit.sbatch',dirs[i])

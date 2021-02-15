import sys
import numpy as np
from scipy.optimize import curve_fit
global intP
def fit(time,k1,k2):
	return	-(k1 + k2) * intP + k2 * time - k2 * time[0] * np.ones(len(time)) 

# Input --------------------
fob = np.loadtxt(sys.argv[1]) 
skip = 5  #5
start = 0     # in au
end   = 4134 # in au #41000
P   = fob[::skip,1][int(start/skip):int(end/skip)]
t   = fob[::skip,0][int(start/skip):int(end/skip)] /41341.37 
#print len(t)
#---------------------------
# Integrate P --------------
#---------------------------
dt   = t[1] - t[0]
intP = np.zeros(len(P)) 
for i in range(len(P)):
	intP[i] = sum(P[:i])  * dt 
#--------------------------- 
param = [0.0,0.0]
fitP = P - P[0] * np.ones(len(P))  
popt, pcov = curve_fit(fit, t, fitP , p0 = param , bounds= (0,[10,10]))                                     
Opt = fit(t,*popt)   
Al = np.array([t,Opt,fitP])
Al = np.transpose(Al)
np.savetxt("%s/k.txt"%(sys.argv[2]),Al)
print popt[0] , popt[1] 

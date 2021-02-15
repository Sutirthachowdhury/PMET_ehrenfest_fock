# Parameters
import math
def marcus(g): 
 Vda = (0.01/27.2114)
 pi  =  math.pi
 lam =  0.65/27.2114
 beta = 1052.85
 A  = (Vda**2) * (pi*beta/lam) **0.5 
 A  = A*41341.37
 dG = -g 
 alpha = beta/(4*lam) 
 Exp = math.exp(-alpha*(dG+lam)**2) 
 K = A*Exp
 print K, dG 
 return K

mar = open("marxus.txt","w+") 
gc = 0.0
while gc<=0.045*2:
        k1 = marcus(gc)
        k2 = marcus(gc-1.0/27.2114)
        mar.write("%s\t%s\n"%(gc*27.2114,k1+k2))
        gc+=0.001
mar.close()

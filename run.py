import random
import shutil
import os


def sbatch(filename):
    """
    Submit a job and get the job id returned
    """
    submit = os.popen("sbatch  %s"%(filename)).read()
    subId = submit.split()[3].replace("\n","")
    return subId


def makeinp(Fold= "./"):
    fob = open(Fold+"/params.in","w+")
    
    nb = 1
    fob.write(str(nb) + "\t" + "nb" + "\n")

    nstep = 4100
    fob.write(str(nstep) +  "\t" + "nstep" + "\n")

    ntraj = 300
    fob.write(str(ntraj) +  "\t" + "ntraj" + "\n")
    
    beta = 1052.85
    fob.write(str(beta) +  "\t" + "beta" + "\n")

    masse = 1.0
    fob.write(str(masse) +  "\t" + "masse" + "\n")

    xi =  678693.24
    fob.write(str(xi) +  "\t" + "xi" + "\n")

    wc =  3.50690243E-4
    fob.write(str(wc) +  "\t" + "wc" + "\n")

    mnuc = 452.9057E6
    fob.write(str(mnuc) +  "\t" + "mnuc" + "\n")

    omega_c = 0.2/27.2114
    fob.write(str(omega_c) +  "\t" + "omega_c" + "\n")
    
    delta = 0.005/27.2114
    fob.write(str(delta) +  "\t" + "delta" + "\n")

    f0 = 55.7
    fob.write(str(f0) +  "\t" + "f0" + "\n")

    shift = 0.00629177
    fob.write(str(shift) +  "\t" + "shift" + "\n")
    
    bias = 0.15/27.2114
    fob.write(str(bias) +  "\t" + "bias" + "\n")

    dt = 10.0
    fob.write(str(dt) +  "\t" + "dt" + "\n")

    ourseed = random.randint(0,1E6)
    fob.write(str(ourseed) +  "\t" + "ourseed" + "\n")
    	

    fob.close()
    

    fob = open(Fold+"/steps.in","w+")

    fob.write("0.5d0\t,0.1d0\t,0.1d0\n")
    fob.write("step(1)\t,step(2)\t,step(3)\n")
    


   # fob.write("413000\t413000\t20\t200\t %s\t %s\n"%(str(G),str(random.randint(0,10E6)))) 
   # fob.write("runtime\tnsteps\tmappingsteps\ttrajectory\tG(eV)\tRandom\n")
   

def md(folder):
    try:
        os.mkdir(folder) 
    except:
        print "Folder exists:", folder
def cp(filename,folder):
    shutil.copy2(filename,folder)

nfold = 100

for i in range(nfold): 
    main = "mapp-" + str(i)
    md(main) 
    makeinp(main)
    cp('rpmd.exe',main)
    cp('submit.sbatch',main)
    os.chdir(main)
    sbatch("submit.sbatch")  
    os.chdir("../")

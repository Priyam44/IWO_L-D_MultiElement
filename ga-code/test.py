import numpy as np
import math
aoa=math.radians(16)
cfd_file = open("/home/fmg/OpenFOAM/fmg-6/run/30P30N/simvalue/CP%i/forces.dat"%1, "r")  #ADD FORCECOEFF FILE PATH HERE
MIN_Line=900
MAX_Line=1000
list=cfd_file.readlines()
#actual code
global forcex
global forcey
forcex=0
forcey=0
for i in range(MIN_Line,MAX_Line):
    forcex+=float(list[i].split('(')[2].split()[0])
    forcey+=float(list[i].split('(')[2].split()[1])
    print(forcex,forcey)

forcex/=(500)
forcey/=(500)
LDcoeffmat=np.array([[math.cos(aoa) ,-math.sin(aoa)],[math.sin(aoa),math.cos(aoa)]])
print(LDcoeffmat)
LDmat=np.dot(LDcoeffmat,np.array([forcey,forcex]).reshape(2,1))
print(LDmat)
cfd_file.close()
LDratio=LDmat[0]/LDmat[1]
print(LDratio.squeeze())

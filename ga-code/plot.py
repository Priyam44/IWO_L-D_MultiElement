import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from stl import mesh
import stl,math,copy
import STLgen as stlgen
#main element
trailing_edge_me=np.array([0.874007,0.030943])
leading_edge_me=np.array([0.043807,-0.016668])
#slat
trailing_edge_slat=np.array([0.0187595,0.00507664])

a=pd.read_csv("flap1.csv")#flap
a=a.values
b=pd.read_csv("me.csv") #main element
b=b.values
c=pd.read_csv("slat1.csv") #slat
c=c.values
c_=np.array([[0.018772,0.005079]])
c=np.append(c,c_,axis=0)

# plt.plot(c.T[0],c.T[1],'-')
# plt.plot(a.T[0],a.T[1],'-')
#gap calculation for flap
min,gap1=10000,10000
for i in range(a.shape[0]):
    min=math.pow((math.pow((a[i][0]-trailing_edge_me[0]),2)+math.pow(a[i][1]-trailing_edge_me[1],2)),0.5)
    if(gap1>min):
        gap1=min
gap1*=100
print(gap1)
# #gap calculation for slat
# min,gap2=10000,10000
# for i in range(b.shape[0]):
#     min=math.pow((math.pow((b[i][0]-trailing_edge_slat[0]),2)+math.pow(b[i][1]-trailing_edge_slat[1],2)),0.5)
#     if(gap2>min):
#         gap2=min
# gap2*=100
# print(gap2)
def cal_gap(point,profile):
    gap=10000
    p=np.around(np.array([0,0]),3)
    for i in range(profile.shape[0]):
        min=math.pow((math.pow((profile[i][0]-point[0]),2)+math.pow(profile[i][1]-point[1],2)),0.5)
        if(gap>min):
            gap=np.around(min,4)
            p=np.around(profile[i],3)
    return gap
def gap_prod(point,profile,target,target_type,tgap):
    tgap/=100
    gap=10000
    p=np.around(np.array([0,0]),3)
    for i in range(profile.shape[0]):
        min=math.pow((math.pow((profile[i][0]-point[0]),2)+math.pow(profile[i][1]-point[1],2)),0.5)
        if(gap>min):
            gap=np.around(min,4)
            p=np.around(profile[i],3)
    print("gap",gap)
    diff=abs(p-point)
    angle=math.atan2(diff[1],diff[0])
    dgap=np.around(gap-tgap,3)
    x=dgap*abs(math.cos(angle))
    y=dgap*abs(math.sin(angle))

    print("point",p)
    if(target_type=="flap"):
        target.T[0]-=x
        target.T[1]+=y
    else:
        target.T[0]-=x
        target.T[1]+=y
    return target
def transform_flap(flap_prof,angle,overhang,vert):
    overhang/=100
    vert/=100
    min,j=100000,0

    for i in range(flap_prof.shape[0]):
        if(flap_prof[i][0]<min):
            min=flap_prof[i][0]
            j=i
    fin=copy.deepcopy(flap_prof[j])
    flap_prof.T[0]-=fin[0] #rotate start
    flap_prof.T[1]-=fin[1]
    flap_prof=np.append(flap_prof,np.zeros((flap_prof.shape[0],1)),axis=1)
    flap_prof=stlgen.rotate(flap_prof,angle-30)
    flap_prof.T[0]+=trailing_edge_me[0]+overhang
    flap_prof.T[1]+=trailing_edge_me[1]-vert #rotate end
    flap_prof=np.delete(flap_prof,2,1)
    return flap_prof
def transform_slat(slat_prof,angle,overhang,vert):
    overhang/=100
    vert/=100
    max=-1000
    j=0
    for i in range(slat_prof.shape[0]):
        if(slat_prof[i][0]>max):
            max=slat_prof[i][0]
            j=i
    fin=copy.deepcopy(slat_prof[i])
    slat_prof.T[0]-=fin[0]
    slat_prof.T[1]-=fin[1]
    slat_prof=np.append(slat_prof,np.zeros((slat_prof.shape[0],1)),axis=1)
    slat_prof=stlgen.rotate(slat_prof,angle+30)
    slat_prof.T[0]+=leading_edge_me[0]+overhang
    slat_prof.T[1]+=leading_edge_me[1]-vert#rotate end
    slat_prof=np.delete(slat_prof,2,1)
    return slat_prof

#flap
#transformation

# a.T[0]-=(0.871401-0.699815*1.01)

a=transform_flap(a,11.42373996,1.09598978,3.39037871)


plt.plot(a.T[0],a.T[1],'-')

#mainelement
plt.plot(b.T[0],b.T[1],'-')
max=-1000
j=0
for i in range(b.shape[0]):
    if(b[i][0]>max):
        max=b[i][0]
        j=i
print("SS",b[j])
#slat

c=transform_slat(c,-29.41167918,-2.25,-2.2)

# print("slat",c)
# c.T[0]+=4.3*2.95/100
# c.T[1]+=1.8*2.5/100

# c.T[0]-=0.11
plt.plot(c.T[0],c.T[1],'-')
plt.axes().set_aspect("equal")
plt.show()
# import STLgen as stlgen
stlgen.STL_Gen(a.T[0],a.T[1],"flap_601.stl")
np.save("flap_601.npy",a/4)
# stlgen.STL_Gen(b.T[0],b.T[1],"me.stl")
stlgen.STL_Gen(c.T[0],c.T[1],"slat_601.stl")
np.save("slat_601.npy",c/4)
# stlgen.plot("flap.stl")
# stlgen.plot("me.stl")
# stlgen.plot("slat.stl")
# flap_mesh = mesh.Mesh.from_file("flap.stl")
# me_mesh = mesh.Mesh.from_file("me.stl")
# slat_mesh = mesh.Mesh.from_file("slat.stl")
# final= mesh.Mesh(np.concatenate([flap_mesh.data, me_mesh.data,slat_mesh.data]))
# final.save("30p30n.stl", mode=stl.Mode.ASCII)
# stlgen.plot("30p30n.stl")
# plt.show()

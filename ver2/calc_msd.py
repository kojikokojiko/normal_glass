import sys
import os
sys.path.append('/home/isobelab2022/build3/hoomd')
import itertools
import math

import gsd.hoomd
import hoomd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patches
from numba import jit, f8, i8, b1, void,njit
import matplotlib.patches as pat
import sys
import os
import fresnel
from PIL import Image
import IPython
import packaging.version
import random
import time
 


nu=float(sys.argv[1])
fixed_percent=float(sys.argv[2])
L=64
 



ver=str(nu)+"_"+str(fixed_percent)
main_dir="./"+ver


traj = gsd.hoomd.open('./'+ver+'/log_pos_'+ver+'.gsd', 'rb')

print(traj[0].log)



pos=[traj[i].particles.position for i in range(len(traj)-1)]

pos=np.array(pos)
print(len(pos))

print(pos.shape)


pos=np.transpose(pos,(2,1,0))

print(pos.shape)

plt.plot(pos[0][63])
plt.savefig(main_dir+"/pos1.png")
plt.clf()
plt.cla()
# plt.clear()

rx=pos[0]
ry=pos[1]
# rz=pos[2]


def adjust_periodic(x,L):
    for t in range(len(x)-1):
        if x[t+1] - x[t] > L/2:
            x[t+1] -= (x[t+1] - x[t]+L/2)//L*L
        if x[t+1] - x[t] < -L/2:
            x[t+1] += (x[t] - x[t+1]+L/2)//L*L

print(rx.shape)

# (N,T)=rx.shape

for i in range(len(rx)):
    adjust_periodic(rx[i],64)
    adjust_periodic(ry[i],64)
    # adjust_periodic(rz[i],64)


def calc_msd_simple(rx,ry):
    n = rx.shape[0]
    T=rx.shape[1]

    msd = []
    for t in range(1,T):
        r2 = 0.0
        for i in range(n-t):
            dx=rx[i+t]-rx[i]
            dy=ry[i+t]-ry[i]
            r2 += dx*dx+dy*dy
        x2 /= (n-t)
        msd.append(x2)
    return msd


N = rx.shape[0]
T=rx.shape[1]

msd_list = []
for t in range(T):
    m2=0.0

    for i in range(N):
        dx=rx[i][t]-rx[i][0]
        dy=ry[i][t]-ry[i][0]
        m2 += dx*dx+dy*dy
    m2 /= N

    msd_list.append(m2)


plt.plot(msd_list)
plt.savefig(main_dir+"/msd_list.png")
plt.cla()
msd_list=np.array(msd_list)
np.save(main_dir+"/msd_list.npy",msd_list)





# lambda_rx=map(adjust_periodic,rx)
# lambda_ry=map(adjust_periodic,ry)
               
plt.plot(pos[0][63])
plt.savefig(main_dir+"/pos2.png")


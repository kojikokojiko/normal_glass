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
import pickle


nu=float(sys.argv[1])
fixed_percent=float(sys.argv[2])
kbT=float(sys.argv[3])
L=64
 



media_dir="/media/isobelab2022/data/normal_glass/ver2"



ver=str(nu)+"_"+str(fixed_percent)+"_"+str(kbT)
main_dir="./"+ver

traj_dir=media_dir+"/"+ver


van_hove_path=traj_dir+"/van_hove_list.npz"

van_hove_npz=np.load(van_hove_path)
van_hove_list=van_hove_npz["van_hove_list"]
zerod_time=van_hove_npz["zerod_time"]
bins_list=van_hove_npz["bins_list"]





watch_around_time_list=[20,30,40]

watch_index_list=[]
for watch_around_time in watch_around_time_list:
    watch_index=(np.abs(zerod_time - watch_around_time)).argmin()
    watch_index_list.append(watch_index)

print(watch_index_list)

plt.figure()
for watch_index in watch_index_list:
    delta_t=zerod_time[watch_index]

    van_hove=van_hove_list[watch_index] 
    bins=bins_list[watch_index]
    plt.plot(bins, van_hove)


plt.xscale("log")
plt.xlabel('delta_r')
plt.ylabel('P(delta_r,delta_t)')
plt.title('van hove plot')
plt.legend()

plt.savefig("van_hove_analy_{0}kbT_{1}nu.png".format(kbT,nu))







print("FiNISH")
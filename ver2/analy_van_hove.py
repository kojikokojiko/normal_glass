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


van_hove_path=traj_dir+"/van_hove.npz"

van_hove_npz=np.load(van_hove_path)
print(van_hove_npz.files)
van_hove_list=van_hove_npz["van_hove_list"]
zerod_time=van_hove_npz["zerod_time"]
bins=van_hove_npz["bins"]




print(zerod_time.shape)
print(van_hove_list.shape)
print(bins.shape)


print(bins[0])

bins_list=np.array(bins)
watch_around_time_list=[0.1,1.0,10,100,300,500]

watch_index_list=[]
for watch_around_time in watch_around_time_list:
    watch_index=(np.abs(zerod_time - watch_around_time)).argmin()
    watch_index_list.append(watch_index)

print(watch_index_list)

plt.figure()
for watch_index in watch_index_list:
    delta_t=zerod_time[watch_index]

    van_hove=van_hove_list[watch_index] 
    # ビンの中心を計算
    bins_centers = (bins[:-1] + bins[1:]) / 2
    # print(bins_centers)


    plt.plot(bins_centers, van_hove,label="delta_t={0}".format(delta_t))

    # plt.hist(van_hove,label="delta_t={0}".format(delta_t),bins=10000)


plt.xlabel('delta_r')
plt.yscale("log")
plt.ylabel('P(delta_r,delta_t)')
plt.title('van hove plot(kbT={0},nu={1})'.format(kbT,nu))
plt.legend()

plt.savefig(main_dir+"/van_hove_analy_{0}kbT_{1}nu.png".format(kbT,nu))







print("FiNISH")
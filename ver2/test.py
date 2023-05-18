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

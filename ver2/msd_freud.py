import sys
import os
sys.path.append('/home/isobelab2022/build3/hoomd')

import matplotlib.pyplot as plt
import matplotlib.patches as pat
import gsd.hoomd
import hoomd
import sys
from PIL import Image
import os
import math
import  numpy as np
import seaborn as sns
import sys
import math
import freud

nu=float(sys.argv[1])
pl=float(sys.argv[2])
fixed_percent=float(sys.argv[3])
ver=str(nu)+"_"+str(pl)+"_"+str(fixed_percent)


main_dir="./"+ver
temp_dir="./"+ver+"/log_pos_"+ver+".gsd"
traj = gsd.hoomd.open(temp_dir, 'rb')

pos=[]
for frame in traj:
    pos.append(frame.particles.position)
# pos=traj.particles.position


msd_analyzer=freud.msd.MSD(mode='direct')
msd_analyzer.compute(pos,reset=False)
print(msd_analyzer.msd)
print(len(msd_analyzer.msd)) 
msd=msd_analyzer.msd
plt.plot(msd)
plt.savefig(main_dir+"/msd.png")
msd_analyzer.plot()
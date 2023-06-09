import numpy as np
import matplotlib.pyplot as plt

nu_list=[0.72,0.74,0.76,0.78 ,0.8]
fixed_percent=0.0
kbT=1.0
# kbT_list=[0.1,0.2,0.3,0.4,0.5,0.6,0.7]


plt.figure()
media_dir="/media/isobelab2022/data/normal_glass/ver2"

for nu in nu_list:

    
    ver=str(nu)+"_"+str(fixed_percent)+"_"+str(kbT)
    main_dir=media_dir="/media/isobelab2022/data/normal_glass/ver2"+"/"+ver


    msd_npz=np.load(main_dir+"/msd.npz")
    msd_list=msd_npz["msd_list"]
    zerod_time_list=msd_npz["zerod_time"]

    plt.plot(zerod_time_list,msd_list,label="nu="+str(nu))


plt.xscale("log")
plt.yscale("log")
plt.xlabel('t')
plt.ylabel('MSD')
plt.title('MSD plot')
plt.legend()

plt.savefig("msd_analy_{0}kbT.png".format(kbT))


import numpy as np
import matplotlib.pyplot as plt

nu_list=[0.7 ,0.72,0.74,0.76,0.78 ]
fixed_percent=0.0
# kbT_list=[0.1,0.2,0.3,0.4,0.5,0.6,0.7]


plt.figure()

for nu in nu_list:

    
    ver=str(nu)+"_"+str(fixed_percent)
    main_dir="./"+ver


    msd_list=np.load(main_dir+"/msd_list.npy")
    plt.plot(msd_list,label="nu="+str(nu))


plt.xscale("log")
plt.yscale("log")
plt.xlabel('t')
plt.ylabel('MSD')
plt.title('MSD plot')
plt.legend()

plt.savefig("msd＿analy.png")


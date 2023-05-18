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





nu=float(sys.argv[1])
pl=float(sys.argv[2])

ver=str(nu)+"_"+str(pl)

main_dir="./"+ver


# temp_dir="./"+ver+"/log_pos_"+ver+".gsd"
temp_dir="./"+ver+"/log_pos_"+ver+".gsd"
figuredir=main_dir+"/figure_scat"
# i_phi_dir=main_dir+"/i_phi"
if not os.path.exists(figuredir): os.makedirs(figuredir)
# if not os.path.exists(i_phi_dir): os.makedirs(i_phi_dir)


print(temp_dir)
traj = gsd.hoomd.open(temp_dir, 'rb')
# traj = gsd.hoomd.open(dir, 'rb')




lx=traj[0].configuration.box[0]
ly=traj[0].configuration.box[1]
ppi=72 # points per inche 
figsize=(10,10*ly/lx)
plt.figure(figsize=figsize)

NP=len(traj[0].particles.position)
print(NP)

diameter=traj[0].particles.diameter


small_index=np.where(diameter==min(diameter))
large_index=np.where(diameter==max(diameter))
# あんま頭よくない処理
size_list=diameter/2.0


bo=(lx/0.4)
l_gx_d=lx/bo
l_gx=l_gx_d
l_gy=l_gx
n_gx=math.ceil(lx/l_gx)
n_gy=math.ceil(ly/l_gy)
n_g_len=n_gx*n_gy






output_list=[]
phi6_output_list=[]
# 描画設定##############################################
# https://qiita.com/stanaka2/items/c40841f858d7083aad4e
fig = plt.figure(figsize=figsize, dpi=100.0)
ax  = fig.add_axes((0.1,  0.1, 0.8, 0.8))

# 枠の範囲指定
xmin,xmax=-lx/2.0, lx/2.0
ymin,ymax=-ly/2.0, ly/2.0

# 描写範囲の長さを取得(dpi単位)
# x軸をベースに計算しているがy軸でも同じ。アスペクト比が違えばおかしくなる
ax_length=ax.bbox.get_points()[1][1]-ax.bbox.get_points()[0][1]

# dpi単位の長さをポイント単位に変換(dpiからインチ単位にし、インチからポイント単位に変換)
ax_point = ax_length*ppi/fig.dpi

# x軸の実スケールをポイント単位に変換
xsize=xmax-xmin
fact=ax_point/xsize

# scatterのマーカーサイズは直径のポイントの二乗を描くため、実スケールの半径をポイントに変換し直径にしておく
size_list*=2*fact

####################################################

# data=np.load(main_dir+"/phi6_2.npy")
data_index=0
for t in range(len(traj)-1,0,-2):

    print(t)
    pos=traj[t].particles.position.T
    rx=pos[0]
    ry=pos[1]


    # extract_index1=np.where(phi6_2>0.0)
    # extract_phi61=phi6_2[extract_index1]
    # extract_rx1=rx[extract_index1]
    # extract_ry1=ry[extract_index1]


    # extract_index2=np.where(0.2<phi6_2)
    # extract_phi62=phi6_2[extract_index2]
    # extract_rx2=rx[extract_index2]
    # extract_ry2=ry[extract_index2]


    # extract_index3=np.where(0.55<phi6_2)
    # extract_phi63=phi6_2[extract_index3]
    # extract_rx3=rx[extract_index3]
    # extract_ry3=ry[extract_index3]

    # extract_index4=np.where(0.6<phi6_2)
    # extract_phi64=phi6_2[extract_index4]
    # extract_rx4=rx[extract_index4]
    # extract_ry4=ry[extract_index4]


    # extract_index5=np.where(0.85<phi6_2)
    # extract_phi65=phi6_2[extract_index5]
    # extract_rx5=rx[extract_index5]
    # extract_ry5=ry[extract_index5]
    
# 描画############################
# Real_phi6#################################
# dpiは任意の値を設定できる。デフォルトは100
    fig = plt.figure(figsize=figsize, dpi=300.0)
    ax  = fig.add_axes((0.1,  0.1, 0.8, 0.8))

    ax.set_aspect("equal") # アスペクト比を等しくする
                        # 今は枠のサイズ設定で最初から等しいのであってもなくてもよい

    # 二乗にして与える 
    # im=ax.scatter(rx,ry,s=size**2,c=phi6_2,cmap="cool", linewidths=0)
    ax.scatter(rx[small_index],ry[small_index],s=size_list[small_index]**2,c="r", linewidths=1,ec="g")
    ax.scatter(rx[large_index],ry[large_index],s=size_list[large_index]**2,c="b", linewidths=1,ec="g") 
    # im1=ax.scatter(extract_rx1,extract_ry1,s=size**2,c="blue", linewidths=0)
    # im2=ax.scatter(extract_rx2,extract_ry2,s=size**2,c="c", linewidths=0)
    # im3=ax.scatter(extract_rx3,extract_ry3,s=size**2,c="m", linewidths=0)
    # im4=ax.scatter(extract_rx4,extract_ry4,s=size**2,c="g", linewidths=0)
    # im5=ax.scatter(extract_rx5,extract_ry5,s=size**2,c="red", linewidths=0)

    # Blue 0.2
    # Cean 0.2<0.55
    # Mazenta 0.55<0.6
    # green 0.6<0.85
    # red　0.85> 

    # for i in range(0,NP,50):
    #     c=pat.Circle(xy=(rx[i],ry[i]),radius=r_sann_1[i],fc="g",alpha=0.2)
    #     ax.add_patch(c)

    # for i in range(100,200):
    #     ax.text(rx[i], ry[i], i,fontsize="small")

    # c=pat.Circle(xy=(lx/4-lx/2,ly/2-ly/2),radius=static_dia/2,fc="b")
    # ax.add_patch(c) 
    # plt.colorbar(im,ticks=[0.0,0.2,0.4 ,0.6,0.75,0.8])

    ax.set_xlim(xmin,xmax)
    ax.set_ylim(ymin,ymax)

    ax.set_xticks(np.linspace(ax.get_xlim()[0], ax.get_xlim()[1], 5))
    ax.set_yticks(np.linspace(ax.get_ylim()[0], ax.get_ylim()[1], 5))

    # ax.grid(which='both', axis='both')

    # plt.colorbar()
    # plt.colorbar(sc, ax=ax)
    plt.title(t)
    plt.savefig(figuredir+"/snap_{0}.png".format(t))
    # plt.savefig(main_dir+"/figure.png")
    # plt.cla()
    plt.close()





# # Real_phi6#################################
#     plt.scatter(rx,ry,s=size**2,c=r_phi6,cmap="seismic", linewidths=0,vmin=-1.0,vmax=1.0)
#     plt.colorbar()
#     plt.savefig(r_phi_dir+"/r_phi_color{0}.png".format(t))
#     # plt.clf()
#     # plt.cla()
#     plt.close()


# output_list=np.array(output_list)
# np.save(main_dir+"/phi6_bin2.npy", output_list)

# phi6_output_list=np.array(phi6_output_list)
# np.save(main_dir+"/phi6_2.npy", phi6_output_list)






# # Imagenary_phi6#################################
# # dpiは任意の値を設定できる。デフォルトは100
#     fig = plt.figure(figsize=(12.5,5.0), dpi=100.0)
#     ax  = fig.add_axes((0.1,  0.1, 0.8, 0.8))

#     ax.set_aspect("equal") # アスペクト比を等しくする
#                         # 今は枠のサイズ設定で最初から等しいのであってもなくてもよい

#         # 二乗にして与える 
#     ax.scatter(rx,ry,s=size**2,c=i_phi6,cmap="seismic", linewidths=0,vmin=-1.0,vmax=1.0)
#     # plt.colorbar()

#     ax.set_xlim(xmin,xmax)
#     ax.set_ylim(ymin,ymax)

#     ax.set_xticks(np.linspace(ax.get_xlim()[0], ax.get_xlim()[1], 5))
#     ax.set_yticks(np.linspace(ax.get_ylim()[0], ax.get_ylim()[1], 5))

#     # ax.grid(which='both', axis='both')

#     # plt.colorbar()
#     # plt.colorbar(sc, ax=ax)
#     plt.title("i_phi6")
#     plt.savefig(i_phi_dir+"/i_phi{0}.png".format(t))

#     # plt.cla()
#     plt.close()

# # get cmap############################################################
#     plt.scatter(rx,ry,s=size**2,c=i_phi6,cmap="seismic", linewidths=0,vmin=-1.0,vmax=1.0)
#     plt.colorbar()
#     plt.savefig(i_phi_dir+"/i_phi_color{0}.png".format(t))
#     # plt.cla()
#     plt.close()


####################################################################


    # local_phi6=0.0
    # for i in range(N):
    #     local_phi6+=math.sqrt(r_phi6[i]**2+i_phi6[i]**2)
    # local_phi6/=N

    # mean_re=np.mean(r_phi6)
    # mean_im=np.mean(i_phi6)
    # global_phi6=math.sqrt(mean_re**2+mean_im**2)

    # global_phi6_list[t]=global_phi6
    # local_phi6_list[t]=local_phi6



# def sun_cut(NP,r_cut_sann,rx,ry,lx,ly,pair_list_g):
#     pair_length_sann=50
#     dr_list_sann=np.zeros(NP,pair_length_sann)
#     pair_index_sann=np.full((NP,pair_length_sann),-1)
#     for i in range(NP):
#         pair_count=0
#         loop_length=pair_list_g[i][-1]

#         for loop_count in range(loop_length):
#             pair_index=G_MAP[i][loop_count]
#             rxji=rx[j]-rx[i]
#             ryji=ry[j]-ry[i]
#             if(rxji>=lx/2):
#                 rxji=rxji-lx
#             elif (rxji<=-lx/2):
#                 rxji=rxji+lx

#             if(ryji>=ly/2):
#                 ryji=ryji-ly
#             elif (ryji<=-ly/2):
#                 rxji=ryij+ly
#             r2=rxji*rxji+ryij*ryij
#             if(r2<=r_cut_sann**2):
#                 pair_index_sann[i][k]=pair_index
#                 dr_list_sann[i][k]=np.sqrt(r2)
#                 pair_count+=1
#         dr_list_sann[i][-1]=pair_count

            
#     return dr_list_sann,pair_index_sann
    
# def sann(dr_list_sann,pair_index_sann,NP):
#     N_MAX=50
#     NN_NUM=np.zeros(N_MAX)
#     NN=np.zeros(NP,N_MAX)

#     for i in range(NP):
#         num_pair=pair_index_sann[i][-1]
#         for j in range(3,num_pair):
#             SANN_M=j
#             r_min=dr_list_sann[i][SANN_M]
#             r_max=dr_list_sann[i][num_pair]
#             for n in range(N_MAX):
#                 r_mid=(r_min+r_max)/2
#                 sum_rm_min=0.0
#                 sum_rm_max=0.0
#                 sum_rm_mid=0.0
#                 for m in range(SANN_M):
#                     sum_rm_min+=math.acos(dr_list_sann[i][m]/r_min)
#                     sum_rm_max+=math.acos(dr_list_sann[i][m]/r_max)
#                     sum_rm_mid+=math.acos(dr_list_sann[i][m]/r_mid)
#                 f_min=sum_rm_min-math.pi
#                 f_max=sum_rm_max-math.pi
#                 f_mid=sum_rm_max-math.pi
#                 if(abs(f_mid)<=1e-10):
#                     r_sann_m=dr_mid
#                     if(r_sann_m<=dr_list_sann[i][SANN_M+1]):
#                         NN_NUM[i]=SANN_M
#                         for k in range(NN_NUM[i]):
#                             NN[i][k]=NN_S[i][k]
                        
#                 if(f_mid*f_min>=0.0):
#                     r_min=r_mid
#                     f_min=f_mid
#                 else:
#                     r_max=r_mid
#                     f_max=f_mid







# def sann_sort(r_tmp,dr,NN_S,NN_NUM_S):
#     for i in range(N):
#         sr_num=int(NN_NUM_S[i])
#         for j in range(NN_MAX):
#             r_tmp=dr[i][j]
#             i_tmp=NN_S[i][j]
#         im_inssrot()
#         for j in range(NN_MAX):
#             dr[i][j]=R_





# def make_neighbor_list(around_grid_num):
#     if around_grid_num==5:
#         neighbor_col=[
#             -5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5,
#             -5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5,
#             -5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5,
#             -5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5,
#             -5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5,
#             -5,-4,-3,-2,-1,    1, 2, 3, 4, 5,
#             -5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5,
#             -5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5,
#             -5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5,
#             -5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5,
#             -5,-4,-3,-2,-1, 0, 1, 2, 3, 4, 5
#         ]
#         neighbor_row=[
#             -5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,
#             -4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,
#             -3,-3,-3,-3,-3,-3,-3,-3,-3,-3,-3,
#             -2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,
#             -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
#              0, 0, 0, 0, 0,    0, 0, 0, 0, 0,
#              1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
#              2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
#              3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
#              4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
#              5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5

#         ]


#     # neighbor_row=[]
#     # neighbor_row_value=-around_grid_num
#     # for i in range(around_grid_num+1):
#     #     for j in range(around_grid_num):
#     #         neighbor_row.append(neighbor_row_value)
#     #     neighbor_row_value+=1
    
#     # for i in range(around_grid_num):
#     #     for j in  range(around_grid_num+1):
#     #         neighbor_row.append(neighbor_row_value)
#     #     neighbor_row_value+=1
    
#     # # yの近接リスト作成
#     # neighbor_col=[]
#     # neighbor_col_value=-around_grid_num
#     # for i in range(0,around_grid_num+1):
#     #     for j in range(around_grid_num):
#     #         neighbor_col.append(neighbor_col_value)
#     #         neighbor_col_value+=1
#     #     neighbor_col_value=-around_grid_num

#     # for i in range(around_grid_num):
#     #     for j in range(around_grid_num+1):
#     #         neighbor_col.append(neighbor_col_value)
#     #         neighbor_col_value+=1
#     #     neighbor_col_value=-around_grid_num

#     return neighbor_row,neighbor_col
    

# def make_pairlist_g(pair_list_g,select_gx,select_gy,G_MAP,select_index,NEIGHBOR_LEN,NEIGHBOR_COL,NEIGHBOR_ROW):
#     particle_counter=0
#     for k in range(NEIGHBOR_LEN):
#         search_gx=select_gx+NEIGHBOR_COL[k]
#         search_gy=select_gy+NEIGHBOR_ROW[k]
#         if(search_gx>=N_GX):
#             search_gx-=N_GX
#         elif(search_gx<0):
#             search_gx+=N_GX
#         if (search_gy>=N_GY):
#             search_gy-=N_GY
#         elif(search_gy<0):
#             search_gy+=N_GY
        
#         search_index=G_MAP[search_gy][search_gx]
#         pair_list_g[select_index][particle_counter]=search_index

#         particle_counter+=1

#     pair_list_g[select_index][-1]=particle_counter





    


# def gmap_create(l_gx,l_gy,pair_length,NP,N_GX,N_GY,rx,ry,neighbor_row,neighbor_col):
#     G_MAP=np.full((N_GY,N_GX),-1)
#     pair_list_g=np.full(pair_length,-1)
#     for i in range(NP):
#         gx_map=int(rx[i]/l_gx)
#         gy_map=int(ry[i]/l_gy)
#         G_MAP[gy_map][gx_map]=i
#     for i in range(N_GY):
#         for j in range(N_GX):
#             select_index=G_MAP[i][j]
#             make_pairlist_g()

#     return pair_list_g

#!/bin/bash



nu_list=(  0.72 0.74 0.76  0.78 0.8 )
# 0.55 0.6 0.65 0.7 0.75 0.8)
fixed_per_list=(0.0)
kbT_list=(1.0 0.5 )

for nu in "${nu_list[@]}"
do 
    for fixed_per in  "${fixed_per_list[@]}"
    do
        for kbT in  "${kbT_list[@]}"
        do
            echo $nu   $fixed_per $kbT
            name="${nu}_${fixed_per}_${kbT}"
            # nohup time python normal_glass.py $nu  $fixed_per $kbT> ${name}.log 2>&1 &
            # nohup time python animation.py $nu  $fixed_per $kbT> ${name}.log 2>&1 &
            # nohup time python calc_msd.py $nu  $fixed_per $kbT> ${name}.log 2>&1 &
            # nohup time python van_hove.py $nu  $fixed_per $kbT> ${name}.log 2>&1 &
            nohup time python analy_van_hove.py $nu  $fixed_per $kbT> ${name}.log 2>&1 &
        done
    done
done

# g++ -std=c++11 -o iabp2 iABP2.cpp
# nohup time ./iabp2 1.0 1.0 10000 # Pe M N
	






# rho=float(sys.argv[1])
# ave_flow=float(sys.argv[2])
# static_dia=float(sys.argv[3])
# reduced_speed=float(sys.argv[5])
# rotational_diffusion=float(sys.argv[6])


#!/bin/bash

set -e

gamma=5e12
N=12

for dt in 10 20
do
	octave -q --eval "diffusionToNativeData2D('/home/other/roald/clusterdata/${gamma}_${N}_${dt}/*/*strand*', 'data/strand_dt${dt}_wait50_time1000_N${N}_g${gamma}');" &
	octave -q --eval "diffusionToNativeData2D('/home/other/roald/clusterdata/${gamma}_${N}_${dt}/*/*particle*', 'data/particle_dt${dt}_wait50_time1000_N${N}_g${gamma}');" &
	wait
done





exit

# for gamma in 3e12 5e12; do
# 	for N in 30 50; do
# 		echo $gamma $N
# 		octave -q --eval "diffusionToNativeData2D('/home/other/roald/clusterdata/${gamma}_${N}/*strand*', 'data/strand_dt40_wait50_time1000_N${N}_g${gamma}');" &
# 		octave -q --eval "diffusionToNativeData2D('/home/other/roald/clusterdata/${gamma}_${N}/*particle*', 'data/particle_dt40_wait50_time1000_N${N}_g${gamma}');" &
# 		wait
# 	done
# done


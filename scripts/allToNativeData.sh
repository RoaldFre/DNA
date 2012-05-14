#!/bin/bash

set -e

for gamma in 3e12 5e12; do
	for N in 30 50; do
		echo $gamma $N
		octave -q --eval "toNativeData2D('/home/other/roald/clusterdata/${gamma}_${N}/*strand*', 'data/strand_dt40_wait50_time1000_N${N}_g${gamma}');" &
		octave -q --eval "toNativeData2D('/home/other/roald/clusterdata/${gamma}_${N}/*particle*', 'data/particle_dt40_wait50_time1000_N${N}_g${gamma}');" &
		wait
	done
done


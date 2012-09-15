#!/bin/bash

set -e

gamma=5e12
N=12
dt=15
wait=100
time=1000
T=296K

for gamma in 4e12 5e12 6e12
do
	dir="$HOME/clusterdata/diffusion/dt${dt}_wait${wait}_time${time}_N${N}_T${T}/g${gamma}/"
	octave -q --eval "diffusionToNativeData2D('$dir/*', '$dir/combinedData');"
done


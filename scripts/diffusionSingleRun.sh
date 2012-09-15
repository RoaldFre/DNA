#!/bin/sh

set -e

destinationDirRoot=$1
gamma=$2
N=$3
timestep=$4
time=$5
wait=$6
suffix=$7 # job label when running on cluster

set -e

main=diffusion #Should be in current working dir, or in $PATH
interval=10
sizeFactor=1e10 #So we don't have problems with periodic boundary condition!
boxes="-b1" #So we don't grind ourselves to a halt with the above (huge) worldsize!

T="296K"

integrator=l

#render="-r -R 3"
fps=10


tempOutputFile=`mktemp`
destinationDir="$destinationDirRoot/dt${timestep}_wait${wait}_time${time}_N${N}_T${T}/g${gamma}/"
destination="$destinationDir/diff_${suffix}"

mkdir -p $destinationDir


seq=`perl -e "print 'A'x$N"`
$main -i $integrator $render -f $fps -t $timestep -T $T -I $interval -P $time -W $wait -s $seq -D $tempOutputFile -S $sizeFactor -g $gamma $boxes
mv ${tempOutputFile} ${destination}


#!/bin/sh

destinationDir=$1
gamma=$2
N=$3
runNumber=$4

set -e

main=diffusion #Should be in current working dir, or in $PATH
timestep=40
interval=10
time=100
wait=10
sizeFactor=1e10 #So we don't have problems with periodic boundary condition!
boxes="-b1" #So we don't grind ourselves to a halt with the above (huge) worldsize!

T=300

integrator=l

#render="-r -R 3"
fps=10


outputBaseName=`mktemp`
nameSuffix="dt${timestep}_wait${wait}_time${time}_N${N}_g${gamma}_${runNumber}"
particleFilenameBase="$destinationDir/diff_particle_${nameSuffix}"
strandFilenameBase="$destinationDir/diff_strand_${nameSuffix}"


seq=`perl -e "print 'A'x$N"`
$main -i $integrator $render -f $fps -t $timestep -T $T -I $interval -P $time -W $wait -s $seq -D $outputBaseName -S $sizeFactor -g $gamma $boxes
mv ${outputBaseName}_particle ${particleFilenameBase}
mv ${outputBaseName}_strand ${strandFilenameBase}


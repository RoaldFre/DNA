#!/bin/sh

set -e

destinationDirRoot=$1
N=$2
sampleTemp=$3 #XXX given to -E as well, but wait=0 so shouldn't matter
startTemp=$4
stepTemp=$5
nSteps=$6
measureTime=$7
allowUnbounded=$8
suffix=$9 # job label when running on cluster

set -e

main=hairpin # Should be in current working dir, or in $PATH
timestep=15
interval=10
time=1450
wait=0

fps=10
#render="-r -R 3 -f $fps"


outputBaseDir=`mktemp -d`
outputBaseName="$outputBaseDir/outputFrom_${suffix}"
destinationDir="$destinationDirRoot/N${N}_dt${timestep}_time${time}_measTime${measureTime}_Tsample${sampleTemp}_Tstart${startTemp}_Tstep${stepTemp}_nSteps${nSteps}"
nameSuffix="wait${wait}_interval${interval}_allowUnb${allowUnbounded}_${suffix}"
basePairFilename="$destinationDir/basePair_${nameSuffix}"
hairpinFilename="$destinationDir/hairpin_${nameSuffix}"

mkdir -p $destinationDir

seq=`perl -e "print 'C'x$N,'G'x$N"`
$main -t $timestep -E $sampleTemp -T $sampleTemp -I $interval -P $time -W $wait -s $seq -D $outputBaseName -A $startTemp -B $stepTemp -C $nSteps -G $measureTime -H $allowUnbounded
mv ${outputBaseName}_basePair ${basePairFilename}
mv ${outputBaseName}_hairpin ${hairpinFilename}

rm -rf $outputBaseDir


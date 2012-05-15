#!/bin/sh

set -e

destinationDirRoot=$1
N=$2
relaxTemp=$3
sampleTemp=$4
startTemp=$5
stepTemp=$6
nSteps=$7
measureTime=$8
allowedUnbound=$9
suffix=${10} # job label when running on cluster

set -e

main=hairpin # Should be in current working dir, or in $PATH
timestep=20
interval=10
time=1000
wait=50

fps=10
#render="-r -R 3 -f $fps"


outputBaseDir=`mktemp -d`
outputBaseName="$outputBaseDir/outputFrom_${suffix}"
destinationDir="$destinationDirRoot/N${N}_dt${timestep}_time${time}_measTime${measureTime}_Tsample${sampleTemp}_Tstart${startTemp}_Tstep${stepTemp}_nSteps${nSteps}"
nameSuffix="wait${wait}_interval${interval}_allowUnb${allowUnbounded}_${suffix}"
basePairFilename="$destinationDir/basePair_${nameSuffix}"
hairpinFilename="$destinationDir/hairpin_${nameSuffix}"

mkdir -p $destinationDir

seq=`perl -e "print 'C'x$N,Gx$N"`
$main -t $timestep -E $relaxTemp -T $sampleTemp -I $interval -P $time -W $wait -s $seq -D $outputBaseName 
mv ${outputBaseName}_basePair ${basePairFilename}
mv ${outputBaseName}_hairpin ${hairpinFilename}

rm -rf $outputBaseDir


#!/bin/sh

set -e

destinationDirRoot=$1
N=$2
startTemp=$3
stepTemp=$4
nSteps=$5
measureTime=$6
relaxTime=$7
suffix=$8 # job label when running on cluster

set -e

main=hairpin # Should be in current working dir, or in $PATH
timestep=15
interval=10
wait=0

fps=10
#render="-r -R 3 -f $fps"


outputBaseDir=`mktemp -d`
outputFile="$outputBaseDir/outputFrom_${suffix}"
destinationDir="$destinationDirRoot/meltingTemp_dt${timestep}_measTime${measureTime}_relaxTime${relaxTime}_Tstart${startTemp}_Tstep${stepTemp}_nSteps${nSteps}/N${N}"
destinationFile="$destinationDir/meltingTemp_${suffix}"

mkdir -p $destinationDir

seq=`perl -e "print 'C'x$N,'G'x$N"`
$main -t $timestep -I $interval -W $wait -s $seq -D $outputFile -A $startTemp -B $stepTemp -C $nSteps -G $measureTime -L $relaxTime -X m
mv ${outputFile} ${destinationFile}

rm -rf $outputBaseDir

#!/bin/sh

set -e

destinationDirRoot=$1
seq=$2
startTemp=$3
stepTemp=$4
nSteps=$5
measureTime=$6
relaxTime=$7
seqDescription=$8 # used as directory name
suffix=$9 # job label when running on cluster

set -e

main=hairpin # Should be in current working dir, or in $PATH
timestep=15
interval=10
wait=0

saltConcentration=115

fps=10
#render="-r -R 3 -f $fps"


outputBaseDir=`mktemp -d`
outputFile="$outputBaseDir/outputFrom_${suffix}"
destinationDir="$destinationDirRoot/meltingTemp_dt${timestep}_measTime${measureTime}_relaxTime${relaxTime}_Tstart${startTemp}_Tstep${stepTemp}_nSteps${nSteps}/${seqDescription}"
destinationFile="$destinationDir/meltingTemp_${suffix}"

mkdir -p $destinationDir

$main -t $timestep -I $interval -W $wait -s $seq -D $outputFile -A $startTemp -B $stepTemp -C $nSteps -G $measureTime -L $relaxTime -N $saltConcentration -X m
mv ${outputFile} ${destinationFile}

rm -rf $outputBaseDir


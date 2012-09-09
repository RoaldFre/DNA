#!/bin/sh

set -e

destinationDirRoot=$1
seq=$2
startTemp=$3
stepTemp=$4
nSteps=$5
measureTime=$6
relaxTime=$7
wait=$8
saltConcentration=$9
seqDescription=${10} # used as directory name
suffix=${11} # job label when running on cluster

set -e

main=hairpin # Should be in current working dir, or in $PATH
timestep=15
interval=20

fps=10
#render="-r -R 3 -f $fps"

verbose="-V"


outputBaseDir=`mktemp -d`
outputFile="$outputBaseDir/outputFrom_${suffix}"
destinationDir="$destinationDirRoot/meltingTemp_dt${timestep}_wait${wait}_measTime${measureTime}_relaxTime${relaxTime}_Tstart${startTemp}_Tstep${stepTemp}_nSteps${nSteps}_salt${saltConcentration}/${seqDescription}"
destinationFile="$destinationDir/meltingTemp_${suffix}"

mkdir -p $destinationDir

$main -t $timestep -I $interval -W $wait -s $seq -D $outputFile -T $startTemp -A $startTemp -B $stepTemp -C $nSteps -G $measureTime -L $relaxTime -N $saltConcentration -X m $verbose
mv ${outputFile} ${destinationFile}

rm -rf $outputBaseDir



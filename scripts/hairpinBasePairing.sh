#!/bin/sh

# just dumps the base pairing info
# Useful to compute the melting temperature yourself

set -e

destinationDirRoot=$1
seq=$2
temperature=$3
time=$4
description=$5 # description of sequence
suffix=$6 # job label when running on cluster

saltConcentration=100

main=hairpin # Should be in current working dir, or in $PATH
timestep=15
interval=100

fps=10
#render="-r -R 3 -f $fps"


outputBaseDir=`mktemp -d`
outputFile="$outputBaseDir/outputFrom_${suffix}"
destinationDir="$destinationDirRoot/dt${timestep}_time${time}_salt${saltConcentration}/$description/T${temperature}"
destinationFile="$destinationDir/meltingTemp_${suffix}"

mkdir -p $destinationDir

$main -t $timestep -I $interval -s $seq -D $outputFile -T $temperature -N $saltConcentration -P $time -p -e -k
mv ${outputFile}_basePairing ${destinationFile}_basePairing
mv ${outputFile}_endToEnd ${destinationFile}_endToEnd
mv ${outputFile}_temperature ${destinationFile}_temperature

rm -rf $outputBaseDir


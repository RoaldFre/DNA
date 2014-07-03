#!/bin/bash

scriptsdir="$HOME/DNA/scripts"

set -e

destinationDirRoot=$1
seq=$2
K=$3
Rref=$4
temperature=$5
time=$6
gamma=$7
suffix=$8 # job label when running on cluster

main=hairpin # Should be in current working dir, or in $PATH
timestep=12
interval=5
wait=0

echo "Making temporary directory"
outputBaseDir=`mktemp -d`
echo "Made temporary directory"

outputFile="$outputBaseDir/outputFrom_${suffix}"
destinationDir="$destinationDirRoot/$seq/dt${timestep}/T${temperature}/K${K}/gamma${gamma}_time${time}/"
destinationFile="$destinationDir/endToEnd_${suffix}"

echo "Making destination directory"
mkdir -p $destinationDir
echo "Made destination directory"

echo "Starting main program"
echo ""
$main -t $timestep -T $temperature -I $interval -P $time -W $wait -s $seq -D $outputFile -u ${K}:${Rref} -g $gamma -e -p
echo ""
echo "End of main program"

echo "Moving output files"
mv "${outputFile}_endToEnd" "${destinationFile}"
mv "${outputFile}_basePairing" "${destinationFile}_basePairing"

echo "Converting to native Octave format"
$scriptsdir/eteFreeEnergyToNative.sh "$destinationFile"

echo "Removing temporary directory"
rm -rf $outputBaseDir

echo "All done! Have an awesome day!"


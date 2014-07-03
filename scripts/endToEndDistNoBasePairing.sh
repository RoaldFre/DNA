#!/bin/bash

scriptsdir="$HOME/DNA/scripts"

set -e


destinationDirRoot=$1
fullSequence=$2
N=$3
temperature=$4
time=$5
suffix=$6 # job label when running on cluster

main=hairpin # Should be in current working dir, or in $PATH
timestep=15
interval=20
wait=0

fps=10
#render="-r -R 3 -f $fps"

echo "Making temporary directory"
outputBaseDir=`mktemp -d`
#outputBaseDir=`mktemp -d /home/other/roald/tmp/cluster.XXXXXXXX`
echo "Made temporary directory"
outputFile="$outputBaseDir/outputFrom_${suffix}"
destinationDir="$destinationDirRoot/$fullSequence/T${temperature}/dt${timestep}_time${time}/N${N}"
destinationFile="$destinationDir/endToEnd_${suffix}"

echo "Making sequence"
seq=`$scriptsdir/hairpinFromStem.py $fullSequence $N`
echo "Made sequence"

echo "Starting main program"
echo ""
$main -t $timestep -T $temperature -I $interval -P $time -W $wait -s $seq -D $outputFile -e -n
echo ""
echo "End of main program"

echo "Making destination directory"
mkdir -p $destinationDir
echo "Made destination directory"

echo "Moving output files"
mv "${outputFile}_endToEnd" "${destinationFile}"

echo "Removing temporary directory"
rm -rf $outputBaseDir

echo "All done! Have an awesome day!"

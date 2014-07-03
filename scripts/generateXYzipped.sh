#!/bin/sh

# Start a simulation with a hairpin strand with XY 'nucleus', wait for it 
# to zip and then dump the resulting world configuration to a file in the 
# given directiory

echo "Starting generateXYzipped.sh"

set -ex

scriptsdir="$HOME/DNA/scripts"

destinationDirRoot=$1
fullSequence=$2
N=$3
temperature=$4
minTime=$5
gamma=$6
suffix=$7 # job label when running on cluster


main=hairpin # Should be in current working dir, or in $PATH
timestep=15
interval=10

NXY=4
NtotStem=`python -c "print($N + $NXY)"` # use this for relaxation time etc


outputBaseDir=`mktemp -d`
outputFile="$outputBaseDir/outputFrom_${suffix}"
destinationDir="$destinationDirRoot/$fullSequence/temp${temperature}_gamma${gamma}_minTime${minTime}_NXY${NXY}/dt${timestep}/N${N}"
destinationFile="$destinationDir/configuration_${suffix}"

echo "Temporary directory is $outputBaseDir"
cd "$outputBaseDir"

echo "Making sequence"

seq=`$scriptsdir/hairpinFromStemXY.py $fullSequence $N $NXY`

echo "Sequence is: $seq"
echo "Starting main program! Only enabling base pairing between XY pairs!"

$main -t $timestep -T $temperature -I $interval -K $minTime -s $seq -w $outputFile -X f -H $NXY:1 -M 10000000 -O $temperature -Q $temperature -U 0 -g $gamma -Y # -D "/dev/null"

echo "Making destination directory: $destinationDir"
mkdir -p $destinationDir

echo "Moving output files!"
mv ${outputFile} ${destinationFile}

echo "Deleting temporary directory"
rm -rf $outputBaseDir

echo "All done! *High five*!"


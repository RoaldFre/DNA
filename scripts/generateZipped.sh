#!/bin/sh

# Start a simulation with a hairpin strand, wait for it to zip and then 
# dump the resulting world configuration to a file in the given directory

echo "Starting generateXYzipped.sh"

set -ex

scriptsdir="$HOME/DNA/scripts"

destinationDirRoot=$1
fullSequence=$2
N=$3
allowedUnboundedFraction=$4
temperature=$5
gamma=$6
suffix=$7 # job label when running on cluster


main=hairpin # Should be in current working dir, or in $PATH
timestep=15
interval=10
zippedRelaxationTime=10

requiredBounded=`python -c "print(int(round($N * (1 - $allowedUnboundedFraction))))"`

outputBaseDir=`mktemp -d`
outputFile="$outputBaseDir/outputFrom_${suffix}"
destinationDir="$destinationDirRoot/$fullSequence/T${temperature}/N${N}"
destinationFile="$destinationDir/configuration_${suffix}"

echo "Temporary directory is $outputBaseDir"
cd "$outputBaseDir"

echo "Making sequence"

seq=`$scriptsdir/hairpinFromStem.py $fullSequence $N`

echo "Sequence is: $seq"
echo "Starting main program!"

$main -t $timestep -T $temperature -I $interval -s $seq -w $outputFile -X f -H $requiredBounded:1 -M 10000000 -O $temperature -Q $temperature -U $zippedRelaxationTime -g $gamma # -D "/dev/null"

echo "Making destination directory: $destinationDir"
mkdir -p $destinationDir

echo "Moving output files!"
mv "${outputFile}" "${destinationFile}"

echo "Deleting temporary directory"
rm -rf $outputBaseDir

echo "All done! *High five*!"



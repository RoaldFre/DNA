#!/bin/sh

# Load file that has fully zipped hairpin, relax that with langevin at low 
# temp to get independent state, then quench to high temperature and record 
# the base pairing.

echo "Starting hairpinUnzippingState.sh wrapper"

set -ex

scriptsdir="$HOME/DNA/scripts"

destinationDirRoot=$1
N=$2
relaxTemp=$3
unzipTemp=$4
samplingTime=$5
gamma=$6
suffix=$7 # job label when running on cluster


main=hairpin # Should be in $PATH
timestep=15
interval=5

stateGlob="/home/other/roald/clusterdata/zipped/CACTCAGAGAGTGACTGACTCTCAGACTCACACAGAGAGTCACTGTCTGACTCTCTCTGAGACACTGAGAGTGAGAGTGACTCTGAGTGAGTCACAGTGA/T${relaxTemp}/N${N}/configuration_*";

echo "Getting wait time"
factor=`python -c "print($gamma / 5e12 * 50)"` # Was originally *10!
wait=`$scriptsdir/hairpinRelaxationTimeVariable.py $N $factor`
echo "Wait time: $wait"


outputBaseDir=`mktemp -d`
outputFile="$outputBaseDir/state_${suffix}"
destinationDir="$destinationDirRoot/relaxT${relaxTemp}_unzipT${unzipTemp}_time${samplingTime}_relaxFact${factor}/dt${timestep}/g${gamma}/N${N}"

relaxFile="${outputFile}_finalRelax"
endFile="${outputFile}_endConfig"

echo "Temporary directory is $outputBaseDir"

initialWorldFile=`find $stateGlob | head -n 1`
echo "Loading world configuration from $initialWorldFile"

relaxationDir="$outputBaseDir/relaxation"
mkdir -p $relaxationDir
cd $relaxationDir
echo
echo "Starting langevin relaxation of $wait nanoseconds, writing to $relaxFile"
time $main -i l -d "$initialWorldFile" -w "$relaxFile" -t $timestep -T $relaxTemp -g $gamma -I $interval -P $wait

echo
echo "Starting actual measurement!"
time $main -d "$relaxFile" -w "$endFile" -t $timestep -T $unzipTemp -I $interval -P $samplingTime -D "$outputFile" -a $unzipTemp -X s -g $gamma -e

echo
echo "Converting output to native Octave format"
octave -p "$scriptsdir" --eval "toNative('${outputFile}')"

echo "Making destination directory: $destinationDir"
mkdir -p "$destinationDir"

echo "Removing relaxation files"
rm -r $relaxationDir

echo "Moving output files!"
mv "${outputBaseDir}"/* "${destinationDir}"

echo "Deleting temporary directory"
rm -rf "$outputBaseDir"

echo "All done! *High five*!"


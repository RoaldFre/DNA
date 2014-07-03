#!/bin/sh

# Load file that already has XY zipped, relax with langevin at low temp and 
# only XY, then enable full pairing and watch how it zips at (the same) low 
# temperature. All the time, a constant end-to-end force is added to 
# emulate arms of optical trap.

echo "Starting hairpinStateWithRelaxationXYlowTempGammaForce.sh wrapper"

set -ex

scriptsdir="$HOME/DNA/scripts"

destinationDirRoot=$1
N=$2
temperature=$3
samplingTime=$4
gamma=$5
force=$6
suffix=$7 # job label when running on cluster


main=hairpin # Should be in current working dir, or in $PATH
timestep=15
interval=10 # was 100 for gamma = 5e12 (maybe base on samplingTime to get fixed number of samples?)

statedir="/home/other/roald/clusterdata/XYzipped_force/CACTCAGAGAGTGACTGACTCTCAGACTCACACAGAGAGTCACTGTCTGACTCTCTCTGAGACACTGAGAGTGAGAGTGACTCTGAGTGAGTCACAGTGA/temp${temperature}_gamma*_minTime*_NXY4/dt15/F${force}"

NXY=4
NtotStem=`python -c "print($N + $NXY)"` # use this for relaxation time etc

echo "Getting wait time"
factor=`python -c "print($gamma / 5e12 * 10)"`
wait=`$scriptsdir/hairpinRelaxationTimeVariable.py $NtotStem $factor`
echo "Wait time: $wait"


outputBaseDir=`mktemp -d`
outputFile="$outputBaseDir/state_${suffix}"
destinationDir="$destinationDirRoot/T${temperature}_time${samplingTime}_NXY${NXY}_relaxFact${factor}/dt${timestep}/g${gamma}/F${force}/N${N}"

relaxFile="${outputFile}_finalRelax"
endFile="${outputFile}_endConfig"

echo "Temporary directory is $outputBaseDir"

initialWorldFile=`find $statedir/N${N}/* | head -n 1`
echo "Loading world configuration from $initialWorldFile"

relaxationDir="$outputBaseDir/relaxation"
mkdir -p $relaxationDir
echo
echo "Starting langevin relaxation of at least $wait nanoseconds, only XY pairs, writing to $relaxFile, requiring $NXY bound base pairs at end"
time $main -i l -d "$initialWorldFile" -w "$relaxFile" -t $timestep -T $temperature -I $interval -K $wait -X f -H $NXY:1 -M 10000000 -O $temperature -Q $temperature -U 0 -D "$relaxationDir/relaxation" -Y -g $gamma -E $force:0:0 -e

echo
echo "Starting actual measurement!"
time $main -d "$relaxFile" -w "$endFile" -t $timestep -T $temperature -I $interval -P $samplingTime -D "$outputFile" -a $temperature -X s -g $gamma -E $force:0:0 -e

echo
echo "Converting output to native Octave format"
octave -p "$scriptsdir" --eval "toNative('${outputFile}')"

echo "Making destination directory: $destinationDir"
mkdir -p "$destinationDir"

# Keep the end-to-end record of the relaxation to check stability of 
# end-to-end vector
mv "$relaxationDir/relaxation_endToEnd" "${outputFile}_relaxEndToEnd"
echo "Removing relaxation files"
rm -r $relaxationDir

echo "Moving output files!"
mv "${outputBaseDir}"/* "${destinationDir}"

echo "Deleting temporary directory"
rm -rf "$outputBaseDir"

echo "All done! *High five*!"

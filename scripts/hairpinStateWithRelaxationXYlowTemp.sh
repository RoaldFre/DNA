#!/bin/sh

# Load file that already has XY zipped, relax with monte carlo at low 
# temperature and only XY pairing enabled, do a short final relaxation with 
# langevin at low temp and only XY, then enable XY pairing and watch how it 
# zips at (the same) low temperature.

echo "Starting hairpinStateWithRelaxationXYlowTemp.sh wrapper"

set -ex

scriptsdir="$HOME/DNA/scripts"
statedir="/home/other/roald/clusterdata/XYzipped/CACTCAGAGAGTGACTGACTCTCAGACTCACACAGAGAGTCACTGTCTGACTCTCTCTGAGACACTGAGAGTGAGAGTGACTCTGAGTGAGTCACAGTGA/temp90C_minTime*_NXY4/dt15/"

destinationDirRoot=$1
N=$2
temperature=$3
samplingTime=$4
monteCarloSweeps=$5
wait=$6
suffix=$7 # job label when running on cluster


main=hairpin # Should be in current working dir, or in $PATH
timestep=15
interval=100

NXY=4
NtotStem=`python -c "print($N + $NXY)"` # use this for relaxation time etc

outputBaseDir=`mktemp -d`
outputFile="$outputBaseDir/state_${suffix}"
destinationDir="$destinationDirRoot/T${temperature}_time${samplingTime}_sweeps${monteCarloSweeps}_wait${wait}_NXY${NXY}/dt${timestep}/N${N}"

MCfile="${outputFile}_MCrelax"
finalRelaxFile="${outputFile}_finalRelax"
endFile="${outputFile}_endConfig"

echo "Temporary directory is $outputBaseDir"

initialWorldFile=`find $statedir/N${N}/* | head -n 1`
echo "Loading world configuration from $initialWorldFile"

echo
echo "Starting initial monte carlo relaxation of $monteCarloSweeps sweeps, only XY pairs, writing to $MCfile"
time $main -T $temperature -i m -m $monteCarloSweeps -d "$initialWorldFile" -w "$MCfile" -Y

relaxationDir="$outputBaseDir/relaxation"
mkdir -p $relaxationDir
echo
echo "Starting final langevin relaxation of at least $wait nanoseconds, only XY pairs, writing to $finalRelaxFile, requiring $NXY bound base pairs at end"
time $main -i l -d "$MCfile" -w "$finalRelaxFile" -t $timestep -T $temperature -I $interval -K $wait -X f -H $NXY:1 -M 10000000 -O $temperature -Q $temperature -U 0 -D "$relaxationDir/relaxation" -Y

echo
echo "Starting actual measurement!"
time $main -d "$finalRelaxFile" -w "$endFile" -t $timestep -T $temperature -I $interval -P $samplingTime -D "$outputFile" -a $temperature -X s -e

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

#!/bin/sh

# Load file that already has XY zipped, relax with monte carlo at low 
# temperature and only XY pairing enabled, do a short final relaxation with 
# langevin at low temp and only XY, then enable XY pairing and watch how it 
# zips at (the same) low temperature.

echo "Starting hairpinStateWithRelaxationXYlowTempGamma.sh wrapper"

set -ex

scriptsdir="$HOME/DNA/scripts"

destinationDirRoot=$1
N=$2
temperature=$3
samplingTime=$4
monteCarloSweeps=$5
wait=$6
gamma=$7
timestep=$8
suffix=$9 # job label when running on cluster


main=hairpin # Should be in current working dir, or in $PATH
interval=5

statedir="/home/other/roald/clusterdata/XYzipped_bpRepulsive/CACTCAGAGAGTGACTGACTCTCAGACTCACACAGAGAGTCACTGTCTGACTCTCTCTGAGACACTGAGAGTGAGAGTGACTCTGAGTGAGTCACAGTGA/temp${temperature}_gamma*_minTime*_NXY4/dt15/"

NXY=4
NtotStem=`python -c "print($N + $NXY)"` # use this for relaxation time etc


outputBaseDir=`mktemp -d`
outputFile="$outputBaseDir/state_${suffix}"
destinationDir="$destinationDirRoot/T${temperature}/g${gamma}/N${N}/dt${timestep}/time${samplingTime}_sweeps${monteCarloSweeps}_wait${wait}_NXY${NXY}"

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
time $main -i l -d "$MCfile" -w "$finalRelaxFile" -t $timestep -T $temperature -I $interval -K $wait -X f -H $NXY:1 -M 10000000 -O $temperature -Q $temperature -U 0 -D "$relaxationDir/relaxation" -Y -g $gamma

echo
echo "Starting actual measurement!"
time $main -d "$finalRelaxFile" -w "$endFile" -t $timestep -T $temperature -I $interval -P $samplingTime -D "$outputFile" -a $temperature -X s -g $gamma -e -p -q -J

echo
echo "Converting output to native Octave format"
$scriptsdir/toNative.sh ${outputFile} s
$scriptsdir/toNative.sh ${outputFile}_basePairingEnergy s
$scriptsdir/toNative.sh ${outputFile}_basePairingDistance s
$scriptsdir/toNative.sh ${outputFile}_gyrationRadius s
$scriptsdir/toNative.sh ${outputFile}_endToEnd s
$scriptsdir/toNative.sh ${outputFile}_forceVelFric s
$scriptsdir/toNative.sh ${outputFile}_forceVelFricP s
$scriptsdir/toNative.sh ${outputFile}_forceVelP s
$scriptsdir/toNative.sh ${outputFile}_forceVelS s

echo "Making destination directory: $destinationDir"
mkdir -p "$destinationDir"

echo "Removing relaxation files"
rm -r $relaxationDir

echo "Moving output files!"
mv "${outputBaseDir}"/* "${destinationDir}"

echo "Deleting temporary directory"
rm -rf "$outputBaseDir"

echo "All done! *High five*!"

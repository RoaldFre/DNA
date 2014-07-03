#!/bin/sh

echo "Starting hairpinStateWithRelaxationXY2.sh wrapper"

set -e

scriptsdir="$HOME/DNA/scripts"
statedir="/home/other/roald/clusterdata/XYzipped/CACTCAGAGAGTGACTGACTCTCAGACTCACACAGAGAGTCACTGTCTGACTCTCTCTGAGACACTGAGAGTGAGAGTGACTCTGAGTGAGTCACAGTGA/temp90C_minTime*_NXY4/dt15/"

destinationDirRoot=$1
N=$2
relaxTemp=$3
samplingTemp=$4
samplingTime=$5
monteCarloSweeps=$6
suffix=$7 # job label when running on cluster


main=hairpin # Should be in current working dir, or in $PATH
timestep=15
interval=100

NXY=4
NtotStem=`python -c "print($N + $NXY)"` # use this for relaxation time etc

echo "Getting wait time"
wait=`$scriptsdir/hairpinRelaxationTimeShort.py $NtotStem`
echo "Wait time: $wait"


outputBaseDir=`mktemp -d`
outputFile="$outputBaseDir/state_${suffix}"
destinationDir="$destinationDirRoot/relaxT${relaxTemp}_sampleT${samplingTemp}_time${samplingTime}_sweeps${monteCarloSweeps}_NXY${NXY}/dt${timestep}/N${N}"

MCfile="${outputFile}_MCrelax"
finalRelaxFile="${outputFile}_finalRelax"
endFile="${outputFile}_endConfig"

echo "Temporary directory is $outputBaseDir"

initialWorldFile=`find $statedir/N${N}/* | head -n 1`
echo "Loading world configuration from $initialWorldFile"

echo
echo "Starting initial monte carlo relaxation of $monteCarloSweeps sweeps, writing to $MCfile"
echo $main -T $relaxTemp -i m -m $monteCarloSweeps -d "$initialWorldFile" -w "$MCfile"
$main -T $relaxTemp -i m -m $monteCarloSweeps -d "$initialWorldFile" -w "$MCfile"

echo
echo "Starting final langevin relaxation of at least $wait nanoseconds, writing to $finalRelaxFile, requiring $NXY bound base pairs at end"
echo $main -i l -d "$MCfile" -w "$finalRelaxFile" -t $timestep -T $relaxTemp -I $interval -K $wait -X f -H $NXY -M 10000000 -O $relaxTemp -Q $relaxTemp -U 0 -D "/dev/null"
$main -i l -d "$MCfile" -w "$finalRelaxFile" -t $timestep -T $relaxTemp -I $interval -K $wait -X f -H $NXY -M 10000000 -O $relaxTemp -Q $relaxTemp -U 0 -D "/dev/null"

echo
echo "Starting actual measurement!"
$main -d "$finalRelaxFile" -w "$endFile" -t $timestep -T $samplingTemp -I $interval -P $samplingTime -D "$outputFile" -a $samplingTemp -X s

echo
echo "Converting output to native Octave format"
octave -p "$scriptsdir" --eval "toNative('${outputFile}')"

echo "Making destination directory: $destinationDir"
mkdir -p "$destinationDir"

echo "Moving output files!"
mv "${outputBaseDir}"/* "${destinationDir}"

echo "Deleting temporary directory"
rm -rf "$outputBaseDir"

echo "All done! *High five*!"

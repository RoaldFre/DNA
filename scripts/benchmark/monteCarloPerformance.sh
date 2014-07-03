#!/bin/bash

set -e

nRuns=10

cmd=./generateMonteCarloPerformance.sh

#outputFile=monteCarloPerformance

rm -f $outputFile

$cmd 320 12 $nRuns   # >> $outputFile
$cmd 269 25 $nRuns   # >> $outputFile
$cmd 226 25 $nRuns   # >> $outputFile
$cmd 190 50 $nRuns   # >> $outputFile
$cmd 160 50 $nRuns   # >> $outputFile
$cmd 135 100 $nRuns  # >> $outputFile
$cmd 113 100 $nRuns  # >> $outputFile
$cmd 95 200 $nRuns   # >> $outputFile
$cmd 80 200 $nRuns   # >> $outputFile
$cmd 67 400 $nRuns   # >> $outputFile
$cmd 57 400 $nRuns   # >> $outputFile
$cmd 48 1000 $nRuns  # >> $outputFile
$cmd 40 1000 $nRuns  # >> $outputFile
$cmd 34 2000 $nRuns  # >> $outputFile
$cmd 28 2000 $nRuns  # >> $outputFile
$cmd 24 4000 $nRuns  # >> $outputFile
$cmd 20 4000 $nRuns  # >> $outputFile
$cmd 17 8000 $nRuns  # >> $outputFile
$cmd 14 8000 $nRuns  # >> $outputFile
$cmd 12 16000 $nRuns # >> $outputFile
$cmd 10 16000 $nRuns # >> $outputFile
$cmd 8 16000 $nRuns # >> $outputFile
$cmd 7 32000 $nRuns # >> $outputFile
$cmd 6 32000 $nRuns # >> $outputFile
$cmd 5 16000 $nRuns # >> $outputFile

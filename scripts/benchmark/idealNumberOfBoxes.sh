#!/bin/bash

set -e

nRuns=1

cmd=./generateIdealNumberOfBoxes.sh

N=320; $cmd $N 0.006 27 $nRuns > ideal_$N &
N=269; $cmd $N 0.01  26 $nRuns > ideal_$N &
N=226; $cmd $N 0.02  25 $nRuns > ideal_$N &
N=190; $cmd $N 0.04  24 $nRuns > ideal_$N &
N=160; $cmd $N 0.05  23 $nRuns > ideal_$N &
N=135; $cmd $N 0.05  22 $nRuns > ideal_$N &
N=113; $cmd $N 0.1  21 $nRuns > ideal_$N &
N=95 ; $cmd $N 0.2  20 $nRuns > ideal_$N &
N=80 ; $cmd $N 0.2  18 $nRuns > ideal_$N &
N=67 ; $cmd $N 0.2  16 $nRuns > ideal_$N &
N=57 ; $cmd $N 0.2  14 $nRuns > ideal_$N &
N=48 ; $cmd $N 0.2  12 $nRuns > ideal_$N &
N=40 ; $cmd $N 0.4  10 $nRuns > ideal_$N &
N=34 ; $cmd $N 0.4  8  $nRuns > ideal_$N &
N=28 ; $cmd $N 0.4  7  $nRuns > ideal_$N &
N=24 ; $cmd $N 0.4  6  $nRuns > ideal_$N &
N=20 ; $cmd $N 0.4  5  $nRuns > ideal_$N &
N=17 ; $cmd $N 0.8  5  $nRuns > ideal_$N &
N=14 ; $cmd $N 0.8  5  $nRuns > ideal_$N &
N=12 ; $cmd $N 0.8  5  $nRuns > ideal_$N &
N=10 ; $cmd $N 0.8  5  $nRuns > ideal_$N &

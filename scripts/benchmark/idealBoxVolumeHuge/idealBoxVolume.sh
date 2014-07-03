#!/bin/bash

set -e

nRuns=6

cmd=./generateIdealBoxVolume.sh

S=1000
Bmin=10
Bstep=10
Bmax=50
sweeps=50

N=1280; $cmd $N 0.005 $sweeps $S $Bmin $Bmax $nRuns > ideal_$N
N=640; $cmd $N 0.01 $sweeps $S $Bmin $Bmax $nRuns > ideal_$N


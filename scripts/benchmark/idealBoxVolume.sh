#!/bin/bash

set -e

nRuns=5

cmd=./generateIdealBoxVolume.sh

S=2000
Bmin=1
Bmax=100
sweeps=100

N=320; $cmd $N 0.2 $sweeps $S $Bmin $Bmax $nRuns > ideal_$N
N=160; $cmd $N 0.4 $sweeps $S $Bmin $Bmax $nRuns > ideal_$N
N=80;  $cmd $N 0.8 $sweeps $S $Bmin $Bmax $nRuns > ideal_$N
N=40;  $cmd $N   2 $sweeps $S $Bmin $Bmax $nRuns > ideal_$N
N=20;  $cmd $N   4 $sweeps $S $Bmin $Bmax $nRuns > ideal_$N
N=10;  $cmd $N   8 $sweeps $S $Bmin $Bmax $nRuns > ideal_$N


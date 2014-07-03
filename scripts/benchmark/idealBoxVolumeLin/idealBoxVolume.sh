#!/bin/bash

set -e

nRuns=10

cmd=./generateIdealBoxVolume.sh

S=2000
Bmin=1
Bstep=2
Bmax=99
sweeps=80

N=400; $cmd $N 0.1  $sweeps $S $Bmin $Bstep $Bmax $nRuns > ideal_$N
N=350; $cmd $N 0.11 $sweeps $S $Bmin $Bstep $Bmax $nRuns > ideal_$N
N=300; $cmd $N 0.13 $sweeps $S $Bmin $Bstep $Bmax $nRuns > ideal_$N
N=250; $cmd $N 0.16 $sweeps $S $Bmin $Bstep $Bmax $nRuns > ideal_$N
N=200; $cmd $N 0.2  $sweeps $S $Bmin $Bstep $Bmax $nRuns > ideal_$N
N=150; $cmd $N 0.27 $sweeps $S $Bmin $Bstep $Bmax $nRuns > ideal_$N
N=100; $cmd $N 0.4  $sweeps $S $Bmin $Bstep $Bmax $nRuns > ideal_$N
N=50;  $cmd $N 0.8  $sweeps $S $Bmin $Bstep $Bmax $nRuns > ideal_$N




#N=160; $cmd $N 0.4 $sweeps $S $Bmin $Bstep $Bmax $nRuns > ideal_$N
#N=80;  $cmd $N 0.8 $sweeps $S $Bmin $Bstep $Bmax $nRuns > ideal_$N
#N=40;  $cmd $N   2 $sweeps $S $Bmin $Bstep $Bmax $nRuns > ideal_$N
#N=20;  $cmd $N   4 $sweeps $S $Bmin $Bstep $Bmax $nRuns > ideal_$N
#N=10;  $cmd $N   8 $sweeps $S $Bmin $Bstep $Bmax $nRuns > ideal_$N


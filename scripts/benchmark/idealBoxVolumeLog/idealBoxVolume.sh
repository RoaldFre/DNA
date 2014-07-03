#!/bin/bash

set -e

cmd=./generateIdealBoxVolume.sh

S=2000
Bmin=1
Bstep=2
Bmax=99
Bfinal=100
sweeps=100

#N=320; $cmd $N 0.2 $sweeps $S $Bmin $Bstep $Bmax $Bfinal 15 > ideal_$N
#N=226; $cmd $N 0.3 $sweeps $S $Bmin $Bstep $Bmax $Bfinal 14 > ideal_$N
#N=160; $cmd $N 0.4 $sweeps $S $Bmin $Bstep $Bmax $Bfinal 13 > ideal_$N
#N=113; $cmd $N 0.6 $sweeps $S $Bmin $Bstep $Bmax $Bfinal 12 > ideal_$N
#N=80;  $cmd $N 0.8 $sweeps $S $Bmin $Bstep $Bmax $Bfinal 11 > ideal_$N
#N=57;  $cmd $N 1.4 $sweeps $S $Bmin $Bstep $Bmax $Bfinal 10 > ideal_$N
N=40;  $cmd $N   2 $sweeps $S $Bmin $Bstep $Bmax $Bfinal 9 > ideal_$N
N=28;  $cmd $N   3 $sweeps $S $Bmin $Bstep $Bmax $Bfinal 8 > ideal_$N
N=20;  $cmd $N   4 $sweeps $S $Bmin $Bstep $Bmax $Bfinal 7 > ideal_$N
N=14;  $cmd $N   6 $sweeps $S $Bmin $Bstep $Bmax $Bfinal 6 > ideal_$N
N=10;  $cmd $N   8 $sweeps $S $Bmin $Bstep $Bmax $Bfinal 5 > ideal_$N


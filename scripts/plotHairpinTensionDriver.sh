#!/bin/bash

scriptsdir="$HOME/DNA/scripts"

nproc=1

xargs -P $nproc -n 1 octave-3.6.4 -q -p $scriptsdir --eval << EOF
"pkg load optim; plotHairpinFvPLowTempGammaDriver2"
"pkg load optim; plotHairpinFvSLowTempGammaDriver2" 
"pkg load optim; plotHairpinFvfLowTempGammaDriver2" 
"pkg load optim; plotHairpinFvfPLowTempGammaDriver2" 
EOF







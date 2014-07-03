#!/bin/bash

scriptsdir="$HOME/DNA/scripts"

nproc=1

xargs -P $nproc -n 1 octave-3.6.4 -q -p $scriptsdir --eval << EOF
"pkg load optim; makeHairpinStatePlots"
"pkg load optim; makeHairpinEnergyPlots"
"pkg load optim; makeHairpinEndToEndPlots"
"pkg load optim; makeHairpinGyradPlots"
EOF


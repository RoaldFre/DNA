#!/bin/bash

scriptsdir="$HOME/DNA/scripts"

nproc=20

xargs -P $nproc -n 1 octave -q -p $scriptsdir --eval << EOF
"plotEndToEndZippingDriver"
"plotHairpinEnergyLowTempGammaDriver"
"plotHairpinGyradLowTempGammaDriver"
"plotHairpinStateLowTempGammaDriver"
EOF
#"plotHairpinFvfLowTempGammaDriver"


#octave -q -p $scriptsdir --eval "plotHairpinFvfLowTempGammaDriver" $suffix
###octave -q -p $scriptsdir --eval "plotHairpinFvfLowTempGammaDriver" &
#echo $!

#wait

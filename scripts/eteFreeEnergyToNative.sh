#!/bin/bash

file=$1

scriptsdir="$HOME/DNA/scripts"

octave -q -p $scriptsdir --eval "eteFreeEnergyToNative('$file')"



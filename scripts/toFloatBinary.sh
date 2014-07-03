#!/bin/bash

set -e

file=$1

scriptsdir="$HOME/DNA/scripts"

octave -q -p $scriptsdir --eval "toFloatBinary('$file');"

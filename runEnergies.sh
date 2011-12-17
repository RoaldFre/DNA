#!/bin/sh

set -e

numMonomers=50
samples=300
measInterval=100
timestep=1
thermostatCoupling=0

./main $numMonomers -E $measInterval -s $samples -t $timestep -c $thermostatCoupling -r

octave --persist -q plotEnergies.m

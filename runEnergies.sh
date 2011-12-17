#!/bin/sh

set -e

numMonomers=30
samples=1000
timestep=3
measInterval=200
thermostatCoupling=2000

./main $numMonomers -E $measInterval -s $samples -t $timestep -c $thermostatCoupling -r

octave --persist -q plotEnergies.m

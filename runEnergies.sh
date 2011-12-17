#!/bin/sh

set -e

numMonomers=20
samples=1000
timestep=2
measInterval=300
thermostatCoupling=10000
measureWait=200000

./main $numMonomers -E $measInterval -s $samples -t $timestep -c $thermostatCoupling -w $measureWait -r

octave --persist -q plotEnergies.m

#!/bin/sh

set -e

numMonomers=20
timestep=1
measInterval=300
temperature=300

info="\$N\$ = $numMonomers, \$\Delta t\$ = $timestep\,fs"


samples=1000
thermostatCoupling=0
measureWait=0
./main $numMonomers -E $measInterval -s $samples -t $timestep -c $thermostatCoupling -T $temperature -w $measureWait -r
filename="fixedEnergy"
caption="Fixed energy. Initial configuration in flat vertical column, spacing according to equilibrium distances. This corresponds to an initial temperature of approximately 9000\,K(!). Parameters: $info"
octave --eval "plotEnergies('$filename', '$caption')"


samples=250
thermostatCoupling=10000
measureWait=0
./main $numMonomers -E $measInterval -s $samples -t $timestep -c $thermostatCoupling -T $temperature -w $measureWait -r
filename="fixedTemperature"
caption="Fixed temperature with Berendsen thermostat, initial relaxation phase. Parameters: \$T\$ = $temperature\,K, \$\\\\tau\$ = $thermostatCoupling\,fs, $info."
octave -q --eval "plotEnergies('$filename', '$caption')"


samples=1000
thermostatCoupling=10000
measureWait=200000
./main $numMonomers -E $measInterval -s $samples -t $timestep -c $thermostatCoupling -T $temperature -w $measureWait -r
filename="fixedTemperatureRelaxed"
caption="Fixed temperature with Berendsen thermostat, after relaxation. Parameters: \$T\$ = $temperature\,K, \$\\\\tau\$ = $thermostatCoupling\,fs, $info."
octave -q --eval "plotEnergies('$filename', '$caption')"


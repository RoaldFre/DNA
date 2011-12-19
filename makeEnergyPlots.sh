#!/bin/sh

set -e


for numMonomers in 4 10 20 50 100
do

	timestep=1
	measInterval=300
	temperature=300

	info="\$N\$ = $numMonomers, \$\Delta t\$ = $timestep\,fs"


	samples=1000
	thermostatCoupling=0
	measureWait=0
	./main $numMonomers -E $measInterval -s $samples -t $timestep -c $thermostatCoupling -T $temperature -w $measureWait -r
	filename="fixedEnergy-N$numMonomers"
	caption="Fixed energy. Initial configuration in flat vertical column, spacing according to equilibrium distances. This corresponds to an initial temperature of approximately 9000\,K(!). Parameters: $info."
	octave --eval "plotEnergies('$filename', '$caption')"


	samples=250
	thermostatCoupling=10000
	measureWait=0
	./main $numMonomers -E $measInterval -s $samples -t $timestep -c $thermostatCoupling -T $temperature -w $measureWait -r
	filename="fixedTemperature-N$numMonomers"
	caption="Fixed temperature with Berendsen thermostat, initial relaxation phase. Parameters: \$T\$ = $temperature\,K, \$\\tau\$ = $thermostatCoupling\,fs, $info."
	octave -q --eval "plotEnergies('$filename', '$caption')"


	samples=1000
	thermostatCoupling=10000
	measureWait=200000
	./main $numMonomers -E $measInterval -s $samples -t $timestep -c $thermostatCoupling -T $temperature -w $measureWait -r
	filename="fixedTemperatureRelaxed-N$numMonomers"
	caption="Fixed temperature with Berendsen thermostat, after relaxation. Parameters: \$T\$ = $temperature\,K, \$\\tau\$ = $thermostatCoupling\,fs, $info."
	octave -q --eval "plotEnergies('$filename', '$caption')"

done

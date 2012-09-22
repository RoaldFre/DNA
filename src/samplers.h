#ifndef _SAMPLERS_H_
#define _SAMPLERS_H_

#include "measure.h"

/* Sampler that averages the temperature (in a naive way, so don't use it 
 * on a huge amount of samples [hundreds of thousands] or you will loose 
 * precision). Prints the result at the end. */
Sampler averageTemperatureSampler(void);
/* Sampler that dumps the kinetic temperature. */
Sampler temperatureSampler(void);

/* Sampler that dumps physics stats. */
Sampler dumpStatsSampler(void);

/* Sampler that dumps the position (in m) of the particle. */
Sampler particlePositionSampler(Particle *p);
/* Sampler that dumps the Center Of Mass position (in m) of the strand. */
Sampler strandCOMSampler(Strand *s);

/* Sampler that dumps the squared displacement (in m^2) of the particle. */
Sampler particleSquaredDisplacementSampler(Particle *p);
/* Sampler that dumps the squared displacement of the Center Of Mass (in 
 * m^2) as well as the full position of the COM of the strand. */
Sampler strandCOMSquaredDisplacementSampler(Strand *s);

/* Sampler that dumps the end-to-end distance of the given strand. The 
 * end-to-end distance is the distance between the Center Of Mass of the 
 * first monomer and the COM of the last monomer. */
Sampler endToEndDistSampler(Strand *strand);


typedef struct {
	/* A pair is 'bound' if its base pair potential is lower than the 
	 * given threshold. */
	real energyThreshold;
} BasePairingConfig;
/* Sampler that counts the number of base pair bindings. */
Sampler basePairingSampler(BasePairingConfig *bpc);

typedef struct {
	/* Temperature at which to perform the time-to-zipping measurement */
	real zippingTemperature;

	/* Temperature at which to perform the time-to-unzipping measurement */
	real unzippingTemperature;

	/* A pair is 'bound' if its base pair potential is lower than the 
	 * given threshold. */
	real energyThreshold;

	/* Time the molecule has to be fully zipped(unzipped) in order for 
	 * it to be considered a stable(melted) hairpin */
	real confirmationTime;

	/* Minimum number of bound basepairs required when determining 
	 * whether a hairpin is fully zipped */
	int requiredBoundBPs; 

	/* Maximum number of bound basepairs to allow when determining 
	 * whether a hairpin is fully unzipped */
	int allowedBoundBPs; 

	/* Time to wait after having confirmed zipping before switching the 
	 * temperature to 'unzippingTemperature' and waiting for the 
	 * hairpin to unzip. If after this time, the hairin is not fully 
	 * zipped, we will wait untill it is and only then start sampling. */
	real zippedRelaxationTime;
} HairpinFormationSamplerConfig;

/* When this sampler starts, it sets the temperature to the given 
 * zippingTemperature in the config file and then waits until the hairpin 
 * in the world gets formed. It then keeps the strand in this hairpin state 
 * for the given relaxation time, after which the temperature gets set to 
 * the unzippingTemperature. Then, the time till unzipping again is 
 * measured.
 * See the documentation of HairpinFormationSamplerConfig for more info. */
Sampler hairpinFormationSampler(HairpinFormationSamplerConfig *hfc);


typedef struct {
	/* A pair is 'bound' if its base pair potential is lower than the 
	 * given threshold. */
	real energyThreshold;

	/* Time to wait for relaxation after bumping temperature */
	real relaxationTime;
	/* Time to let the system thermally equilibrate before starting to 
	 * measure */
	real measureTime;

	/* Start measuring at 'Tstart'. Then change the temperature by a 
	 * difference 'Tstep' and measure again. Do this until you have 
	 * measured a total of 'numSteps' times. */
	real Tstart;
	real Tstep;
	int numSteps;

	/* If 'verbose' is true: dump full configuration of hairpin at 
	 * every sample tick. If it is false: only dump the average number 
	 * of bound base pairs (also dumped when 'verbose' is enabled). */
	bool verbose;
} HairpinMeltingTempSamplerConfig;

Sampler hairpinMeltingTempSampler(HairpinMeltingTempSamplerConfig *hsc);

/* A trivial sampler that does nothing. Useful for debugging purposes. */
Sampler trivialSampler(void);

#endif


#ifndef _SAMPLERS_H_
#define _SAMPLERS_H_

#include "measure.h"

/* Sampler that averages the temperature (in a naive way, so don't use it 
 * on a huge amount of samples or you will loose precision). Prints the 
 * result at the end. */
Sampler averageTemperatureSampler(void);

/* Sampler that dumps physics stats. */
Sampler dumpStatsSampler(void);

/* Sampler that dumps the position (in fm) of the particle. */
Sampler particlePositionSampler(Particle *p);
/* Sampler that dumps the Center Of Mass position (in fm) of the strand. */
Sampler strandCOMSampler(Strand *s);

/* Sampler that dumps the squared displacement (in fm^2) of the particle. */
Sampler particleSquaredDisplacementSampler(Particle *p);
/* Sampler that dumps the squared displacement of Center Of Mass (in fm^2) 
 * the strand. */
Sampler strandCOMSquaredDisplacementSampler(Strand *s);

typedef struct {
	/* A pair is 'bound' if its base pair potential is lower than the 
	 * given threshold. */
	double energyThreshold;

	/* When starting to sample, set the thermostat to this temperature 
	 * value. */
	double T;
} BasePairingConfig;
/* Sampler that counts the number of base pair bindings. */
Sampler basePairingSampler(BasePairingConfig *bpc);

typedef struct {
	double energyThreshold;
	double confirmationTime; /* Time the molecule has to be fully 
				    zipped in order for it to be considered 
				    a stable hairpin */
	int allowedUnboundBPs; /* Number of unbound basepairs to allow */
} HairpinFormationSamplerConfig;

Sampler hairpinFormationSampler(HairpinFormationSamplerConfig *hfc);


typedef struct {
	/* A pair is 'bound' if its base pair potential is lower than the 
	 * given threshold. */
	double energyThreshold;

	/* Time to wait for relaxation after bumping temperature */
	double relaxationTime;
	/* Time to let the system thermally equilibrate before starting to 
	 * measure */
	double measureTime;

	/* Start measuring at 'Tstart'. Then change the temperature by a 
	 * difference 'Tstep' and measure again. Do this until you have 
	 * measured a total of 'numSteps' times. */
	double Tstart;
	double Tstep;
	int numSteps;

	/* If 'verbose' is true: dump full configuration of hairpin at 
	 * every sample tick. If it is false: only dump the average number 
	 * of bound base pairs (also dumped when 'verbose' is enabled). */
	bool verbose;
} HairpinMeltingTempSamplerConfig;

Sampler hairpinMeltingTempSampler(HairpinMeltingTempSamplerConfig *hsc);

#endif

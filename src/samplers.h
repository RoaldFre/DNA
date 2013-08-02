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

/* Sampler that dumps the gyration radius of the given strand. */
Sampler gyrationRadiusSampler(Strand *strand);


typedef enum
{
	BPSM_BOOLEAN,  /* Dump bound/unbound data, based on energyThreshold */
	BPSM_DISTANCE, /* Dump the distance between corresponding base pairs */
	BPSM_ENERGY,   /* Dump the energy of corresponding base pairs */
	/* TODO< Currently, the boolean mode is the only one that works for 
	 * arbitrary configs. For hairpin-type configurations, all three 
	 * modes work. */
} BasePairingSamplerMode;
typedef struct {
	BasePairingSamplerMode mode;

	/* A pair is 'bound' if its base pair potential is lower than the 
	 * given threshold. */
	double energyThreshold;
} BasePairingConfig;
/* Sampler that counts the number of base pair bindings. */
Sampler basePairingSampler(BasePairingConfig *bpc);

typedef struct {
	/* Temperature at which to perform the time-to-zipping measurement */
	double zippingTemperature;

	/* Temperature at which to perform the time-to-unzipping measurement */
	double unzippingTemperature;

	/* A pair is 'bound' if its base pair potential is lower than the 
	 * given threshold. */
	double energyThreshold;

	/* Minimum number of bound basepairs required when determining 
	 * whether a hairpin is fully zipped */
	int requiredBoundBPs; 

	/* Maximum number of bound basepairs to allow when determining 
	 * whether a hairpin is fully unzipped */
	int allowedBoundBPs; 

	/* Minimum number of bound basepairs required when determining 
	 * whether there has been nucleation for the zipping. Nucleation is 
	 * defined as the first time the number of bound base pairs reaches 
	 * this threshold, and stays above (or equals) this threshold all 
	 * the way up to zipping confirmation. */
	int nucleationBoundBPs; 

	/* Time the molecule has to be fully zipped (according to the 
	 * requiredBoundBPs criterium) in order for it to be considered a 
	 * stable hairpin */
	double zipConfirmationTime;

	/* Time the molecule has to be fully unzipped (according to the 
	 * allowedBoundBPs criterium) in order for it to be considered a 
	 * stable melted hairpin */
	double unzipConfirmationTime;

	/* Time to wait after having confirmed zipping before switching the 
	 * temperature to 'unzippingTemperature' and waiting for the 
	 * hairpin to unzip. If after this time, the hairin is not fully 
	 * zipped, we will wait untill it is and only then start sampling. */
	double zippedRelaxationTime;

	/* In case the zipping happenend very fast, extend the 
	 * zippedRelaxationTime in order to have sampled the 
	 * zipping+relaxation process for at least minZippingSamplingTime. */
	double minZippingSamplingTime;

	/* Write the state of the world at the time of zipping confirmation 
	 * to this file. Do nothing if this is NULL. */
	const char *zippedStateFile;
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


typedef struct {
	/* A pair is 'bound' if its base pair potential is lower than the 
	 * given threshold. */
	double energyThreshold;

	/* Set the temperature to this value at the start of sampling. Use 
	 * a negative value here if you don't want to change the 
	 * temperature. */
	double temperature;

	/* If this is true: dump the distance between base pairs instead of 
	 * just whether they are bound or not according to energyThreshold. */
	/* TODO this is redundant with basePairingSampler! Use that for the 
	 * sampling? */
	bool dumpDistances;
} HairpinStateSamplerConfig;

Sampler hairpinStateSampler(HairpinStateSamplerConfig *hssc);

/* A sampler that, for each monomer, dumps the magnitude of the force |F| 
 * and the velocity |v| of (the COM of) the monomer, as well as an estimate 
 * of the 'friction' on the monomer equal to (F dot v)/(|F|^2).  */
Sampler forceVelFricSampler(Strand *strand);
/* Analogous, but for the backbone phosphates only */
Sampler forceVelFricPSampler(Strand *strand);

/* A trivial sampler that does nothing. Useful for debugging purposes. */
Sampler trivialSampler(void);

/* A trivial sampler that does nothing, except for requesting a quit after 
 * the given time. Useful for debugging purposes. */
Sampler trivialStoppingSampler(double stopTime);

#endif


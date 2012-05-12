#include "samplers.h"
#include "physics.h"
#include "spgrid.h"

/* Simple sampler start that just passes the configuration data as the 
 * state pointer. */
static void *passConf(void *conf)
{
	return conf;
}
/* Simple sampler stop that just frees the state. */
static void freeState(SamplerData *sd, void *state)
{
	UNUSED(sd);
	free(state);
}


/* TEMPERATURE */

static void *avgTempStart(void *conf)
{
	UNUSED(conf);
	double *accum = malloc(sizeof *accum); //yes, this is a silly malloc :P
	return accum;
}
static bool avgTempSample(SamplerData *sd, void *data)
{
	UNUSED(sd);
	double *accum = (double*) data;
	*accum += temperature();
	return true;
}
static void avgTempStop(SamplerData *sd, void *data)
{
	double *accum = (double*) data;
	printf("Average temperature: %f\n", *accum / sd->sample);
	free(accum);
}
Sampler averageTemperatureSampler(void)
{
	Sampler sampler;
	sampler.samplerConf = NULL;
	sampler.start  = &avgTempStart;
	sampler.sample = &avgTempSample;
	sampler.stop   = &avgTempStop;
	return sampler;
}




/* STATS / VERBOSE */

static bool dumpStatsSample(SamplerData *sd, void *data)
{
	UNUSED(sd);
	UNUSED(data);
	dumpStats();
	return true;
}
Sampler dumpStatsSampler(void)
{
	Sampler sampler;
	sampler.samplerConf = NULL;
	sampler.start  = NULL;
	sampler.sample = &dumpStatsSample;
	sampler.stop   = NULL;
	return sampler;
}



/* POSITION / CENTER OF MASS */

typedef struct
{
	Particle *ps;
	int num;
} ParticlesCOMSamplerConf;
static bool particlesCOMSample(SamplerData *sd, void *data)
{
	UNUSED(sd);
	ParticlesCOMSamplerConf *pcsc = (ParticlesCOMSamplerConf*) data;
	Vec3 COM = getCOM(pcsc->ps, pcsc->num);
	printVectorExp(COM);
	printf("\n");
	return true;
}
static Sampler particlesCOMSampler(Particle *ps, int num)
{
	ParticlesCOMSamplerConf *pcsc = malloc(sizeof(*pcsc));
	pcsc->ps = ps;
	pcsc->num = num;

	Sampler sampler;
	sampler.samplerConf = pcsc;
	sampler.start  = &passConf;
	sampler.sample = &particlesCOMSample;
	sampler.stop   = &freeState;
	return sampler;
}

Sampler particlePositionSampler(Particle *p)
{
	return particlesCOMSampler(p, 1);
}
Sampler strandCOMSampler(Strand *s)
{
	return particlesCOMSampler(s->all, 3 * s->numMonomers);
}



/* SQUARED DISPLACEMENT */

typedef struct
{
	Particle *ps;
	int num;
	Vec3 initialPos;
	double initialTime;
} SquaredDisplacementConf;

static void *particlesSquaredDisplacementStart(void *conf)
{
	SquaredDisplacementConf *sdc = (SquaredDisplacementConf*) conf;
	sdc->initialPos = getCOM(sdc->ps, sdc->num);
	sdc->initialTime = getTime();
	return sdc;
}
static bool particlesSquaredDisplacementSample(SamplerData *sd, void *state)
{
	UNUSED(sd);
	SquaredDisplacementConf *sdc = (SquaredDisplacementConf*) state;

	Vec3 COM = getCOM(sdc->ps, sdc->num);
	Vec3 displacement = sub(COM, sdc->initialPos);
	double squaredDisplacement = length2(displacement);
	//printf("%e %e", getTime() - sdc->initialTime, squaredDisplacement);
	printf("%e ",getTime());
	printVectorExp(COM);
	printf("\n");

	if (squaredDisplacement > 0.25 * SQUARE(world.worldSize)) {
		fprintf(stderr, "WARNING! displacement > 0.5*worldsize. "
				"Probably errors due to periodic boundary "
				"conditions! Bailing out!\n");
		return false;
	}


	return true;
}
static Sampler particlesSquaredDisplacementSampler(Particle *ps, int num)
{
	SquaredDisplacementConf *sdc = malloc(sizeof(*sdc));
	sdc->ps = ps;
	sdc->num = num;

	Sampler sampler;
	sampler.samplerConf = sdc;
	sampler.start  = &particlesSquaredDisplacementStart;
	sampler.sample = &particlesSquaredDisplacementSample;
	sampler.stop   = &freeState;
	return sampler;
}

Sampler particleSquaredDisplacementSampler(Particle *p)
{
	return particlesSquaredDisplacementSampler(p, 1);
}
Sampler strandCOMSquaredDisplacementSampler(Strand *s)
{
	return particlesSquaredDisplacementSampler(s->all, 3 * s->numMonomers);
}



/* BASE PAIRING */

typedef struct
{
	int count;
	double threshold;
} BasePairingCounterData;
static void basePairingCounter(Particle *p1, Particle *p2, void *data)
{
	BasePairingCounterData *bpcd = (BasePairingCounterData*) data;
	double V = VbasePair(p1, p2);
	if (V < bpcd->threshold)
		bpcd->count++;
}
static bool basePairingSample(SamplerData *sd, void *state)
{
	double *threshold = (double*) state;

	/* All base pairs */
	BasePairingCounterData bpcd;
	bpcd.count = 0;
	bpcd.threshold = *threshold;
	forEveryPairD(&basePairingCounter, &bpcd);

	/* Matching base pairs if two equal-length strands in the world 
	 * (assumed to be complementary) */
	int correctlyBound = -1; /* guard */
	if (world.numStrands == 2  &&  (world.strands[0].numMonomers
					== world.strands[1].numMonomers)) {
		correctlyBound = 0;
		int n = world.strands[0].numMonomers;
		for (int i = 0; i < n; i++) {
			int j = n - 1 - i; //TODO
			double V = VbasePair(&world.strands[0].Bs[i],
					     &world.strands[1].Bs[j]);
			if (V < *threshold) {
				correctlyBound++;
				printf("1 ");
			} else {
				printf("0 ");
			}
		}
	}

	if (correctlyBound >= 0) {
		printf("%e\t%d\t%d\n", getTime(), bpcd.count, correctlyBound);
		if(sd->string != NULL)
			snprintf(sd->string, sd->strBufSize,
					"All BPs: %d, Correct BPs: %d",
					bpcd.count, correctlyBound);
	} else {
		printf("%e\t%d\n", getTime(), bpcd.count);
		if(sd->string != NULL)
			snprintf(sd->string, sd->strBufSize,
					"All BPs: %d", bpcd.count);
	}

	return true;
}
Sampler basePairingSampler(double energyThreshold)
{
	Sampler sampler;
	double *thresholdCpy = malloc(sizeof(*thresholdCpy));
	*thresholdCpy = energyThreshold;

	sampler.samplerConf = thresholdCpy;
	sampler.start  = &passConf;
	sampler.sample = &basePairingSample;
	sampler.stop   = &freeState;
	return sampler;
}


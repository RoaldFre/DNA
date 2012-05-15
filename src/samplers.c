#include "samplers.h"
#include "physics.h"
#include "spgrid.h"
#include <string.h>

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
static SamplerSignal avgTempSample(SamplerData *sd, void *data)
{
	UNUSED(sd);
	double *accum = (double*) data;
	*accum += temperature();
	return SAMPLER_OK;
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

static SamplerSignal dumpStatsSample(SamplerData *sd, void *data)
{
	UNUSED(sd);
	UNUSED(data);
	dumpStats();
	return SAMPLER_OK;
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
static SamplerSignal particlesCOMSample(SamplerData *sd, void *data)
{
	UNUSED(sd);
	ParticlesCOMSamplerConf *pcsc = (ParticlesCOMSamplerConf*) data;
	Vec3 COM = getCOM(pcsc->ps, pcsc->num);
	printVectorExp(COM);
	printf("\n");
	return SAMPLER_OK;
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
static SamplerSignal particlesSquaredDisplacementSample(SamplerData *sd, void *state)
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
		return SAMPLER_ERROR;
	}


	return SAMPLER_OK;
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
static void* basePairingStart(void *conf)
{
	BasePairingConfig *bpc = (BasePairingConfig*) conf;
	config.thermostatTemp = bpc->T;
	return bpc;
}
static SamplerSignal basePairingSample(SamplerData *sd, void *state)
{
	BasePairingConfig *bpc = (BasePairingConfig*) state;

	/* All base pairs */
	BasePairingCounterData bpcd;
	bpcd.count = 0;
	bpcd.threshold = bpc->energyThreshold;
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
			if (V < bpc->energyThreshold) {
				correctlyBound++;
				printf("1 ");
			} else {
				printf("0 ");
			}
		}
	}

	/* Matching base pairs if one strands in the world (assumed to be a 
	 * hairpin) */
	if (world.numStrands == 1) {
		correctlyBound = 0;
		int n = world.strands[0].numMonomers;
		for (int i = 0; i < n/2; i++) {
			int j = n - 1 - i;
			double V = VbasePair(&world.strands[0].Bs[i],
					     &world.strands[0].Bs[j]);
			if (V < bpc->energyThreshold) {
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

	return SAMPLER_OK;
}
Sampler basePairingSampler(BasePairingConfig *bpc)
{
	Sampler sampler;
	BasePairingConfig *bpcCopy = malloc(sizeof(*bpcCopy));
	memcpy(bpcCopy, bpc, sizeof(*bpcCopy));

	sampler.samplerConf = bpcCopy;
	sampler.start  = &basePairingStart;
	sampler.sample = &basePairingSample;
	sampler.stop   = &freeState;
	return sampler;
}





/* HAIR PINS */
typedef struct {
	enum {
		WAITING_TO_FORM,
		WAITING_FOR_CONFIRMATION,
		START_RELAXATION,
		WAITING_TO_RELAX,
		MEASURING,
	} status;
	double confirmationStartTime;
	double relaxStartTime;
	double measureStartTime;
	int currentStep;
	HairpinSamplerConfig conf;
} HairpinSamplerData;

static void* hairpinStart(void *conf)
{
	HairpinSamplerConfig *hsc = (HairpinSamplerConfig*) conf;
	HairpinSamplerData *hsd = malloc(sizeof(*hsd));
	hsd->conf = *hsc;
	hsd->status = WAITING_TO_FORM;
	return hsd;
}
static SamplerSignal hairpinSample(SamplerData *sd, void *state)
{
	UNUSED(sd);
	HairpinSamplerData *hsd = (HairpinSamplerData*) state;
	HairpinSamplerConfig *hsc = &hsd->conf;

	int correctlyBound = 0;
	int n = world.strands[0].numMonomers;
	for (int i = 0; i < n/2; i++) {
		int j = n - 1 - i;
		double V = VbasePair(&world.strands[0].Bs[i],
				     &world.strands[0].Bs[j]);
		if (V < hsc->energyThreshold)
			correctlyBound++;
	}

	int requiredBounds = n/2 - hsc->allowedUnboundBPs;
	double time = getTime();

	switch (hsd->status) {
	case WAITING_TO_FORM:
		if (correctlyBound >= requiredBounds) {
			hsd->confirmationStartTime = time;
			hsd->status = WAITING_FOR_CONFIRMATION;
			printf("# Reached binding threshold of %d base pairs at %e\n",
					requiredBounds, time);
		}
		break;
	case WAITING_FOR_CONFIRMATION:
		if (correctlyBound < requiredBounds) {
			/* It was a dud! */
			hsd->status = WAITING_TO_FORM;
			printf("# Could not confirm, waiting to form again at %e\n",
					time);
			break;
		}
		if (time - hsd->confirmationStartTime
						< hsc->confirmationTime)
			break; /* need to wait for confirmation */

		printf("# Confirmd, at %e\n", time);
		/* We have confirmation: start relaxation phase */
		/* FALL THROUGH! */
	case START_RELAXATION:
		hsd->relaxStartTime = time;
		hsd->status = WAITING_TO_RELAX;
		/* Set temperature */
		config.thermostatTemp = hsc->Tstart
				+ hsd->currentStep * hsc->Tstep;
		printf("# Setting temperature to %f, step %d\n",
				config.thermostatTemp, hsd->currentStep);
		break;
	case WAITING_TO_RELAX:
		if (time - hsd->relaxStartTime < hsc->relaxationTime)
			break;
		printf("# Relaxation done at %e, step %d\n",
				time, hsd->currentStep);
		hsd->measureStartTime = time;
		hsd->status = MEASURING;
		break;
	case MEASURING:
		if (time - hsd->measureStartTime >= hsc->measureTime) {
			/* End of this measurement run */
			hsd->currentStep++;
			if (hsd->currentStep >= hsc->numSteps) {
				printf("# Sampled to the end! :-)");
				return SAMPLER_STOP; /* All done! */
			}

			hsd->status = START_RELAXATION;
			break;
		}
		printf("%e %d ", time, correctlyBound);
		/* print (half of) the config */
		for (int i = 0; i < n/2; i++) {
			int j = n - 1 - i;
			double V = VbasePair(&world.strands[0].Bs[i],
					     &world.strands[0].Bs[j]);
			if (V < hsc->energyThreshold) {
				correctlyBound++;
				printf("1 ");
			} else {
				printf("0 ");
			}
		}
		printf("\n");
		break;
	}
	return SAMPLER_OK;
}
Sampler hairpinSampler(HairpinSamplerConfig *hsc)
{
	Sampler sampler;
	HairpinSamplerConfig *hpcCopy = malloc(sizeof(*hpcCopy));
	memcpy(hpcCopy, hsc, sizeof(*hpcCopy));

	sampler.samplerConf = hpcCopy;
	sampler.start  = &hairpinStart;
	sampler.sample = &hairpinSample;
	sampler.stop   = &freeState;
	return sampler;
}


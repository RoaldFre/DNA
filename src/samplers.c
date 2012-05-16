#include <string.h>
#include "samplers.h"
#include "physics.h"
#include "spgrid.h"
#include "octave.h"

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
	*accum = 0;
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
	memset(pcsc, 0, sizeof(*pcsc));
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
	memset(sdc, 0, sizeof(*sdc));
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
	double startTime;
	double confirmationStartTime;
	double relaxStartTime;
	double measureStartTime;
	int measureIteration;
	int numMeasureIterations;
	int currentStep;
	HairpinFormationSamplerConfig conf;
} HairpinFormationSamplerData;

static void* hairpinFormationStart(void *conf)
{
	HairpinFormationSamplerConfig *hfsc = (HairpinFormationSamplerConfig*) conf;
	HairpinFormationSamplerData *hfsd = malloc(sizeof(*hfsd));
	memset(hfsd, 0, sizeof(*hfsd));

	hfsd->conf = *hfsc; /* struct copy */
	hfsd->status = WAITING_TO_FORM;
	hfsd->currentStep = 0;
	hfsd->startTime = getTime();

	octaveScalar("sampleStartTime", getTime());

	return hfsd;
}
static SamplerSignal hairpinFormationSample(SamplerData *sd, void *state)
{
	HairpinFormationSamplerData *hfsd = (HairpinFormationSamplerData*) state;
	HairpinFormationSamplerConfig *hfsc = &hfsd->conf;

	int correctlyBound = 0;
	int n = world.strands[0].numMonomers;
	for (int i = 0; i < n/2; i++) {
		int j = n - 1 - i;
		double V = VbasePair(&world.strands[0].Bs[i],
				     &world.strands[0].Bs[j]);
		if (V < hfsc->energyThreshold)
			correctlyBound++;
	}

	int requiredBounds = n/2 - hfsc->allowedUnboundBPs;
	double time = getTime();

	switch (hfsd->status) {
	case WAITING_TO_FORM:
		if (correctlyBound >= requiredBounds) {
			hfsd->confirmationStartTime = time;
			hfsd->status = WAITING_FOR_CONFIRMATION;
			octaveComment("Reached binding threshold of %d base "
					"pairs at %e after a time %e\n",
					requiredBounds, time, time - hfsd->startTime);
		}
		break;
	case WAITING_FOR_CONFIRMATION:
		if (correctlyBound < requiredBounds) {
			/* It was a dud! */
			hfsd->status = WAITING_TO_FORM;
			octaveComment("Could not confirm, waiting to form "
					"again at %e\n", time);
			break;
		}
		if (time - hfsd->confirmationStartTime
						< hfsc->confirmationTime)
			break; /* need to wait for confirmation */

		/* We have confirmation: start relaxation phase */

		/* Set this once here, so all steps have the exact same number of 
		 * iterations. It's a PITA to do data processing otherwise. */
		hfsd->numMeasureIterations = hfsc->measureTime / sd->sampleInterval;

		/* Dump all info before starting 3D matrix */
		octaveComment("Confirmed at %e after a time since initial "
				"threshold %e\n", time,
				time - hfsd->startTime - hfsc->confirmationTime);

		octaveScalar("timeTillZipping",
				time - hfsd->startTime - hfsc->confirmationTime);
		octaveScalar("temperature", config.thermostatTemp);

		octaveScalar("confirmationTime",  hfsc->confirmationTime);
		octaveScalar("relaxationTime",    hfsc->relaxationTime);
		octaveScalar("measureTime",       hfsc->measureTime);
		octaveScalar("allowedUnboundBPs", hfsc->allowedUnboundBPs);
		octaveScalar("timestep",          config.timeStep);
		octaveScalar("numMonomers",       n);

		octaveMatrixHeader("temperatures", hfsc->numSteps, 1);
		for (int i = 0; i < hfsc->numSteps; i++)
			printf("%e\n", hfsc->Tstart + i * hfsc->Tstep);

		/* Start the header for the 3D matrix of data */
		octaveComment("Series of 2D matrices of the form");
		octaveComment("data(:,:,i)' =");
		octaveComment(" [ time(1)  numBound(1)  lastBasePair(1)  ...  middleBasePair(1)");
		octaveComment("   time(2)  numBound(2)  lastBasePair(2)  ...  middleBasePair(2)");
		octaveComment("    ...        ...            ...         ...        ...        ");
		octaveComment("   time(n)  numBound(n)  lastBasePair(n)  ...  middleBasePair(n) ]");
		octaveComment("NOTE THE TRANSPOSED SIGN at the end of data(:,:,)'  <-- !!");
		octaveComment("Each such 2D matix i is measured with the corresponding");
		octaveComment("temperature in the array 'temperatures': temperatures(i).");
		octave3DMatrixHeader("data",
				2 + n/2,
				hfsd->numMeasureIterations,
				hfsc->numSteps);

		/* INTENTIONAL FALL THROUGH! */
	case START_RELAXATION:
		hfsd->relaxStartTime = time;
		hfsd->status = WAITING_TO_RELAX;
		/* Set temperature */
		config.thermostatTemp = hfsc->Tstart
				+ hfsd->currentStep * hfsc->Tstep;
		break;
	case WAITING_TO_RELAX:
		if (time - hfsd->relaxStartTime < hfsc->relaxationTime)
			break;
		/* Relaxation done */
		hfsd->measureStartTime = time;
		hfsd->measureIteration = 0;
		hfsd->status = MEASURING;
		break;
	case MEASURING:
		hfsd->measureIteration++;
		if (hfsd->measureIteration > hfsd->numMeasureIterations) {
			/* End of this measurement run */
			hfsd->currentStep++;
			if (hfsd->currentStep >= hfsc->numSteps) {
				octaveComment("# Sampled to the end! :-)");
				return SAMPLER_STOP; /* All done! */
			}

			hfsd->status = START_RELAXATION;
			break;
		}
		printf("%e\n%d\n", time, correctlyBound);
		/* print (half of) the config */
		for (int i = 0; i < n/2; i++) {
			int j = n - 1 - i;
			double V = VbasePair(&world.strands[0].Bs[i],
					     &world.strands[0].Bs[j]);
			if (V < hfsc->energyThreshold) {
				correctlyBound++;
				printf("1\n");
			} else {
				printf("0\n");
			}
		}
		printf("\n");
		break;
	}
	return SAMPLER_OK;
}
Sampler hairpinFormationSampler(HairpinFormationSamplerConfig *hfsc)
{
	Sampler sampler;
	HairpinFormationSamplerConfig *hpcCopy = malloc(sizeof(*hpcCopy));
	memcpy(hpcCopy, hfsc, sizeof(*hpcCopy));

	sampler.samplerConf = hpcCopy;
	sampler.start  = &hairpinFormationStart;
	sampler.sample = &hairpinFormationSample;
	sampler.stop   = &freeState;
	return sampler;
}


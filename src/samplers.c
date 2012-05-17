#include <string.h>
#include "samplers.h"
#include "physics.h"
#include "spgrid.h"
#include "octave.h"

/* Simple sampler start that just passes the configuration data as the 
 * state pointer. */
static void *passConf(SamplerData *sd, void *conf)
{
	UNUSED(sd);
	return conf;
}
/* Simple sampler stop that just frees the state. */
static void freeState(SamplerData *sd, void *state)
{
	UNUSED(sd);
	free(state);
}


/* TEMPERATURE */

static void *avgTempStart(SamplerData *sd, void *conf)
{
	UNUSED(sd);
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

static void *particlesSquaredDisplacementStart(SamplerData *sd, void *conf)
{
	UNUSED(sd);
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
static void* basePairingStart(SamplerData *sd, void *conf)
{
	UNUSED(sd);
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
	int n = world.strands[0].numMonomers;

	/* Matching base pairs if two equal-length strands in the world 
	 * (assumed to be complementary) */
	int correctlyBound = -1; /* guard */
	if (world.numStrands == 2  &&  (world.strands[0].numMonomers
					== world.strands[1].numMonomers)) {
		correctlyBound = 0;
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





/* HAIR PIN FORMATION */
static int getCorrectlyBoundHairpinBasePairs(Strand *s, double energyThreshold)
{
	int correctlyBound = 0;
	int n = s->numMonomers;
	for (int i = 0; i < n/2; i++) {
		int j = n - 1 - i;
		double V = VbasePair(&s->Bs[i], &s->Bs[j]);
		if (V < energyThreshold)
			correctlyBound++;
	}
	return correctlyBound;
}

typedef struct {
	enum {
		WAITING_TO_FORM,
		WAITING_FOR_CONFIRMATION,
	} status;
	double startTime;
	double confirmationStartTime;
	HairpinFormationSamplerConfig conf;
} HairpinFormationSamplerData;

static void* hairpinFormationStart(SamplerData *sd, void *conf)
{
	UNUSED(sd);
	HairpinFormationSamplerConfig *hfc =
			(HairpinFormationSamplerConfig*) conf;
	HairpinFormationSamplerData *hfd = malloc(sizeof(*hfd));
	memset(hfd, 0, sizeof(*hfd));

	hfd->conf = *hfc; /* struct copy */
	hfd->status = WAITING_TO_FORM;
	hfd->startTime = getTime();

	octaveScalar("sampleStartTime", getTime());

	return hfd;
}
static SamplerSignal hairpinFormationSample(SamplerData *sd, void *state)
{
	HairpinFormationSamplerData *hfd =
			(HairpinFormationSamplerData*) state;
	HairpinFormationSamplerConfig *hfc = &hfd->conf;

	int correctlyBound = getCorrectlyBoundHairpinBasePairs(
				&world.strands[0], hfc->energyThreshold);
	int n = world.strands[0].numMonomers;
	int requiredBounds = n/2 - hfc->allowedUnboundBPs;

	if(sd->string != NULL)
		snprintf(sd->string, sd->strBufSize,
				"Correct hairpin BPs: %d, required: %d",
				correctlyBound, requiredBounds);

	double time = getTime();

	switch (hfd->status) {
	case WAITING_TO_FORM:
		if (correctlyBound >= requiredBounds) {
			hfd->confirmationStartTime = time;
			hfd->status = WAITING_FOR_CONFIRMATION;
			octaveComment("Reached binding threshold of %d base "
					"pairs at %e after a time %e\n",
					requiredBounds, time, time - hfd->startTime);
		}
		break;
	case WAITING_FOR_CONFIRMATION:
		if (correctlyBound < requiredBounds) {
			/* It was a dud! */
			hfd->status = WAITING_TO_FORM;
			octaveComment("Could not confirm, waiting to form "
					"again at %e\n", time);
			break;
		}
		if (time - hfd->confirmationStartTime
						< hfc->confirmationTime)
			break; /* need to wait for confirmation */

		/* We have confirmation! */
		octaveComment("Confirmed at %e after a time since initial "
				"threshold %e\n", time,
				time - hfd->startTime - hfc->confirmationTime);

		octaveScalar("timeTillZipping",
				time - hfd->startTime - hfc->confirmationTime);
		octaveScalar("timeAtZipping",
				time - hfc->confirmationTime);
		octaveScalar("temperature",       config.thermostatTemp);
		octaveScalar("allowedUnboundBPs", hfc->allowedUnboundBPs);
		octaveScalar("numMonomers",       n);
		octaveScalar("timestep",          config.timeStep);

		octaveComment("Successful end! :-)");
		return SAMPLER_STOP;
	default:
		fprintf(stderr, "Unknown status in hairpinFormationSample!\n");
		assert(false); return SAMPLER_ERROR;
	}
	return SAMPLER_OK;
}
Sampler hairpinFormationSampler(HairpinFormationSamplerConfig *hfc)
{
	Sampler sampler;
	HairpinFormationSamplerConfig *hfcCopy = malloc(sizeof(*hfcCopy));
	memcpy(hfcCopy, hfc, sizeof(*hfcCopy));

	sampler.samplerConf = hfcCopy;
	sampler.start  = &hairpinFormationStart;
	sampler.sample = &hairpinFormationSample;
	sampler.stop   = &freeState;
	return sampler;
}



/* HAIR PIN MELTING TEMPERATURE */
typedef struct {
	enum {
		START_RELAXATION,
		WAITING_TO_RELAX,
		MEASURING,
	} status;
	double startTime;
	double relaxStartTime;
	double measureStartTime;
	int measureIteration;
	int numMeasureIterations;
	int currentStep;
	long accumulatedBoundBasePairs;
	double *averageBPsPerStep;
	HairpinMeltingTempSamplerConfig conf;
} HairpinMeltingTempSamplerData;

static void* hairpinMeltingTempStart(SamplerData *sd, void *conf)
{
	HairpinMeltingTempSamplerConfig *hmtc = (HairpinMeltingTempSamplerConfig*) conf;
	HairpinMeltingTempSamplerData *hmtd = malloc(sizeof(*hmtd));
	memset(hmtd, 0, sizeof(*hmtd));

	hmtd->conf = *hmtc; /* struct copy */
	hmtd->status = WAITING_TO_FORM;
	hmtd->currentStep = 0;
	hmtd->startTime = getTime();

	/* Set this once here, so all steps have the exact same number of 
	 * iterations. It's a PITA to do data processing otherwise. */
	hmtd->numMeasureIterations = hmtc->measureTime / sd->sampleInterval;

	/* We need to accumulate these and only dump them at the end, 
	 * because we cannot interfere with the matrix stream that we dump 
	 * when 'hmtc->verbose' is true. */
	hmtd->averageBPsPerStep = calloc(sizeof(*hmtd->averageBPsPerStep), 
					hmtc->numSteps);
	int n = world.strands[0].numMonomers;
	octaveScalar("sampleStartTime",   getTime());
	octaveScalar("temperature",       config.thermostatTemp);
	octaveScalar("relaxationTime",    hmtc->relaxationTime);
	octaveScalar("measureTime",       hmtc->measureTime);
	octaveScalar("numMonomers",       n);
	octaveScalar("timestep",          config.timeStep);

	octaveMatrixHeader("temperatures", hmtc->numSteps, 1);
	for (int i = 0; i < hmtc->numSteps; i++)
		printf("%e\n", hmtc->Tstart + i * hmtc->Tstep);

	if (hmtc->verbose) {
		/* Start the header for the 3D matrix of data */
		octaveComment("Series of 2D matrices of the form");
		octaveComment("data(:,:,i)' =");
		octaveComment(" [  time(t1) = t1      time(t2) = t2     ...   time(tn) = tn");
		octaveComment("    numBounds(t1)      numBounds(t2)     ...   numBounds(tn)");
		octaveComment("   lastBasePair(t1)   lastBasePair(t2)   ...  lastBasePair(tn)");
		octaveComment("         ...               ...           ...         ...      ");
		octaveComment("   firstBasePair(tn)  firstBasePair(t2)  ...  firstBasePair(tn) ]");
		octaveComment("Each such 2D matix i is measured with the corresponding");
		octaveComment("temperature in the array 'temperatures': temperatures(i).");
		octave3DMatrixHeader("data",
				2 + n/2,
				hmtd->numMeasureIterations,
				hmtc->numSteps);
	}

	return hmtd;
}
static SamplerSignal hairpinMeltingTempSample(SamplerData *sd, void *state)
{
	HairpinMeltingTempSamplerData *hmtd = (HairpinMeltingTempSamplerData*) state;
	HairpinMeltingTempSamplerConfig *hmtc = &hmtd->conf;

	int correctlyBound = getCorrectlyBoundHairpinBasePairs(
				&world.strands[0], hmtc->energyThreshold);

	if(sd->string != NULL)
		snprintf(sd->string, sd->strBufSize,
				"Correct hairpin BPs: %d", correctlyBound);

	double time = getTime();

	switch (hmtd->status) {
	case START_RELAXATION:
		hmtd->relaxStartTime = time;
		hmtd->status = WAITING_TO_RELAX;
		/* Set temperature */
		config.thermostatTemp = hmtc->Tstart
				+ hmtd->currentStep * hmtc->Tstep;
		break;
	case WAITING_TO_RELAX:
		if (time - hmtd->relaxStartTime < hmtc->relaxationTime)
			break;
		/* Relaxation done */
		hmtd->measureStartTime = time;
		hmtd->measureIteration = 0;
		hmtd->accumulatedBoundBasePairs = 0;
		hmtd->status = MEASURING;
		/* Intentional fall through */
	case MEASURING:
		hmtd->measureIteration++;
		if (hmtd->measureIteration > hmtd->numMeasureIterations) {
			/* End of this measurement run */
			hmtd->averageBPsPerStep[hmtd->currentStep] = 
					hmtd->accumulatedBoundBasePairs 
						/ hmtd->numMeasureIterations;

			hmtd->currentStep++;
			if (hmtd->currentStep >= hmtc->numSteps) {
				/* We are finished! Don't forget to dump 
				 * the accumulated bound base pairs! */
				octaveMatrixHeader("averageBoundBasePairs", 
							hmtc->numSteps, 1);
				for (int i = 0; i < hmtc->numSteps; i++)
					printf("%e\n", hmtd->averageBPsPerStep[i]);

				octaveComment("Successful end! :-)");
				return SAMPLER_STOP; /* All done! */
			}

			hmtd->status = START_RELAXATION;
			break;
		}

		/* Not at the end of this step yet: just sample */
		hmtd->accumulatedBoundBasePairs += correctlyBound;
		if (hmtc->verbose) {
			printf("%e\n%d\n", time, correctlyBound);
			/* print (half of) the config */
			int n = world.strands[0].numMonomers;
			for (int i = 0; i < n/2; i++) {
				int j = n - 1 - i;
				double V = VbasePair(&world.strands[0].Bs[i],
						     &world.strands[0].Bs[j]);
				if (V < hmtc->energyThreshold) {
					printf("1\n");
				} else {
					printf("0\n");
				}
			}
			printf("\n");
		}
		break;
	}
	return SAMPLER_OK;
}
Sampler hairpinMeltingTempSampler(HairpinMeltingTempSamplerConfig *hmtc)
{
	Sampler sampler;
	HairpinMeltingTempSamplerConfig *hmtcCopy = malloc(sizeof(*hmtcCopy));
	memcpy(hmtcCopy, hmtc, sizeof(*hmtcCopy));

	sampler.samplerConf = hmtcCopy;
	sampler.start  = &hairpinMeltingTempStart;
	sampler.sample = &hairpinMeltingTempSample;
	sampler.stop   = &freeState;
	return sampler;
}


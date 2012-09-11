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

static SamplerSignal tempSample(SamplerData *sd, void *data)
{
	UNUSED(sd);
	UNUSED(data);
	printf("%e %f\n", getTime(), getKineticTemperature());
	return SAMPLER_OK;
}
Sampler temperatureSampler(void)
{
	Sampler sampler;
	sampler.samplerConf = NULL;
	sampler.start  = NULL;
	sampler.sample = &tempSample;
	sampler.stop   = NULL;
	return sampler;
}


/* AVERAGE TEMPERATURE */

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
	*accum += getKineticTemperature();
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



/* END TO END DISTANCE */

typedef struct
{
	Strand *strand;
} EndToEndDistConf;

static SamplerSignal endToEndDistSample(SamplerData *sd, void *state)
{
	UNUSED(sd);
	EndToEndDistConf *etedc = (EndToEndDistConf*) state;
	Strand *s = etedc->strand;

	double endToEndDist = distance(
			getMonomerCOM(s, 0),
			getMonomerCOM(s, s->numMonomers - 1));
	printf("%e\t%e\n", getTime(), endToEndDist);

	return SAMPLER_OK;
}
Sampler endToEndDistSampler(Strand *strand)
{
	EndToEndDistConf *etedc = malloc(sizeof(*etedc));
	memset(etedc, 0, sizeof(*etedc));
	etedc->strand = strand;

	Sampler sampler;
	sampler.samplerConf = etedc;
	sampler.start  = &passConf;
	sampler.sample = &endToEndDistSample;
	sampler.stop   = &freeState;
	return sampler;
}





/* HAIR PIN FORMATION */
/* Print (half of) the hairpin config and return  the number of correctly 
 * bound base pairs. */
static int dumpHairpinState(Strand *s, double energyThreshold)
{
	int correctlyBound = 0;
	int n = s->numMonomers;
	for (int i = 0; i < n/2; i++) {
		int j = n - 1 - i;
		double V = VbasePair(&s->Bs[i], &s->Bs[j]);
		if (V < energyThreshold) {
			correctlyBound++;
			printf("1 ");
		} else {
			printf("0 ");
		}
	}
	return correctlyBound;
}
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
		WAITING_TO_ZIP,
		WAITING_FOR_ZIPPING_CONFIRMATION,
		RELAXATION_IN_ZIPPED_STATE,
		WAITING_TO_UNZIP,
		WAITING_FOR_UNZIPPING_CONFIRMATION,
	} status;
	double zippingPhaseStartTime;
	double unzippingPhaseStartTime;
	double confirmationStartTime; /* Used for both zipping & unzipping */
	double timeOfConfirmation; /* Used for both zipping & unzipping */
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
	hfd->status = WAITING_TO_ZIP;
	hfd->zippingPhaseStartTime = getTime();

	/* Set temperature for zipping */
	config.thermostatTemp = hfc->zippingTemperature;

	/* Dump info */
	int n = world.strands[0].numMonomers;
	octaveScalar("sampleStartTime",         getTime());
	octaveScalar("zippingPhaseStartTime",   getTime());
	octaveScalar("zippingTemperature",      hfc->zippingTemperature);
	octaveScalar("unzippingTemperature",    hfc->unzippingTemperature);
	octaveScalar("requiredBoundBPs",        hfc->requiredBoundBPs);
	octaveScalar("allowedBoundBPs",         hfc->allowedBoundBPs);
	octaveScalar("numMonomers",             n);
	octaveScalar("timestep",                config.timeStep);
	octaveScalar("sampleInterval",          sd->sampleInterval);

	return hfd;
}
static SamplerSignal hairpinFormationSample(SamplerData *sd, void *state)
{
	HairpinFormationSamplerData *hfd =
			(HairpinFormationSamplerData*) state;
	HairpinFormationSamplerConfig *hfc = &hfd->conf;

	int correctlyBound = getCorrectlyBoundHairpinBasePairs(
				&world.strands[0], hfc->energyThreshold);
	int requiredBounds = hfc->requiredBoundBPs; /* for zipping */
	int allowedBounds = hfc->allowedBoundBPs; /* for unzipping */

	if(sd->string != NULL)
		snprintf(sd->string, sd->strBufSize,
				"Correct hairpin BPs: %d, thresholds: %d and %d",
				correctlyBound, allowedBounds, requiredBounds);

	double time = getTime();

	switch (hfd->status) {
	case WAITING_TO_ZIP:
		octaveStartComment();
		printf("[waiting to zip] %e %d ", time, correctlyBound);
		dumpHairpinState(&world.strands[0], hfc->energyThreshold);
		octaveEndComment();
		if (correctlyBound < requiredBounds)
			break; /* Keep waiting */
		
		/* We have detected initial zipping! */
		hfd->confirmationStartTime = time;
		hfd->status = WAITING_FOR_ZIPPING_CONFIRMATION;
		octaveComment("Reached zipping binding threshold of %d base "
				"pairs at %e after a time %e", requiredBounds,
				time, time - hfd->zippingPhaseStartTime);
		/* Intentional fall through */
	case WAITING_FOR_ZIPPING_CONFIRMATION:
		octaveStartComment();
		printf("[waiting to zip] %e %d ", time, correctlyBound);
		dumpHairpinState(&world.strands[0], hfc->energyThreshold);
		octaveEndComment();
		if (correctlyBound < requiredBounds) {
			/* It was a dud! */
			hfd->status = WAITING_TO_ZIP;
			octaveComment("Could not confirm zipping, waiting to "
					"zip again at %e", time);
			break;
		}
		if (time - hfd->confirmationStartTime
						< hfc->confirmationTime)
			break; /* Need to wait for confirmation */

		/* We have zipping confirmation! */
		double timeTillZipping = time - hfd->zippingPhaseStartTime
						- hfc->confirmationTime;
		octaveComment("Confirmed zipping at %e", time);
		octaveScalar("timeTillZipping", timeTillZipping);
		octaveComment("Starting relaxation phase in zipped state "
				"at %e", time);
		hfd->timeOfConfirmation = time;
		hfd->status = RELAXATION_IN_ZIPPED_STATE;
		/* Intentional fall through */
	case RELAXATION_IN_ZIPPED_STATE:
		octaveStartComment();
		printf("[zipped relaxation] %e %d ", time, correctlyBound);
		dumpHairpinState(&world.strands[0], hfc->energyThreshold);
		octaveEndComment();
		if (time - hfd->timeOfConfirmation < hfc->zippedRelaxationTime)
			break; /* Relax further */
		/* End of relaxation phase. Go to WAITING_TO_UNZIP if we 
		 * are still zipped, otherwise, wait until we are zipped 
		 * again! */
		if (correctlyBound < requiredBounds) {
			/* Not zipped! */
			octaveComment("Relaxation in zipped state: passed "
					"zippedRelaxationTime but not zipped "
					"anymore at: %e -- total relax time: %e",
					time, time - hfd->timeOfConfirmation);
			break; /* Wait to fully zip again */
		}

		/* We are still zipped! Go to next phase. */
		config.thermostatTemp = hfc->unzippingTemperature;
		octaveScalar("unzippingPhaseStartTime", time);
		hfd->unzippingPhaseStartTime = time;
		hfd->status = WAITING_TO_UNZIP;
		break;
	case WAITING_TO_UNZIP:
		octaveStartComment();
		printf("[waiting to unzip] %e %d ", time, correctlyBound);
		dumpHairpinState(&world.strands[0], hfc->energyThreshold);
		octaveEndComment();
		if (correctlyBound > allowedBounds)
			break;

		/* We have detected initial unzipping! */
		hfd->confirmationStartTime = time;
		hfd->status = WAITING_FOR_UNZIPPING_CONFIRMATION;
		octaveComment("Reached unzipping binding threshold of %d base "
				"pairs at %e after a time %e", allowedBounds,
				time, time - hfd->unzippingPhaseStartTime);
		/* Intentional fall through */
	case WAITING_FOR_UNZIPPING_CONFIRMATION:
		octaveStartComment();
		printf("[waiting to unzip] %e %d ", time, correctlyBound);
		dumpHairpinState(&world.strands[0], hfc->energyThreshold);
		octaveEndComment();
		if (correctlyBound > allowedBounds) {
			/* It was a dud! */
			hfd->status = WAITING_TO_UNZIP;
			octaveComment("Could not confirm unzipping, waiting "
					"to unzip again at %e", time);
			break;
		}
		if (time - hfd->confirmationStartTime
						< hfc->confirmationTime)
			break; /* need to wait for confirmation */

		/* We have unzipping confirmation! */
		double timeTillUnzipping = time - hfd->unzippingPhaseStartTime
						- hfc->confirmationTime;
		octaveComment("Confirmed unzipping at %e", time);
		octaveScalar("timeTillUnzipping", timeTillUnzipping);

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
	hmtd->status = START_RELAXATION;
	hmtd->currentStep = 0;

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
	octaveScalar("sampleInterval",    sd->sampleInterval);

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
		octaveComment("Each such 2D matrix i is measured with the corresponding");
		octaveComment("temperature in the array 'temperatures': temperatures(i).");
		octave3DMatrixHeader("hairpinState",
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
				((double) hmtd->accumulatedBoundBasePairs)
					/ (double) hmtd->numMeasureIterations;

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
			printf("%e ", time);
			dumpHairpinState(&world.strands[0], 
					hmtc->energyThreshold);
			printf("\n");
		}
		break;
	default:
		fprintf(stderr, "Unknown status in hairpinMeltingTempSample!\n");
		assert(false); return SAMPLER_ERROR;
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




/* BASE PAIRING */

typedef struct
{
	int count;
	double threshold;
} BasePairingCounterData;
static int dumpDualStrandState(Strand *s1, Strand *s2, double energyThreshold)
{
	assert(s1->numMonomers == s2->numMonomers);
	int n = s1->numMonomers;
	int correctlyBound = 0;
	for (int i = 0; i < n; i++) {
		int j = n - 1 - i; //TODO
		double V = VbasePair(&s1->Bs[i],
				     &s2->Bs[j]);
		if (V < energyThreshold) {
			correctlyBound++;
			printf("1 ");
		} else {
			printf("0 ");
		}
	}
	return correctlyBound;
}
static void basePairingCounter(Particle *p1, Particle *p2, void *data)
{
	BasePairingCounterData *bpcd = (BasePairingCounterData*) data;
	double V = VbasePair(p1, p2);
	if (V < bpcd->threshold)
		bpcd->count++;
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
		correctlyBound = dumpDualStrandState(&world.strands[0], 
				&world.strands[1], bpc->energyThreshold);
	}
	/* Matching base pairs if one strands in the world (assumed to be a 
	 * hairpin) */
	if (world.numStrands == 1) {
		correctlyBound = dumpHairpinState(&world.strands[0], 
				bpc->energyThreshold);
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
	sampler.start  = &passConf;
	sampler.sample = &basePairingSample;
	sampler.stop   = &freeState;
	return sampler;
}



/* TRIVIAL SAMPLER */

Sampler trivialSampler(void) {
	Sampler sampler = {
			.samplerConf = NULL,
			.start = NULL,
			.sample = NULL,
			.stop = NULL
	};
	return sampler;
}


#include "main.h" //TODO for LENGTH_FACTOR
#include "samplers.h"
#include "physics.h"

/* Simple sampler start that just passes the configuration data as the 
 * state pointer. */
static void *passConf(void *conf)
{
	return conf;
}
/* Simple sampler stop that just frees the state. */
static void freeState(long n, void *state)
{
	UNUSED(n);
	free(state);
}

void *avgTempStart(void *conf);
bool avgTempSample(long i, void *data);
void avgTempStop(long n, void *data);

Sampler averageTemperatureSampler(void)
{
	Sampler sampler;
	sampler.samplerConf = NULL;
	sampler.start  = &avgTempStart;
	sampler.sample = &avgTempSample;
	sampler.stop   = &avgTempStop;
	return sampler;
}

void *avgTempStart(void *conf)
{
	UNUSED(conf);
	double *accum = malloc(sizeof *accum); //yes, this is a silly malloc :P
	return accum;
}
bool avgTempSample(long i, void *data)
{
	UNUSED(i);
	double *accum = (double*) data;
	*accum += temperature();
	return true;
}
void avgTempStop(long n, void *data)
{
	double *accum = (double*) data;
	printf("Average temperature: %f\n", *accum / n);
	free(accum);
}



bool dumpStatsSample(long i, void *data);

Sampler dumpStatsSampler(void)
{
	Sampler sampler;
	sampler.samplerConf = NULL;
	sampler.start  = NULL;
	sampler.sample = &dumpStatsSample;
	sampler.stop   = NULL;
	return sampler;
}
bool dumpStatsSample(long i, void *data)
{
	UNUSED(i);
	UNUSED(data);
	dumpStats();
	return true;
}



static Sampler particlesCOMSampler(Particle *ps, int num);
static bool particlesCOMSample(long i, void *data);

Sampler particlePositionSampler(Particle *p)
{
	return particlesCOMSampler(p, 1);
}
Sampler strandCOMSampler(Strand *s)
{
	return particlesCOMSampler(s->all, 3*s->numMonomers);
}

typedef struct particlesCOMSamplerConf
{
	Particle *ps;
	int num;
} ParticlesCOMSamplerConf;

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

static bool particlesCOMSample(long i, void *data)
{
	UNUSED(i);
	ParticlesCOMSamplerConf *pcsc = (ParticlesCOMSamplerConf*) data;
	Vec3 COM = getCOM(pcsc->ps, pcsc->num);
	scale(&COM, 1 / LENGTH_FACTOR, &COM);
	printVector(&COM);
	printf("\n");
	return true;
}


#include "samplers.h"
#include "physics.h"

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


#include "measure.h"
#include "main.h" //for TIME_FACTOR, TODO put somewhere else
#include <string.h>

void *measStart(void *initialData);
bool measTick(void *state);
void measStop(void *state);

void *samplerStart(Sampler *sampler);
bool samplerSample(Sampler *sampler, long i, void *state);
void samplerStop(Sampler *sampler, long n, void *state);

Task measurementTask(Measurement *measurement)
{
	/* We must make a copy because given pointer is not guaranteed to 
	 * remain valid. */
	Measurement *measCpy = malloc(sizeof(*measCpy));
	memcpy(measCpy, measurement, sizeof(*measCpy));

	/* Make task */
	Task task;
	task.initialData = measCpy;
	task.start = &measStart;
	task.tick  = &measTick;
	task.stop  = &measStop;

	return task;
}

typedef struct measTaskState
{
	Sampler sampler;
	void *samplerState;
	long sample;
	enum {RELAXING, SAMPLING} measStatus;
	MeasurementConf measConf;
	double intervalTime; /* Time since last sample (or start). */
} MeasTaskState;

void *measStart(void *initialData)
{
	Measurement *meas = (Measurement*) initialData;
	assert(meas != NULL);
	Sampler *sampler = &meas->sampler;

	MeasTaskState *state = malloc(sizeof(*state));

	state->intervalTime = 0;
	state->sample = 0;
	state->sampler = meas->sampler; /* struct copy */
	state->samplerState = samplerStart(sampler);
	state->measConf = meas->measConf; /* struct copy */
	state->measStatus = (meas->measConf.measureWait > 0 ?
				RELAXING : SAMPLING);

	free(initialData);
	return state;
}

bool measTick(void *state)
{
	MeasTaskState *measState = (MeasTaskState*) state;
	assert(measState != NULL);
	MeasurementConf *measConf = &measState->measConf;
	Sampler *sampler = &measState->sampler;
	double measWait = measConf->measureWait;
	double measInterval = measConf->measureInterval;
	double time = getTime();
	bool ret = true;

	if (measInterval < 0)
		return ret;

	switch (measState->measStatus) {
	case RELAXING:
		if (fmod(time, measWait / 100) < config.timeStep) {
			printf("\rRelax time %13f of %f",
					(time + measWait/100) / TIME_FACTOR, 
					measWait / TIME_FACTOR);
			fflush(stdout);
		}
		if (time >= measWait) {
			measState->intervalTime = time - measWait;
			printf("\nStarting measurement.\n");
			measState->measStatus = SAMPLING;
		}
		break;
	case SAMPLING:
		measState->intervalTime += config.timeStep; //TODO nicer?

		if (measState->intervalTime < measInterval)
			break;

		measState->intervalTime -= measInterval;
		ret = samplerSample(sampler, measState->sample, 
				measState->samplerState);
		measState->sample++;
		
		if (measConf->measureSamples < 0)
			break; /* Go on indefinitely, don't print anything */

		printf("\rMeasured sample %ld/%ld", measState->sample,
						    measConf->measureSamples);
		fflush(stdout);
		if (measState->sample >= measConf->measureSamples) {
			printf("\nFinished sampling!\n");
			ret = false; /* Stop running */
		}
		break;
	default:
		fprintf(stderr, "Unknown measurement status!\n");
		assert(false);
	}
	return ret;
}

void measStop(void *state)
{
	MeasTaskState *measState = (MeasTaskState*) state;
	Sampler *sampler = &measState->sampler;

	sampler->stop(measState->sample, measState->samplerState);

	free(measState);
}




/* Some wrappers for samplers: */

/* Returns the state pointer that gets returned from sampler.start(), or 
 * NULL if sampler.start == NULL */
void *samplerStart(Sampler *sampler)
{
	if (sampler->start == NULL)
		return NULL;
	return sampler->start(sampler->samplerConf);
}

/* Sample and return sampler.sample(), or do nothing if sampler.sample == 
 * NULL and return true. Note that this would be a pretty useless sampler 
 * in the latter case... */
bool samplerSample(Sampler *sampler, long i, void *state)
{
	if (sampler->sample == NULL)
		return true;
	return sampler->sample(i, state);
}

/* Stop the sampler, or do nothing if sampler.stop == NULL */
void samplerStop(Sampler *sampler, long n, void *state)
{
	if (sampler->stop == NULL)
		return;
	sampler->stop(n, state);
}



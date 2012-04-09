#ifndef _MEASURE_H_
#define _MEASURE_H_

#include "physics.h"
#include <stdlib.h>

/* Configuration of a generic measurement */
typedef struct measurementConf
{
	/* Total number of samples to accumulate. Negative to go on 
	 * indefinitely. */
	long measureSamples;

	/* Time between samples. Negative to disable measurement. */
	double measureInterval; 

	/* Time to wait before starting measurement */
	double measureWait;
} MeasurementConf;

/* A Sampler does the actual measurement by sampling the state of the 
 * world. */
typedef struct sampler
{
	/* Data pointer that gets passed to the start() function below. Can 
	 * be used to pass configuration data. Make sure that this gets 
	 * allocated on the heap, so it will stay resident when we start 
	 * the sampling! (It can be freed in start() and should at least 
	 * get freed in stop()) */
	void *samplerConf;

	/* Called at the start of the measurement. Returns the state 
	 * pointer that gets passed to the functions below. */
	void *(*start)(void *samplerConf);

	/* Called after sample step i (starts at 0). Gets passed the 
	 * data pointer that start() returned. Returns false if the 
	 * measurement has to be stopped, true if everything can continue. */
	bool (*sample)(long i, void *state);

	/* Called at the end of the measurement. The total number of 
	 * samples that were actually measured is given in n.*/
	void (*stop)  (long n, void *state);
} Sampler;

/* Everything we need to know about a measurement. */
typedef struct measurement
{
	MeasurementConf measConf;
	Sampler sampler;
} Measurement;

/* Generate a task that performs the given measurement. */
Task measurementTask(Measurement *measurement);

#endif

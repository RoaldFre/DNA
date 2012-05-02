#ifndef _MEASURE_H_
#define _MEASURE_H_

#include "physics.h"
#include <stdlib.h>

/* Configuration of a generic measurement */
typedef struct {
	/* Total number of samples to accumulate. Negative to go on 
	 * indefinitely. */
	long measureSamples;

	/* Time between samples. Negative to disable measurement. */
	double measureInterval; 

	/* Time to wait before starting measurement */
	double measureWait;

	/* Path to a file to dump the measurement in. NULL means dump to 
	 * stdout. */
	const char *measureFile;

	/* Size of the string buffer to allocate for rendering the output 
	 * of the sampler to screen. Set to 0 or less to disable rendering 
	 * of the sampler output. */
	int renderStrBufSize;

	/* Pixel coordinates to render the above string at */
	int x, y;
} MeasurementConf;

typedef struct {
	/* The current sample, or the total number of samples performed if 
	 * the sampling is stopped. */
	long sample;

	/* String that can be written to, and the result will be rendered 
	 * if rendering is enabled. */
	char *string;

	/* Size of the string buffer above. */
	int strBufSize;
} SamplerData;

/* A Sampler does the actual measurement by sampling the state of the 
 * world. */
typedef struct {
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
	bool (*sample)(SamplerData *sd, void *state);

	/* Called at the end of the measurement. The total number of 
	 * samples that were actually measured is given in n.*/
	void (*stop)  (SamplerData *sd, void *state);
} Sampler;

/* Everything we need to know about a measurement. */
typedef struct {
	MeasurementConf measConf;
	Sampler sampler;
} Measurement;

/* Generate a task that performs the given measurement. If the measurement 
 * config has a valid renderStrBufSize, this task wil also render the 
 * string created by the sampler. */
Task measurementTask(Measurement *measurement);

#endif

#ifndef _SYSTEM_H_
#define _SYSTEM_H_

#include <stdbool.h>
#include "vmath.h"

typedef struct
{
	double timeStep;        /* The timestep (dt) in the simulation */
	double thermostatTemp;  /* Thermostat temperature */
	double thermostatTau;   /* Thermostat relaxation time */
	double langevinGamma;	/* Friction coefficient for Langevin dynamics */
	double truncationLen;	/* Length at which potentials are truncated */
	double saltConcentration; /* Na+ concentration in the environment */
} Config;

extern Config config;



/* Tasks are used to do stuff during the simulation, such as (accumulating) 
 * measurements. */
typedef enum
{
	TASK_OK,	/* Task ticked correctly and does not want to stop. */
	TASK_STOP,	/* Regular stop request. Everything went according to plan. */
	TASK_ERROR,	/* Error stop request. Something went wrong! */
	/* Note: the ordering of these enums is such that lower enums (ie 
	 * higher numbers) get precedence over previous ones in the list 
	 * (ie lower numbers). */
} TaskSignal;
typedef struct
{
	/* Data pointer that gets passed to the start() function below. Can 
	 * be used to pass configuration data. Make sure that this gets 
	 * allocated on the heap, so it will stay resident when we start 
	 * the simulation! */
	void *initialData;

	/* Called at the start of the simulation run. Returns the state 
	 * pointer that gets passed to the functions below. */
	void *(*start)(void *initialData);

	/* Called at every iteration step. Gets passed the data pointer 
	 * that start() returned.
	 * Time gets updated after the call to tick(), ie the simulation 
	 * time when tick() gets called first is always 0. */
	TaskSignal (*tick)(void *state);

	/* Called at the end of the simulation run.  Gets passed the data 
	 * pointer that start() returned. */
	void (*stop)(void *state);
} Task;

/* Run the given task in the simulation. Returns true if everything went 
 * according to plan. False if the task requested to stop because of an 
 * error. */
bool run(Task *task);

/* Get the time in the simulation after the last completed iteration. */
double getTime(void);

/* Get the number of the last completed iteration. */
long getIteration(void);

#endif

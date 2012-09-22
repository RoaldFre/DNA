#ifndef _INTEGRATOR_H_
#define _INTEGRATOR_H_

#include "system.h"

/* Integrator stuff */

typedef enum {
	LANGEVIN,
	VERLET
} IntegratorType;

typedef struct {
	real tau; 	/* Berendsen thermostat relaxation time */
} VerletSettings;

typedef struct {
	real gamma; 	/* Friction coefficient */
} LangevinSettings;

typedef struct {
	IntegratorType type;
	union {
		LangevinSettings langevin;
		VerletSettings verlet;
	} settings;
} Integrator;

typedef struct {
	Integrator integrator;
	real timeStep;      /* The timestep (dt) in the simulation */
	int numBoxes;         /* For space partition grid */
	real reboxInterval; /* Time interval for reboxing particles */
} IntegratorConf;

Task makeIntegratorTask(IntegratorConf *conf);

/* Returns a string with information about the integrator. Each line is 
 * prefixed with '#'.  You need to free the string pointer afterwards! */
char *integratorInfo(IntegratorConf *conf);

real getTimeStep(void);
void setTimeStep(real dt);

/* Get the simulated time. */
real getTime(void);



/* Heat bath stuff */

real getHeatBathTemperature(void);
void setHeatBathTemperature(real temperature);

typedef struct {
	/* Set the heat bath temperature to <temperature> at time <time>. */
	real time;
	real temperature;
} TemperatureSetpoint;

typedef struct {
	/* <setpoints> is a list that is SORTED ON INCREASING TIME and 
	 * consists in total of <numSetpoints> temperature setpoints. The 
	 * pointer must remain valid throughout the simulation run. */
	TemperatureSetpoint* setpoints;
	int numSetpoints;
} TemperatureTable;

Task makeTemperatureTask(TemperatureTable table);

#endif

#ifndef _INTEGRATOR_H_
#define _INTEGRATOR_H_

#include "system.h"

typedef enum {
	LANGEVIN,
	VERLET
} IntegratorType;

typedef struct {
	double tau; 	/* Berendsen thermostat relaxation time */
} VerletSettings;

typedef struct {
	double gamma; 	/* Friction coefficient */
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
	double timeStep;      /* The timestep (dt) in the simulation */
	int numBoxes;         /* For space partition grid */
	double reboxInterval; /* Time interval for reboxing particles */
} IntegratorConf;

Task makeIntegratorTask(IntegratorConf *conf);

/* Returns a string with information about the integrator. Each line is 
 * prefixed with '#'.  You need to free the string pointer afterwards! */
char *integratorInfo(IntegratorConf *conf);

double getTimeStep(void);
void setTimeStep(double dt);

/* Get the simulated time. */
double getTime(void);

double getHeatBathTemperature(void);
void setHeatBathTemperature(double temperature);

#endif

#define _GNU_SOURCE

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/resource.h>
#include "main.h"
#include "vmath.h"
#include "system.h"
#include "render.h"

static void parseArguments(int argc, char **argv);

/* Defaults */
#define DEF_TIMESTEP 			1.0
#define DEF_TEMPERATURE 		300.0
#define DEF_COUPLING_TIMESTEP_FACTOR 	1000
#define DEF_MONOMER_WORLDSIZE_FACTOR    5.5
#define DEF_MONOMERS_PER_RENDER 	2000
#define DEF_RENDER_RADIUS 		1.5

static void printUsage(void)
{
	printf("Usage: main <number of monomers> [flags]\n");
	printf("Flags:\n");
	printf(" -t <flt>  length of Time steps (in femtoseconds)\n");
	printf("             default: %f\n", DEF_TIMESTEP);
	printf(" -T <flt>  Temperature\n");
	printf("             default: %f\n", DEF_TEMPERATURE);
	printf(" -c <flt>  thermal bath Coupling: relaxation time (zero to disable)\n");
	printf("             default: %d * timestep\n", DEF_COUPLING_TIMESTEP_FACTOR);
	printf(" -j <int>  perform <int> physics steps between rendering frames.\n");
	printf("             default: %d/(number of monomers)\n", DEF_MONOMERS_PER_RENDER);
	printf(" -r        Render\n");
	printf(" -R <flt>  Radius (in Angstrom) of the particles when rendering\n");
	printf("             default: %f\n", DEF_RENDER_RADIUS);
	printf(" -S <flt>  Size of world (in Angstrom) (for rendering).\n");
	printf("             default: (number of monomers) * %f\n", DEF_MONOMER_WORLDSIZE_FACTOR);
	printf(" -v <int>  Verbose: dump statistics every <int> iterations\n");
	printf(" -E <flt>  dump Energy statistics every <flt> femtoseconds\n");
	printf(" -s <int>  accumulate <int> measurement Samples\n");
	printf("             default: Don't sample, loop forever\n");
}

static void parseArguments(int argc, char **argv)
{
	int c;

	/* defaults */
	config.verbose          = 0;
	config.measureSamples   = -1; /* loop indefinitely */
	config.timeStep 	= DEF_TIMESTEP * TIME_FACTOR;
	config.thermostatTemp	= DEF_TEMPERATURE;
	config.radius		= DEF_RENDER_RADIUS * LENGTH_FACTOR;

	/* guards */
	config.renderSteps = -1;
	config.worldSize = -1;
	config.thermostatTau = -1;
	config.measureInterval = -1;

	while ((c = getopt(argc, argv, ":t:T:c:j:rR:S:v:E:s:h")) != -1)
	{
		switch (c)
		{
		case 't':
			config.timeStep = atof(optarg) * TIME_FACTOR;
			if (config.timeStep <= 0)
				die("Invalid timestep %s\n", optarg);
			break;
		case 'T':
			config.thermostatTemp = atof(optarg);
			if (config.thermostatTemp < 0)
				die("Invalid temperature %s\n", optarg);
			break;
		case 'c':
			config.thermostatTau = atof(optarg);
			if (config.thermostatTau < 0)
				die("Invalid thermostat relaxation time %s\n",
						optarg);
			break;
		case 'j':
			config.renderSteps = atol(optarg);
			if (config.renderSteps < 0)
				die("Invalid number of renderer steps %s\n",
						optarg);
			break;
		case 'r':
			config.render = true;
			break;
		case 'R':
			config.radius = atof(optarg) * LENGTH_FACTOR;
			if (config.radius <= 0)
				die("Invalid radius %s\n", optarg);
			break;
		case 'S':
			config.worldSize = atof(optarg) * 1e-10;
			if (config.worldSize <= 0)
				die("Invalid world size %s\n", optarg);
			break;
		case 'v':
			config.verbose = atoi(optarg);
			if (config.verbose <= 0)
				die("Verbose: invalid number of iterations %s\n",
						optarg);
			break;
		case 'E':
			config.measureInterval = atoi(optarg) * TIME_FACTOR;
			if (config.measureInterval <= 0)
				die("Verbose: invalid measurement interval %s\n",
						optarg);
			break;
		case 's':
			config.measureSamples = atol(optarg);
			if (config.measureSamples < 0)
				die("Invalid number of samples %d\n",
						config.measureSamples);
			break;
		case 'h':
			printUsage();
			exit(0);
			break;
		case ':':
			printUsage();
			die("Option -%c requires an argument\n", optopt);
			break;
		case '?':
			printUsage();
			die("Option -%c not recognized\n", optopt);
			break;
		default:
			/* XXX */
			break;
		}
	}

	argc -= optind;
	argv += optind;

	if (argc < 1) {
		printUsage();
		die("\nNot enough required arguments!\n");
	}

	config.numMonomers = atoi(argv[0]);
	if (config.numMonomers <= 0)
		die("Invalid number of monomers!\n");

	if (config.worldSize < 0)
		config.worldSize = 1e-10 * config.numMonomers 
					* DEF_MONOMER_WORLDSIZE_FACTOR;
	if (config.thermostatTau < 0)
		config.thermostatTau = DEF_COUPLING_TIMESTEP_FACTOR
						* config.timeStep;
	if (config.renderSteps < 0)
		config.renderSteps = 1 + DEF_MONOMERS_PER_RENDER 
						/ config.numMonomers;

	if (config.measureInterval > 0 
			&& config.measureInterval < config.timeStep)
		die("Measurement interval %f smaller than timestep of %f!\n",
				config.measureInterval / TIME_FACTOR,
				config.timeStep / TIME_FACTOR);
	if (config.measureInterval > 0 
			&& config.measureInterval < 5 * config.timeStep)
		printf("WARNING: Measurement interval %f not much larger than "
				"timestep of %f! Exact timing will be off!\n",
				config.measureInterval / TIME_FACTOR,
				config.timeStep / TIME_FACTOR);
}

void die(const char *fmt, ...)
{
	va_list args;

	va_start(args, fmt);
	vfprintf(stderr, fmt, args);
	va_end(args);

	exit(1);
}

/* Advance the simulation by one time step. Render and/or dump statistics 
 * if neccesary. Return false if the user wants to quit. */
static bool stepSimulation(void) {
	static int stepsSinceRender = 0;
	static int stepsSinceVerbose = 0;
	static double timeSinceMeasurement = 0;
	static long samples = 0;

	stepWorld();
	assert(physicsCheck());

	if (config.verbose > 0) {
		stepsSinceVerbose++;
		if (stepsSinceVerbose > config.verbose) {
			stepsSinceVerbose = 0;
			dumpStats();
		}
	}

	if (config.render) {
		stepsSinceRender++;
		if (stepsSinceRender > config.renderSteps) {
			stepsSinceRender = 0;
			return stepGraphics();
		}
	}

	if (config.measureInterval > 0) {
		timeSinceMeasurement += config.timeStep;
		if (timeSinceMeasurement >= config.measureInterval) {
			timeSinceMeasurement -= config.measureInterval;
			dumpEnergies(stdout);

			samples++;
			if (samples >= config.measureSamples)
				return false;
		}
	}

	return true;
}

int main(int argc, char **argv)
{
	srand(2); //seed random generator

	parseArguments(argc, argv);

	allocWorld();
	fillWorld();
	
	if (config.render)
		initRender();

	while (stepSimulation());
}


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
#define DEF_TIMESTEP 			0.001
#define DEF_TEMPERATURE 		2.0
#define DEF_MONOMER_WORLDSIZE_FACTOR    8.0
#define DEF_PARTICLES_PER_RENDER 	10000
#define DEF_RENDER_RADIUS 		1.5

static void printUsage(void)
{
	printf("Usage: main <number of monomers> [flags]\n");
	printf("Flags:\n");
	printf(" -t <flt>  length of Time steps\n");
	printf("             default: %f\n", DEF_TIMESTEP);
	printf(" -T <flt>  Temperature\n");
	printf("             default: %f\n", DEF_TEMPERATURE);
	printf(" -j <int>  perform <int> physics steps between rendering frames.\n");
	printf("             default: %d/(number of monomers)\n", DEF_PARTICLES_PER_RENDER);
	printf(" -r        Render\n");
	printf(" -R <flt>  Radius (in Angstrom) of the particles when rendering\n");
	printf("             default: %f\n", DEF_RENDER_RADIUS);
	printf(" -S <flt>  Size of world (in Angstrom) (for rendering).\n");
	printf("             default: (number of monomers) * %f\n", DEF_MONOMER_WORLDSIZE_FACTOR);
	printf(" -v <int>  Verbose: dump statistics every <int> iterations\n");
}

static void parseArguments(int argc, char **argv)
{
	int c;

	/* defaults */
	config.verbose          = 0;
	config.timeStep 	= DEF_TIMESTEP;
	config.temperature	= DEF_TEMPERATURE;
	config.radius		= DEF_RENDER_RADIUS * 1e-10;

	/* guards */
	config.renderSteps = -1;
	config.worldSize = -1;

	while ((c = getopt(argc, argv, ":t:T:j:rR:S:v:h")) != -1)
	{
		switch (c)
		{
		case 't':
			config.timeStep = atof(optarg);
			if (config.timeStep <= 0)
				die("Invalid timestep %f\n",
						config.timeStep);
			break;
		case 'T':
			config.temperature = atof(optarg);
			if (config.temperature < 0)
				die("Invalid temperature %f\n",
						config.temperature);
			break;
		case 'j':
			config.renderSteps = atol(optarg);
			if (config.renderSteps < 0)
				die("Invalid number of renderer steps %d\n",
						config.renderSteps);
			break;
		case 'r':
			config.render = true;
			break;
		case 'R':
			config.radius = atof(optarg) * 1e-10;
			if (config.radius <= 0)
				die("Invalid radius %f\n", config.radius);
			break;
		case 'S':
			config.worldSize = atof(optarg) * 1e-10;
			if (config.worldSize <= 0)
				die("Invalid world size %f\n", config.worldSize);
			break;
		case 'v':
			config.verbose = atoi(optarg);
			if (config.verbose <= 0)
				die("Verbose: invalid number of iterations %d\n",
						config.verbose);
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

	if (config.renderSteps < 0)
		config.renderSteps = 1 + DEF_PARTICLES_PER_RENDER 
						/ config.numMonomers;
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

	stepWorld();

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

	return true;
}

int main(int argc, char **argv)
{
	bool keepGoing = true;

	parseArguments(argc, argv);

	allocWorld();
	fillWorld();
	
	if (config.render)
		initRender();

	while (stepSimulation());
}


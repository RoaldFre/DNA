#define _GNU_SOURCE

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <string.h>
#include "main.h"
#include "world.h"
#include "render.h"
#include "samplers.h"
#include "math.h"

#define DEF_DATA_PATH "data"

/* Defaults */
#define DEF_BASE_SEQUENCE		"ACCAATTTTTTTTTTTTTTGGG" /* T_12 in Bonnet */
#define DEF_TIMESTEP 			1.0
#define DEF_TEMPERATURE 		300.0
#define DEF_LANGEVIN_GAMMA		5e12 //TODO sane?
#define DEF_COUPLING_TIMESTEP_FACTOR 	1000
#define DEF_TRUNCATION_LENGTH		20.0 //TODO sane?
#define DEF_MONOMER_WORLDSIZE_FACTOR    80.0
#define DEF_MONOMERS_PER_RENDER 	2000
#define DEF_MEASUREMENT_WAIT 		4e4
#define DEF_RENDER_FRAMERATE 		30.0
#define DEF_RENDER_RADIUS 		0.5
#define DEF_INTEGRATOR			LANGEVIN


/* Static global configuration variables */

static MeasurementConf verboseConf =
{
	.measureTime = -1, /* loop forever */
	.measureInterval = -1, /* default: disable */
	.measureWait = 0,
	.measureFile = NULL,
	.verbose = false,
};
/* generic measurement config for whatever measurement we will do */
static MeasurementConf measurementConf =
{
	.measureTime = -1,
	.measureInterval = -1,
	.measureWait = 0,
	.measureFile = DEF_DATA_PATH,
	.verbose = true,
	.renderStrBufSize = 64,
	.renderStrX = 10,
	.renderStrY = 80,
};
static RenderConf renderConf =
{
	.framerate = DEF_RENDER_FRAMERATE,
	.radius    = DEF_RENDER_RADIUS * LENGTH_FACTOR,
	.drawForces = false,
};
static bool render;
static IntegratorConf integratorConf =
{
	.integrator = DEF_INTEGRATOR,
	.numBoxes   = -1, /* guard */
	.interactionSettings = {
			.enableBond	= true,
			.enableAngle	= true,
			.enableDihedral	= true,
			.enableStack	= true,
			.enableExclusion= true,
			.enableBasePair	= false,
			.enableCoulomb	= false,
			.mutuallyExclusivePairForces = true,
			.basePairInteraction = BASE_PAIR_ALL,
	},
};
static const char* baseSequence = DEF_BASE_SEQUENCE;
static bool buildCompStrand = false;
static double worldSizeFactor = DEF_MONOMER_WORLDSIZE_FACTOR;
static double worldSize = -1;


static void printUsage(void)
{
	printf("Usage: main [flags]\n");
	printf("\n");
	printf("Flags:\n");
	printf(" -s <str>  base Sequence of the DNA strand to simulate\n");
	printf("             default: %s\n", DEF_BASE_SEQUENCE);
	printf(" -d        build a Double helix with complementary strand as well\n");
	printf(" -t <flt>  length of Time steps (in femtoseconds)\n");
	printf("             default: %f\n", DEF_TIMESTEP);
	printf(" -T <flt>  Temperature\n");
	printf("             default: %f\n", DEF_TEMPERATURE);
	printf(" -r        Render\n");
	printf(" -f <flt>  desired Framerate when rendering.\n");
	printf("             default: %f)\n", DEF_RENDER_FRAMERATE);
	printf(" -R <flt>  Radius (in Angstrom) of the particles when rendering\n");
	printf("             default: %f\n", DEF_RENDER_RADIUS);
	printf(" -F        draw Forces on particles when rendering\n");
	printf(" -l <flt>  truncation Length of potentials (in Angstrom).\n");
	printf("             negative value: sets truncation to worldsize/2\n");
	printf("             default: %f\n", DEF_TRUNCATION_LENGTH);
	printf(" -S <flt>  worldSize factor per monomer (in Angstrom).\n");
	printf("             default: %f\n", DEF_MONOMER_WORLDSIZE_FACTOR);
	printf(" -b <num>  number of Boxes per dimension\n");
	printf("             default: max so that boxsize >= potential truncation length\n");
	printf(" -v <int>  Verbose: dump statistics every <flt> femtoseconds\n");
	printf(" -i <type> Integrator to use. Values for <type>:\n");
	printf("             l: Langevin (velocity BBK) [default]\n");
	printf("             v: velocity Verlet with Berendsen thermostat\n");
	printf("\n");
	printf("Parameters for measurements:\n");
	printf(" -W <flt>  Waiting time before starting the measurement (in nanosectonds)\n");
	printf("             default: 0 (ie, no relaxation phase)\n");
	printf(" -I <flt>  sample Interval (in picosectonds)\n");
	printf("             default: don't measure\n");
	printf(" -P <flt>  measurement Period: total time to sample the system (in nanosectonds)\n");
	printf("             default: sample indefinitely\n");
	printf(" -D <path> Data file to Dump measurement output\n");
	printf("             default: %s\n", DEF_DATA_PATH);
	printf("\n");
	printf("Parameters for Langevin integrator:\n");
	printf(" -g <flt>  Gamma: friction coefficient for Langevin dynamics\n");
	printf("             default: %e\n", DEF_LANGEVIN_GAMMA);
	printf("\n");
	printf("Parameters for velocity Verlet integrator + Berendsen termostat:\n");
	printf(" -c <flt>  thermal bath Coupling: relaxation time (zero to disable)\n");
	printf("             default: %d * timestep\n", DEF_COUPLING_TIMESTEP_FACTOR);
}

static void parseArguments(int argc, char **argv)
{
	int c;

	/* defaults */
	config.timeStep 	= DEF_TIMESTEP * TIME_FACTOR;
	config.thermostatTemp	= DEF_TEMPERATURE;
	config.truncationLen    = DEF_TRUNCATION_LENGTH * LENGTH_FACTOR;
	config.langevinGamma	= DEF_LANGEVIN_GAMMA;

	/* guards */
	config.thermostatTau = -1;

	while ((c = getopt(argc, argv, ":s:dt:T:g:c:f:rR:Fl:S:b:v:i:W:I:P:D:h")) != -1)
	{
		switch (c)
		{
		case 's':
			baseSequence = optarg;
			break;
		case 'd':
			buildCompStrand = true;
			break;
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
		case 'g':
			config.langevinGamma = atof(optarg);
			if (config.langevinGamma < 0)
				die("Invalid friction coefficient %s\n", optarg);
			break;
		case 'c':
			config.thermostatTau = atof(optarg) * TIME_FACTOR;
			if (config.thermostatTau < 0)
				die("Invalid thermostat relaxation time %s\n",
						optarg);
			break;
		case 'f':
			renderConf.framerate = atof(optarg);
			if (renderConf.framerate < 0)
				die("Invalid framerate %s\n", optarg);
			break;
		case 'r':
			render = true;
			break;
		case 'R':
			renderConf.radius = atof(optarg) * LENGTH_FACTOR;
			if (renderConf.radius <= 0)
				die("Invalid radius %s\n", optarg);
			break;
		case 'F':
			renderConf.drawForces = true;
			break;
		case 'l':
			config.truncationLen = atof(optarg) * LENGTH_FACTOR;
			break;
		case 'S':
			worldSizeFactor = atof(optarg) * A;
			if (worldSizeFactor <= 0)
				die("Invalid world size factor %s\n", optarg);
			break;
		case 'b':
			integratorConf.numBoxes = atoi(optarg);
			if (integratorConf.numBoxes <= 0)
				die("Invalid number of boxes %s\n",
						optarg);
			break;
		case 'v':
			verboseConf.measureInterval = atof(optarg) * TIME_FACTOR;
			if (verboseConf.measureInterval <= 0)
				die("Verbose: invalid verbose interval %s\n",
						optarg);
			break;
		case 'i':
			if (optarg[0] == '\0' || optarg[1] != '\0')
				die("Integrator: badly formatted integrator type\n");
			switch(optarg[0]) {
			case 'l': integratorConf.integrator = LANGEVIN; break;
			case 'v': integratorConf.integrator = VERLET; break;
			default: die("Unknown integrator type '%s'\n", optarg);
				 break;
			}
			break;
		case 'W':
			measurementConf.measureWait = atof(optarg) * NANOSECONDS;
			if (measurementConf.measureWait < 0)
				die("Invalid relaxation time %s\n", optarg);
			break;
		case 'I':
			measurementConf.measureInterval = atof(optarg) * PICOSECONDS;
			if (measurementConf.measureInterval<= 0)
				die("Invalid measurement interval %s\n", optarg);
			break;
		case 'P':
			measurementConf.measureTime = atof(optarg) * NANOSECONDS;
			if (measurementConf.measureTime <= 0)
				die("Invalid measurement time %s\n", optarg);
			break;
		case 'D':
			measurementConf.measureFile = optarg;
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
			die("Error parsing options!");
			break;
		}
	}

	argc -= optind;
	argv += optind;

	if (argc > 0) {
		printUsage();
		die("\nFound unrecognised option(s) at the command line!\n");
	}

	worldSize = strlen(baseSequence) * worldSizeFactor;

	if (config.thermostatTau < 0)
		config.thermostatTau = DEF_COUPLING_TIMESTEP_FACTOR
						* config.timeStep;

	if (config.truncationLen < 0) {
		printf("Disabling space partitioning\n");
		/* Disable truncation -> no space partitioning */
		integratorConf.numBoxes = 1;
		config.truncationLen = worldSize / 2.0;
		/* Due to (cubic) periodicity, we still need to truncate at 
		 * worldsize/2 to have correct energy conservation and to 
		 * make sure every particle 'sees' the same (spherical) 
		 * potential, no matter where it is within the cube. */
	} else if (config.truncationLen > worldSize / 2.0) {
		printf("Truncation (%e) > worldSize/2 (%e)\n"
				"   => Disabling space partitioning\n",
				config.truncationLen, worldSize/2);
		config.truncationLen = worldSize / 2.0; /* same reason */
	}

	if (integratorConf.numBoxes == -1) {
		integratorConf.numBoxes = worldSize / config.truncationLen;
		if (integratorConf.numBoxes  < 1)
			integratorConf.numBoxes = 1;
		printf("Number of boxes per dimension: %d\n",
				integratorConf.numBoxes);
	}

	if (worldSize / integratorConf.numBoxes < config.truncationLen)
		die("The boxsize (%e) is smaller than the potential "
			"truncation radius (%e)!\n",
			worldSize / integratorConf.numBoxes / LENGTH_FACTOR,
			config.truncationLen / LENGTH_FACTOR);
}

void die(const char *fmt, ...)
{
	va_list args;

	va_start(args, fmt);
	vfprintf(stderr, fmt, args);
	va_end(args);

	exit(1);
}

int main(int argc, char **argv)
{
	seedRandom();

	parseArguments(argc, argv);

	allocWorld(1, worldSize);
	fillStrand(&world.strands[0], baseSequence);

	killMomentum();

	assert(worldSanityCheck());

	Measurement verbose;
	verbose.measConf = verboseConf;
	verbose.sampler = dumpStatsSampler();
	Task verboseTask = measurementTask(&verbose);

	const char *filenameBase = measurementConf.measureFile;
	char *particleFile, *strandFile;
	if (0 >	asprintf(&particleFile, "%s_particle", filenameBase))
		die("Could not create particle file name string!\n");
	if (0 > asprintf(&strandFile, "%s_strand", filenameBase))
		die("Could not create strand file name string!\n");

	Measurement particleDiff;
	particleDiff.sampler = particleSquaredDisplacementSampler(
						&world.strands[0].Ps[0]);
	particleDiff.measConf = measurementConf;
	particleDiff.measConf.measureFile = particleFile;
	Task particleDiffTask = measurementTask(&particleDiff);


	Measurement strandDiff;
	strandDiff.sampler = strandCOMSquaredDisplacementSampler(
						&world.strands[0]);
	strandDiff.measConf = measurementConf;
	strandDiff.measConf.measureFile = strandFile;
	strandDiff.measConf.verbose = false; //particleDiff already does it
	Task strandDiffTask = measurementTask(&strandDiff);


	Task renderTask = makeRenderTask(&renderConf);

	Task integratorTask = makeIntegratorTask(&integratorConf);

	Task *tasks[5];
	tasks[0] = (render ? &renderTask : NULL);
	tasks[1] = &integratorTask;
	tasks[2] = &verboseTask;
	tasks[3] = &particleDiffTask;
	tasks[4] = &strandDiffTask;
	Task task = sequence(tasks, 5);

	bool everythingOK = run(&task);

	free(strandFile);
	free(particleFile);

	if (!everythingOK)
		return 1;
	return 0;
}

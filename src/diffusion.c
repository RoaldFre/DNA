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
#define END_TO_END_DIST_FILE_SUFFIX "_endToEnd"
#define BASE_PAIRING_FILE_SUFFIX "_basePairing"
#define TEMPERATURE_FILE_SUFFIX "_temperature"

/* Defaults */
#define DEF_BASE_SEQUENCE		"AAAAAAAAAAAA" /* 12xA */
#define DEF_TIMESTEP 			15.0
#define DEF_REBOX_INTERVAL		200.0
#define DEF_INITIAL_TEMPERATURE		"20C"
#define DEF_SALT_CONCENTRATION		50.0  /* mol/m^3 */
#define DEF_LANGEVIN_GAMMA		5e12
#define DEF_COUPLING_TIMESTEP_FACTOR 	1000
#define DEF_TRUNCATION_LENGTH		20.0
#define DEF_MONOMER_WORLDSIZE_FACTOR    4.3
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
	.reboxInterval = DEF_REBOX_INTERVAL * FEMTOSECONDS,
	.interactionSettings = {
			.enableBond	= true,
			.enableAngle	= true,
			.enableDihedral	= true,
			.enableStack	= true,
			.enableExclusion= true,
			.enableBasePair	= true,
			.enableCoulomb	= true,
			.mutuallyExclusivePairForces = true,
			.basePairInteraction = BASE_PAIR_HAIRPIN,
	},
};
static const char* baseSequence = DEF_BASE_SEQUENCE;
static double worldSize = -1; /* guard */


static void printUsage(void)
{
	printf("Usage: main [flags]\n");
	printf("\n");
	printf("Flags:\n");
	printf(" -s <str>  base Sequence of the DNA strand to simulate\n");
	printf("             default: %s\n", DEF_BASE_SEQUENCE);
	printf(" -t <flt>  length of Time steps (in femtoseconds)\n");
	printf("             default: %f\n", DEF_TIMESTEP);
	printf(" -T <flt><C|K>  initial Temperature (example: 20C or 300K)\n");
	printf("             note: some measurements set their own temperature when sampling!\n");
	printf("             default: %s\n", DEF_INITIAL_TEMPERATURE);
	printf(" -N <flt>  concentration of Na+ in the environment (in mol/m^3)\n");
	printf("             default: %f\n", DEF_SALT_CONCENTRATION);
	printf(" -r        Render\n");
	printf(" -f <flt>  desired Framerate when rendering.\n");
	printf("             default: %f)\n", DEF_RENDER_FRAMERATE);
	printf(" -R <flt>  Radius (in Angstrom) of the particles when rendering\n");
	printf("             default: %f\n", DEF_RENDER_RADIUS);
	printf(" -F        draw Forces on particles when rendering\n");
	printf(" -l <flt>  truncation Length of potentials (in Angstrom).\n");
	printf("             negative value: sets truncation to worldsize/2\n");
	printf("             default: %f\n", DEF_TRUNCATION_LENGTH);
	printf(" -S <flt>  Size of world (in Angstrom).\n");
	printf("             default: (number of monomers + 2) * %f\n", DEF_MONOMER_WORLDSIZE_FACTOR);
	printf(" -b <num>  number of Boxes per dimension\n");
	printf("             default: max so that boxsize >= potential truncation length\n");
	printf(" -x <flt>  time interval between reboXing (in femtoseconds)\n");
	printf("             default: %f\n", DEF_REBOX_INTERVAL);
	printf(" -v <int>  Verbose: dump statistics every <flt> picoseconds\n");
	printf(" -i <type> Integrator to use. Values for <type>:\n");
	printf("             l: Langevin (velocity BBK) [default]\n");
	printf("             v: velocity Verlet with Berendsen thermostat\n");
	printf("\n");
	printf("Parameters for Langevin integrator:\n");
	printf(" -g <flt>  Gamma: friction coefficient for Langevin dynamics\n");
	printf("             default: %e\n", DEF_LANGEVIN_GAMMA);
	printf("\n");
	printf("Parameters for velocity Verlet integrator + Berendsen termostat:\n");
	printf(" -c <flt>  thermal bath Coupling: relaxation time (zero to disable)\n");
	printf("             default: %d * timestep\n", DEF_COUPLING_TIMESTEP_FACTOR);
	printf("\n");
	printf("\n");
	printf("Parameters for measurements:\n");
	printf(" -W <flt>  Waiting time before starting the measurement (in nanosectonds)\n");
	printf("             default: 0 (ie, no relaxation phase)\n");
	printf(" -I <flt>  sample Interval (in picosectonds)\n");
	printf("             default: don't measure\n");
	printf(" -P <flt>  measurement Period: total time to sample the system (in nanosectonds)\n");
	printf("             default: sample indefinitely\n");
	printf(" -D <path> Data file to Dump measurement output. The directory must exist.\n");
	printf("             default: %s\n", DEF_DATA_PATH);
}

static void parseArguments(int argc, char **argv)
{
	int c;

	/* defaults */
	config.timeStep 	 = DEF_TIMESTEP * FEMTOSECONDS;
	config.thermostatTemp	 = parseTemperature(DEF_INITIAL_TEMPERATURE);
	config.saltConcentration = DEF_SALT_CONCENTRATION;
	config.truncationLen     = DEF_TRUNCATION_LENGTH * LENGTH_FACTOR;
	config.langevinGamma	 = DEF_LANGEVIN_GAMMA;

	/* guards */
	config.thermostatTau = -1;

	while ((c = getopt(argc, argv, ":s:t:T:N:g:c:f:rR:Fl:S:b:x:v:i:W:I:P:D:")) != -1)
	{
		switch (c)
		{
		case 's':
			baseSequence = optarg;
			break;
		case 't':
			config.timeStep = atof(optarg) * FEMTOSECONDS;
			if (config.timeStep <= 0)
				die("Invalid timestep %s\n", optarg);
			break;
		case 'T':
			config.thermostatTemp = parseTemperature(optarg);
			if (config.thermostatTemp < 0)
				die("Invalid temperature!\n");
			break;
		case 'N':
			config.saltConcentration = atof(optarg);
			if (config.saltConcentration < 0)
				die("Invalid salt concentration %s\n", optarg);
			break;
		case 'g':
			config.langevinGamma = atof(optarg);
			if (config.langevinGamma < 0)
				die("Invalid friction coefficient %s\n", optarg);
			break;
		case 'c':
			config.thermostatTau = atof(optarg) * FEMTOSECONDS;
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
			worldSize = atof(optarg) * A;
			if (worldSize <= 0)
				die("Invalid world size %s\n", optarg);
			break;
		case 'b':
			integratorConf.numBoxes = atoi(optarg);
			if (integratorConf.numBoxes <= 0)
				die("Invalid number of boxes %s\n",
						optarg);
			break;
		case 'x':
			integratorConf.reboxInterval = atof(optarg) * FEMTOSECONDS;
			if (integratorConf.reboxInterval <= 0)
				die("Invalid rebox interval %s\n", optarg);
			break;
		case 'v':
			verboseConf.measureInterval = atof(optarg) * NANOSECONDS;
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

	if (worldSize < 0)
		worldSize = LENGTH_FACTOR * (strlen(baseSequence) + 2)
					* DEF_MONOMER_WORLDSIZE_FACTOR;

	if (config.thermostatTau < 0)
		config.thermostatTau = DEF_COUPLING_TIMESTEP_FACTOR
						* config.timeStep;

	if (config.truncationLen < 0) {
		/* Disable truncation -> no space partitioning */
		printf("Disabling space partitioning, "
				"maximizing truncation length.\n");
		integratorConf.numBoxes = 1;
		config.truncationLen = worldSize / 2.0;
		/* Due to (cubic) periodicity, we still need to truncate at 
		 * worldsize/2 to have correct energy conservation and to 
		 * make sure every particle 'sees' the same (spherical) 
		 * potential, no matter where it is within the cube. */
	} else if (config.truncationLen > worldSize / 2.0) {
		/* World is too small, extend it so we correctly compute 
		 * the potentials and forces up to the requested truncation 
		 * length. We need twice the truncation length for the same 
		 * reason as above. */
		printf("Truncation (%e) > worldSize/2 (%e)\n   => "
				"Extending worldSize to 2*Truncation (%e).\n",
				config.truncationLen, worldSize/2.0,
				2.0 * config.truncationLen);
		worldSize = 2.0 * config.truncationLen;
		integratorConf.numBoxes = 1; /* No space partitioning */
	} else if (integratorConf.numBoxes == -1) {
		/* Automatically determine ideal number of boxes.
		 * Formula is fitted for N between 100 and 1000 on a Core2 
		 * Duo E8500.
		 * (+ 0.5 for correct rounding to int) */
		int ideal = 0.5 + 2.73 * pow(strlen(baseSequence), 0.446);
		integratorConf.numBoxes = MIN(ideal,
					worldSize / config.truncationLen);
		if (integratorConf.numBoxes < 1)
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

	/* Verbose task */
	Measurement verbose;
	verbose.measConf = verboseConf;
	verbose.sampler = dumpStatsSampler();
	Task verboseTask = measurementTask(&verbose);

	/* Diffusion task */
	Measurement diffusion;
	diffusion.sampler = strandCOMSquaredDisplacementSampler(&world.strands[0]);
	diffusion.measConf = measurementConf;
	Task diffusionTask = measurementTask(&diffusion);

	/* Render task */
	Task renderTask = makeRenderTask(&renderConf);

	/* Integrator task */
	Task integratorTask = makeIntegratorTask(&integratorConf);

	/* Combined task */
	Task *tasks[4];
	tasks[0] = (render ? &renderTask : NULL);
	tasks[1] = &integratorTask;
	tasks[2] = &verboseTask;
	tasks[3] = &diffusionTask;;
	Task task = sequence(tasks, 4);

	bool everythingOK = run(&task);

	if (!everythingOK)
		return 1;

	return 0;
}


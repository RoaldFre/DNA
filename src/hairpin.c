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
#define DEF_BASE_SEQUENCE		"GCCTATTTTTTAATAGGC" /* N=4 in Kuznetsov nov 2001 */
#define DEF_TIMESTEP 			20.0
#define DEF_REBOX_INTERVAL		200.0
#define DEF_INITIAL_TEMPERATURE		70.0
#define DEF_SAMPLING_TEMPERATURE	20.0
#define DEF_SALT_CONCENTRATION		100.0  /* mol/m^3 */
#define DEF_LANGEVIN_GAMMA		5e12
#define DEF_COUPLING_TIMESTEP_FACTOR 	1000
#define DEF_TRUNCATION_LENGTH		20.0
#define DEF_MONOMER_WORLDSIZE_FACTOR    4.5
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
static BasePairingConfig bpc =
{
	.energyThreshold = -0.1 * EPSILON,
	.T = CELSIUS(DEF_SAMPLING_TEMPERATURE),
};
static HairpinFormationSamplerConfig hfc =
{
	.energyThreshold = -0.1 * EPSILON,
	.confirmationTime = 1 * NANOSECONDS,
	.allowedUnboundBPs = 1,
	.allowedBoundBPs = 1,
	.zippingTemperature = CELSIUS(20),
	.unzippingTemperature = CELSIUS(90),
	.zippedRelaxationTime = 5 * NANOSECONDS,
};
static HairpinMeltingTempSamplerConfig hmtc =
{
	.energyThreshold = -0.1 * EPSILON,
	.Tstart = 20,
	.Tstep = 10,
	.numSteps = 10,
	.relaxationTime = 5 * NANOSECONDS,
	.measureTime = 100 * NANOSECONDS,
	.verbose = false,
};
static enum {
	HAIRPIN_MELTING_TEMPERATURE,
	HAIRPIN_FORMATION_TIME,
	HAIRPIN_NO_MEASUREMENT,
} hairpinMeasurementType = HAIRPIN_NO_MEASUREMENT;

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
static bool buildCompStrand = false;
static double worldSize = -1; /* guard */


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
	printf(" -E <flt>  initial Equilibration temperature during relaxation phase (Celsius)\n");
	printf("             default: %f\n", DEF_INITIAL_TEMPERATURE);
	printf(" -T <flt>  Temperature during sampling (Celsius)\n");
	printf("             default: %f\n", DEF_SAMPLING_TEMPERATURE);
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
	printf("             default: (number of monomers) * %f\n", DEF_MONOMER_WORLDSIZE_FACTOR);
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
	printf(" -D <path> Data file to Dump measurement output\n");
	printf("             default: %s\n", DEF_DATA_PATH);
	printf(" -X <m|f>  measurement to perform:\n");
	printf("             m: hairpin Melting temperature\n");
	printf("             f: hairpin Formation time\n");
	printf("\n");
	printf("For more info and default values of params below: see code :P\n");
	//TODO parameter names are a mess
	printf("Parameters for hairpin melting temperature measurement:\n");
	printf(" -A <flt>  startTemp (Celsius)\n");
	printf(" -B <flt>  stepTemp\n");
	printf(" -C <int>  nSteps\n");
	printf(" -G <flt>  measureTime per step (nanoseconds)\n");
	printf(" -L <flt>  reLaxTime per step (nanoseconds)\n");
	printf(" -V        Be verbose: also dump hairpin state to file\n");
	printf("Parameters for hairpin formation measurement:\n");
	printf(" -H <int>  allowed unbounded base pairs\n");
	printf(" -M <int>  allowed bounded base pairs\n");
	printf(" -O <flt>  zipping temperature (Celsius)\n");
	printf(" -Q <flt>  unzipping temperature (Celsius)\n");
	printf(" -U <flt>  zipped relaxation phase duration (nanoseconds)\n");
	printf("\n");
}

static void parseArguments(int argc, char **argv)
{
	int c;

	/* defaults */
	config.timeStep 	 = DEF_TIMESTEP * FEMTOSECONDS;
	config.thermostatTemp	 = CELSIUS(DEF_INITIAL_TEMPERATURE);
	config.saltConcentration = DEF_SALT_CONCENTRATION;
	config.truncationLen     = DEF_TRUNCATION_LENGTH * LENGTH_FACTOR;
	config.langevinGamma	 = DEF_LANGEVIN_GAMMA;

	/* guards */
	config.thermostatTau = -1;

	while ((c = getopt(argc, argv, ":s:dt:E:T:N:g:c:f:rR:Fl:S:b:x:v:i:W:I:P:D:X:hA:B:C:G:L:VH:M:O:Q:U:")) != -1)
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
			config.timeStep = atof(optarg) * FEMTOSECONDS;
			if (config.timeStep <= 0)
				die("Invalid timestep %s\n", optarg);
			break;
		case 'E':
			config.thermostatTemp = CELSIUS(atof(optarg));
			if (config.thermostatTemp < 0)
				die("Invalid equilibration temperature %s\n", optarg);
			break;
		case 'T':
			bpc.T = CELSIUS(atof(optarg));
			if (bpc.T < 0)
				die("Invalid sampling temperature %s\n", optarg);
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
		/* A,B,C,G,H: because I'm running out of sensible letters ... */
		case 'A':
			hmtc.Tstart = CELSIUS(atof(optarg));
			break;
		case 'B':
			hmtc.Tstep = atof(optarg);
			break;
		case 'C':
			hmtc.numSteps = atoi(optarg);
			break;
		case 'G':
			hmtc.measureTime = atof(optarg) * NANOSECONDS;
			break;
		case 'L':
			hmtc.relaxationTime = atof(optarg) * NANOSECONDS;
			break;
		case 'V':
			hmtc.verbose = true;
			break;
		case 'H':
			hfc.allowedUnboundBPs = atoi(optarg);
			break;
		case 'M':
			hfc.allowedBoundBPs = atoi(optarg);
			break;
		case 'O':
			hfc.zippingTemperature = CELSIUS(atof(optarg));
			break;
		case 'Q':
			hfc.unzippingTemperature = CELSIUS(atof(optarg));
			break;
		case 'U':
			hfc.zippedRelaxationTime = atof(optarg) * NANOSECONDS;
			break;
		case 'X':
			switch(optarg[0]) {
			case 'f': 
				hairpinMeasurementType = HAIRPIN_FORMATION_TIME;
				break;
			case 'm':
				hairpinMeasurementType = HAIRPIN_MELTING_TEMPERATURE;
				break;
			default:
				die("Unknown hairpin measurement type '%s'\n", optarg);
				break;
			}
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
		worldSize = LENGTH_FACTOR * strlen(baseSequence)
					* DEF_MONOMER_WORLDSIZE_FACTOR;

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

	if (buildCompStrand)
		allocWorld(2, worldSize);
	else
		allocWorld(1, worldSize);

	fillStrand(&world.strands[0], baseSequence);
	if (buildCompStrand)
		fillComplementaryStrand(&world.strands[1], baseSequence);

	killMomentum();

	assert(worldSanityCheck());

	Measurement verbose;
	verbose.measConf = verboseConf;
	verbose.sampler = dumpStatsSampler();
	Task verboseTask = measurementTask(&verbose);

	const char *filenameBase = measurementConf.measureFile;
	char *basePairFile;
	if (0 >	asprintf(&basePairFile, "%s_basePair", filenameBase))
		die("Could not create base pair file name string!\n");

#if 0
	Measurement basePairing;
	basePairing.sampler = basePairingSampler(&bpc);
	basePairing.measConf = measurementConf;
	basePairing.measConf.measureFile = basePairFile;
	basePairing.measConf.verbose = false; /* Let output come from 
						 hairpin sampler */
	Task basePairingTask = measurementTask(&basePairing);
#endif

	Measurement hairpin;
	switch (hairpinMeasurementType) {
	case HAIRPIN_MELTING_TEMPERATURE:
		hairpin.sampler = hairpinMeltingTempSampler(&hmtc);
		break;
	case HAIRPIN_FORMATION_TIME:
		hairpin.sampler = hairpinFormationSampler(&hfc);
		break;
	case HAIRPIN_NO_MEASUREMENT:
		/* This is probably an error, so die */
		die("You should give a hairpin measurement type!\n");
	default:
		assert(false); die("Unknown hairpin measurement type!\n");
	}
	hairpin.measConf = measurementConf;
	Task hairpinTask = measurementTask(&hairpin);

	Task renderTask = makeRenderTask(&renderConf);

	Task integratorTask = makeIntegratorTask(&integratorConf);

	Task *tasks[5];
	tasks[0] = (render ? &renderTask : NULL);
	tasks[1] = &integratorTask;
	tasks[2] = &verboseTask;
	//tasks[3] = &basePairingTask;
	tasks[3] = NULL; /* Disabled atm */
	tasks[4] = &hairpinTask;;
	Task task = sequence(tasks, 5);
	bool everythingOK = run(&task);

	free(basePairFile);

	if (!everythingOK)
		return 1;
	return 0;
}


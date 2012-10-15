#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <string.h>
#include "system.h"
#include "world.h"
#include "integrator.h"
#include "render.h"
#include "samplers.h"
#include "math.h"

#define DEF_DATA_PATH "hairpinData"
#define END_TO_END_DIST_FILE_SUFFIX "_endToEnd"
#define BASE_PAIRING_FILE_SUFFIX "_basePairing"
#define TEMPERATURE_FILE_SUFFIX "_temperature"

/* Defaults */
#define DEF_BASE_SEQUENCE		"GCCTATTTTTTAATAGGC" /* N=4 in Kuznetsov nov 2001 */
#define DEF_TIMESTEP 			15.0
#define DEF_REBOX_INTERVAL		100.0
#define DEF_INITIAL_TEMPERATURE		"20C"
#define DEF_SALT_CONCENTRATION		100.0  /* mol/m^3 */
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
	.renderStrBufSize = 0,
	.renderStrX = 10,
	.renderStrY = 80,
};
static BasePairingConfig bpc =
{
	.energyThreshold = -0.1 * EPSILON,
};
static HairpinFormationSamplerConfig hfc =
{
	.energyThreshold = -0.1 * EPSILON,
	.requiredBoundBPs = -1, /* guard */
	.zipConfirmationTime = 0,
	.allowedBoundBPs = 0,
	.unzipConfirmationTime = 1 * NANOSECONDS,
	.zippingTemperature = CELSIUS(20),
	.unzippingTemperature = CELSIUS(90),
	.zippedRelaxationTime = 5 * NANOSECONDS,
};
static HairpinMeltingTempSamplerConfig hmtc =
{
	.energyThreshold = -0.1 * EPSILON,
	.Tstart = CELSIUS(20),
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

static bool measureEndToEndDistance = false;
static bool measureBasePairing = false;
static bool measureTemperature = false;

static RenderConf renderConf =
{
	.framerate = DEF_RENDER_FRAMERATE,
	.radius    = DEF_RENDER_RADIUS * ANGSTROM,
	.drawForces = false,
};
static bool render;

static IntegratorType integratorType = DEF_INTEGRATOR;
static VerletSettings verletSettings = {
	.tau = -1, /* guard */
};
static LangevinSettings langevinSettings = {
	.gamma = DEF_LANGEVIN_GAMMA,
};
static IntegratorConf integratorConf = {
	.timeStep      = DEF_TIMESTEP * FEMTOSECONDS,
	.numBoxes      = -1, /* guard */
	.reboxInterval = DEF_REBOX_INTERVAL * FEMTOSECONDS,
};

static double temperature;

static InteractionSettings interactionSettings = {
	.enableBond	= true,
	.enableAngle	= true,
	.enableDihedral	= true,
	.enableStack	= true,
	.enableExclusion= true,
	.enableBasePair	= true,
	.enableCoulomb	= true,
	.mutuallyExclusivePairForces = true,
	.basePairInteraction = BASE_PAIR_HAIRPIN,
	.saltConcentration   = DEF_SALT_CONCENTRATION,
	.truncationLen       = DEF_TRUNCATION_LENGTH * ANGSTROM,
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
	printf(" -W <flt>  Waiting time before starting the measurement (in nanoseconds)\n");
	printf("             default: 0 (ie, no relaxation phase)\n");
	printf(" -I <flt>  sample Interval (in picosectonds)\n");
	printf("             default: don't measure\n");
	printf(" -P <flt>  measurement Period: total time to sample the system (in nanoseconds)\n");
	printf("             default: sample indefinitely\n");
	printf(" -D <path> Data file to Dump measurement output. The directory must exist.\n");
	printf("             default: %s\n", DEF_DATA_PATH);
	printf(" -X <m|f>  measurement to eXecute:\n");
	printf("             m: hairpin Melting temperature\n");
	printf("             f: hairpin Formation time\n");
	printf(" -e        also measure End-to-end distance of the strand\n");
	printf("             output: the data filename (see -D) with suffix: '%s'\n",
							END_TO_END_DIST_FILE_SUFFIX);
	printf(" -p        also measure raw hairpin base Pairing state\n");
	printf("             output: the data filename (see -D) with suffix: '%s'\n",
							BASE_PAIRING_FILE_SUFFIX);
	printf(" -k        also measure the Kinetic temperature\n");
	printf("             output: the data filename (see -D) with suffix: '%s'\n",
							TEMPERATURE_FILE_SUFFIX);
	printf("\n");
	printf("For more info and default values of params below: see code :P\n");
	//TODO parameter names are a mess
	printf("Parameters for hairpin melting temperature measurement:\n");
	printf(" -A <flt><C|K> startTemp\n");
	printf("             default: %fC\n", TO_CELSIUS(hmtc.Tstart));
	printf(" -B <flt>  stepTemp\n");
	printf("             default: %f\n", hmtc.Tstep);
	printf(" -C <int>  nSteps\n");
	printf("             default: %d\n", hmtc.numSteps);
	printf(" -G <flt>  measureTime per step (nanoseconds)\n");
	printf("             default: %f\n", hmtc.measureTime / NANOSECONDS);
	printf(" -L <flt>  reLaxTime per step (nanoseconds)\n");
	printf("             default: %f\n", hmtc.relaxationTime / NANOSECONDS);
	printf(" -V        Be verbose: also dump hairpin state to file\n");
	printf("Parameters for hairpin formation measurement:\n");
	printf(" -H <int>  min required bounded base pairs for confirming 'zipped' state\n");
	printf(" -M <int>  max allowed bounded base pairs for confirming 'unzipped' state\n");
	printf("             default: %d\n", hfc.allowedBoundBPs);
	printf(" -O <flt><C|K> zipping temperature\n");
	printf("             default: %fC\n", TO_CELSIUS(hfc.zippingTemperature));
	printf(" -Q <flt><C|K> unzipping temperature\n");
	printf("             default: %fC\n", TO_CELSIUS(hfc.unzippingTemperature));
	printf(" -U <flt>  zipped relaxation phase duration (nanoseconds)\n");
	printf("             default: %f\n", hfc.zippedRelaxationTime / NANOSECONDS);
	printf("\n");
}

static void parseArguments(int argc, char **argv)
{
	int c;

	/* defaults */
	temperature = parseTemperature(DEF_INITIAL_TEMPERATURE);

	while ((c = getopt(argc, argv, ":s:dt:T:N:g:c:f:rR:Fl:S:b:x:v:i:W:I:P:D:X:epkhA:B:C:G:L:VH:M:O:Q:U:")) != -1)
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
			integratorConf.timeStep = atof(optarg) * FEMTOSECONDS;
			if (integratorConf.timeStep <= 0)
				die("Invalid timestep %s\n", optarg);
			break;
		case 'T':
			temperature = parseTemperature(optarg);
			if (temperature < 0)
				die("Invalid temperature!\n");
			break;
		case 'N':
			interactionSettings.saltConcentration = atof(optarg);
			if (interactionSettings.saltConcentration < 0)
				die("Invalid salt concentration %s\n", optarg);
			break;
		case 'g':
			langevinSettings.gamma = atof(optarg);
			if (langevinSettings.gamma < 0)
				die("Invalid friction coefficient %s\n", optarg);
			break;
		case 'c':
			verletSettings.tau = atof(optarg) * FEMTOSECONDS;
			if (verletSettings.tau < 0)
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
			renderConf.radius = atof(optarg) * ANGSTROM;
			if (renderConf.radius <= 0)
				die("Invalid radius %s\n", optarg);
			break;
		case 'F':
			renderConf.drawForces = true;
			break;
		case 'l':
			interactionSettings.truncationLen = atof(optarg) * ANGSTROM;
			break;
		case 'S':
			worldSize = atof(optarg) * ANGSTROM;
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
			verboseConf.measureInterval = atof(optarg) * PICOSECONDS;
			if (verboseConf.measureInterval <= 0)
				die("Verbose: invalid verbose interval %s\n",
						optarg);
			break;
		case 'i':
			if (optarg[0] == '\0' || optarg[1] != '\0')
				die("Integrator: badly formatted integrator type\n");
			switch(optarg[0]) {
			case 'l': integratorType = LANGEVIN; break;
			case 'v': integratorType = VERLET; break;
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
			hmtc.Tstart = parseTemperature(optarg);
			if (hmtc.Tstart < 0)
				die("Invalid starting temperature '%s'\n", optarg);
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
			hfc.requiredBoundBPs = atoi(optarg);
			break;
		case 'M':
			hfc.allowedBoundBPs = atoi(optarg);
			break;
		case 'O':
			hfc.zippingTemperature = parseTemperature(optarg);
			if (hfc.zippingTemperature < 0)
				die("Invalid zipping temperature '%s'\n", optarg);
			break;
		case 'Q':
			hfc.unzippingTemperature = parseTemperature(optarg);
			if (hfc.unzippingTemperature < 0)
				die("Invalid unzipping temperature '%s'\n", optarg);
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
		case 'e':
			measureEndToEndDistance = true;
			break;
		case 'p':
			measureBasePairing = true;
			break;
		case 'k':
			measureTemperature = true;
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
		worldSize = ((strlen(baseSequence) + 2)
					* DEF_MONOMER_WORLDSIZE_FACTOR) * ANGSTROM;

	if (verletSettings.tau < 0)
		verletSettings.tau = DEF_COUPLING_TIMESTEP_FACTOR
						* integratorConf.timeStep;

	if (interactionSettings.truncationLen < 0) {
		/* Disable truncation -> no space partitioning */
		printf("Disabling space partitioning, "
				"maximizing truncation length.\n");
		integratorConf.numBoxes = 1;
		interactionSettings.truncationLen = worldSize / 2.0;
		/* Due to (cubic) periodicity, we still need to truncate at 
		 * worldsize/2 to have correct energy conservation and to 
		 * make sure every particle 'sees' the same (spherical) 
		 * potential, no matter where it is within the cube. */
	} else if (interactionSettings.truncationLen > worldSize / 2.0) {
		/* World is too small, extend it so we correctly compute 
		 * the potentials and forces up to the requested truncation 
		 * length. We need twice the truncation length for the same 
		 * reason as above. */
		printf("Truncation (%e) > worldSize/2 (%e)\n   => "
				"Extending worldSize to 2*Truncation (%e).\n",
				interactionSettings.truncationLen, worldSize/2.0,
				2.0 * interactionSettings.truncationLen);
		worldSize = 2.0 * interactionSettings.truncationLen;
		integratorConf.numBoxes = 1; /* No space partitioning */
	} else if (integratorConf.numBoxes == -1) {
		/* Automatically determine ideal number of boxes.
		 * Formula is fitted for N between 100 and 1000 on a Core2 
		 * Duo E8500.
		 * (+ 0.5 for correct rounding to int) */
		int ideal = 0.5 + 2.73 * pow(strlen(baseSequence), 0.446);
		integratorConf.numBoxes = MIN(ideal,
					worldSize / interactionSettings.truncationLen);
		if (integratorConf.numBoxes < 1)
			integratorConf.numBoxes = 1;
		printf("Number of boxes per dimension: %d\n",
				integratorConf.numBoxes);
	}

	if (worldSize / integratorConf.numBoxes < interactionSettings.truncationLen)
		die("The boxsize (%e) is smaller than the potential "
			"truncation radius (%e)!\n",
			worldSize / integratorConf.numBoxes / ANGSTROM,
			interactionSettings.truncationLen / ANGSTROM);

	renderConf.numBoxes = integratorConf.numBoxes;
}

int main(int argc, char **argv)
{
	seedRandom();

	parseArguments(argc, argv);
	const char *filenameBase = measurementConf.measureFile;

	if (buildCompStrand)
		allocWorld(2, worldSize);
	else
		allocWorld(1, worldSize);

	fillStrand(&world.strands[0], baseSequence);
	if (buildCompStrand)
		fillComplementaryStrand(&world.strands[1], baseSequence);

	killMomentum();

	assert(worldSanityCheck());

	/* Integrator config */
	Integrator integrator;
	integrator.type = integratorType;
	switch (integratorType) {
	case VERLET:
		integrator.settings.verlet = verletSettings;
		break;
	case LANGEVIN:
		integrator.settings.langevin = langevinSettings;
		break;
	default:
		assert(false); die("Unknown integrator type!\n");
	}
	integratorConf.integrator = integrator;

	/* Measurement header */
	char *measHeaderStrings[2];
	measHeaderStrings[0] = getWorldInfo();
	measHeaderStrings[1] = integratorInfo(&integratorConf);
	char *measHeader = asprintfOrDie("%s%s",
				measHeaderStrings[0], measHeaderStrings[1]);
	free(measHeaderStrings[0]); free(measHeaderStrings[1]);
	measurementConf.measureHeader = measHeader;

	/* Measurement config for additional measurements (end to end, base 
	 * pairing, temperature, ...) */
	MeasurementConf additionalMeasConf = measurementConf; /* struct copy */
	additionalMeasConf.verbose = false; /* Let output come from main 
					       measurement. */
	additionalMeasConf.measureFile = NULL; /* Don't overwrite! */
	additionalMeasConf.measureWait = 0; /* Start sampling immediately 
					       for additional measurements. */

	/* Integrator task */
	Task integratorTask = makeIntegratorTask(&integratorConf);

	/* Verbose task */
	Measurement verbose;
	verbose.measConf = verboseConf;
	verbose.sampler = dumpStatsSampler();
	Task verboseTask = measurementTask(&verbose);

	/* Base pairing task */
	char *basePairFile = asprintfOrDie("%s%s", filenameBase,
						BASE_PAIRING_FILE_SUFFIX);
	Measurement basePairing;
	basePairing.sampler = basePairingSampler(&bpc);
	basePairing.measConf = additionalMeasConf; /* struct copy */
	basePairing.measConf.measureFile = basePairFile;
	Task basePairingTask = measurementTask(&basePairing);

	/* End to end task */
	char *endToEndFile = asprintfOrDie("%s%s", filenameBase,
						END_TO_END_DIST_FILE_SUFFIX);
	Measurement endToEnd;
	endToEnd.sampler = endToEndDistSampler(&world.strands[0]);
	endToEnd.measConf = additionalMeasConf; /* struct copy */
	endToEnd.measConf.measureFile = endToEndFile;
	Task endToEndTask = measurementTask(&endToEnd);	

	/* Temperature task */
	char *temperatureFile = asprintfOrDie("%s%s", filenameBase,
						TEMPERATURE_FILE_SUFFIX);
	Measurement tempMeas;
	tempMeas.sampler = temperatureSampler();
	tempMeas.measConf = additionalMeasConf; /* struct copy */
	tempMeas.measConf.measureFile = temperatureFile;
	Task temperatureTask = measurementTask(&tempMeas);	


	/* Hairpin task */
	Measurement hairpin;
	switch (hairpinMeasurementType) {
	case HAIRPIN_MELTING_TEMPERATURE:
		hairpin.sampler = hairpinMeltingTempSampler(&hmtc);
		break;
	case HAIRPIN_FORMATION_TIME:
		if (hfc.requiredBoundBPs < 0) /* guard */
			die("You need to give the minimum required number "
					"of bound base pairs to measure "
					"hairpin formation time!\n");
		hairpin.sampler = hairpinFormationSampler(&hfc);
		break;
	case HAIRPIN_NO_MEASUREMENT:
		hairpin.sampler = trivialSampler();
		break;
	default:
		assert(false); die("Unknown hairpin measurement type!\n");
	}
	hairpin.measConf = measurementConf;
	Task hairpinTask = measurementTask(&hairpin);

	/* Render task */
	Task renderTask = makeRenderTask(&renderConf);

	/* Combined task */
	Task *tasks[7];
	tasks[0] = (render ? &renderTask : NULL);
	tasks[1] = &integratorTask;
	tasks[2] = &verboseTask;
	tasks[3] = &hairpinTask;;
	tasks[4] = (measureEndToEndDistance ? &endToEndTask : NULL);
	tasks[5] = (measureBasePairing ? &basePairingTask : NULL);
	tasks[6] = (measureTemperature ? &temperatureTask : NULL);
	Task task = sequence(tasks, 7);

	setHeatBathTemperature(temperature);
	registerInteractionSettings(interactionSettings);
	initPhysics();

	bool everythingOK = run(&task);

	freeWorld();
	free(basePairFile);
	free(endToEndFile);
	free(temperatureFile);
	free(measHeader);

	if (!everythingOK)
		return 1;
	return 0;
}


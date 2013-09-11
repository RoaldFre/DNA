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
#include "spgrid.h"
#include "integrator.h"
#include "extraInteractions.h"
#include "monteCarlo.h"
#include "render.h"
#include "samplers.h"
#include "math.h"

#define DEF_DATA_PATH "hairpinData"
#define END_TO_END_DIST_FILE_SUFFIX "_endToEnd"
#define BASE_PAIRING_FILE_SUFFIX "_basePairing"
#define DISTANCE_BASE_PAIRING_FILE_SUFFIX "_basePairingDistance"
#define ENERGY_BASE_PAIRING_FILE_SUFFIX "_basePairingEnergy"
#define TEMPERATURE_FILE_SUFFIX "_temperature"
#define GYRATION_FILE_SUFFIX "_gyrationRadius"
#define FORCE_VEL_FRIC_FILE_SUFFIX "_forceVelFric"
#define FORCE_VEL_FRIC_P_FILE_SUFFIX "_forceVelFricP"
#define FORCE_VEL_P_FILE_SUFFIX "_forceVelP"
#define FORCE_VEL_S_FILE_SUFFIX "_forceVelS"
#define STATE_SUFFIX "_state"
#define ZIPPED_STATE_FILE_SUFFIX "_zippedState"

//#define ENERGY_THRESHOLD (-0.1 * EPSILON)
#define ENERGY_THRESHOLD (-EPSILON)

/* Defaults */
#define DEF_BASE_SEQUENCE		"GCCTATTTTTTAATAGGC" /* N=4 in Kuznetsov nov 2001 */
#define DEF_TIMESTEP 			15.0
#define DEF_REBOX_INTERVAL		100.0
#define DEF_INITIAL_TEMPERATURE		"20C"
#define DEF_SALT_CONCENTRATION		100.0  /* mol/m^3 */
#define DEF_LANGEVIN_GAMMA		5e12
#define DEF_COUPLING_TIMESTEP_FACTOR 	1000
#define DEF_TRUNCATION_LENGTH		20.0
#define DEF_MONOMER_WORLDSIZE_FACTOR    5.0
#define DEF_RENDER_FRAMERATE 		30.0
#define DEF_RENDER_RADIUS 		0.5
#define DEF_INTEGRATOR			LANGEVIN


/* Static global configuration variables */

static const char* initialStateFile = NULL; /* Read world state from here. */
static const char* finalStateFile = NULL; /* Dump world state here at end. */
static const char* randomSeedSource = NULL;
static const char* randomSeedDest = NULL;
static MeasurementConf verboseConf =
{
	.maxMeasureTime = -1, /* loop forever */
	.minMeasureTime = 0,
	.measureInterval = -1, /* default: disable */
	.measureWait = 0,
	.measureFile = NULL,
	.verbose = false,
};
/* generic measurement config for whatever measurement we will do */
static MeasurementConf measurementConf =
{
	.maxMeasureTime = -1,
	.minMeasureTime = 0,
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
	.energyThreshold = ENERGY_THRESHOLD,
	.mode = BPSM_BOOLEAN,
};
static HairpinFormationSamplerConfig hfc =
{
	.energyThreshold = ENERGY_THRESHOLD,
	.requiredBoundBPs = -1, /* guard */
	.nucleationBoundBPs = -1, /* guard */
	.zipConfirmationTime = 0,
	.allowedBoundBPs = 0,
	.unzipConfirmationTime = 1 * NANOSECONDS,
	.zippingTemperature = CELSIUS(20),
	.unzippingTemperature = CELSIUS(90),
	.zippedRelaxationTime = 5 * NANOSECONDS,
	.minZippingSamplingTime = 0,
	.zippedStateFile = NULL,
};
static HairpinMeltingTempSamplerConfig hmtc =
{
	.energyThreshold = ENERGY_THRESHOLD,
	.Tstart = CELSIUS(20),
	.Tstep = 10,
	.numSteps = 10,
	.relaxationTime = 5 * NANOSECONDS,
	.measureTime = 100 * NANOSECONDS,
	.verbose = false,
};
static HairpinStateSamplerConfig hsc =
{
	.energyThreshold = ENERGY_THRESHOLD,
	.temperature = -1.0,
	.dumpDistances = false,
};
static HairpinConfigSamplerConfig hcsc =
{
	.energyThreshold = ENERGY_THRESHOLD,
};
static enum {
	HAIRPIN_MELTING_TEMPERATURE,
	HAIRPIN_FORMATION_TIME,
	HAIRPIN_STATE,
	HAIRPIN_NO_MEASUREMENT,
} hairpinMeasurementType = HAIRPIN_NO_MEASUREMENT;

static bool measureEndToEndDistance = false;
static bool measureBasePairing = false;
static bool measureDistanceBasePairing = false;
static bool measureEnergyBasePairing = false;
static bool measureTemperature = false;
static bool measureGyrationRadius = false;
static bool measureForceVelFric = false;

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
	.reboxInterval = DEF_REBOX_INTERVAL * FEMTOSECONDS,
};
static MonteCarloMove moves[] = {
		{.weight = 1.0, .m = &jiggleMover},
		{.weight = 5.0, .m = &pivotMover},
};
static MonteCarloConfig monteCarloConfig = {
	.sweeps = -1,
	.verbose = true,
	.moves = {
		  .numMoves = 2,
	          .moves = moves,
	         }
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
	.onlyXYbasePairing = false,
	.onlyRepulsiveBasePairing = false,
	.basePairFactor = 1.0,
	.basePairInteraction = BASE_PAIR_HAIRPIN,
	.saltConcentration   = DEF_SALT_CONCENTRATION,
	.truncationLen       = DEF_TRUNCATION_LENGTH * ANGSTROM,
};
static HarmonicEndToEndInt hetei = {
	.K = -1, /* guard */
	.Rref = -1,
	.s = NULL,
};
static bool harmonicEndToEndInteraction = false;
static EndToEndForceInt etefi = {
	.F = {0, 0, 0}, /* guard */
	.s = NULL,
};
static bool endToEndForceInteraction = false;

static const char* baseSequence = DEF_BASE_SEQUENCE;
static double worldSize = -1; /* guard */
static int numBoxes = -1; /* guard */
static bool useMonteCarlo = false;
static bool noDynamics = false;


static void printUsage(void)
{
	printf("Usage: main [flags]\n");
	printf("\n");
	printf("Flags:\n");
	printf(" -s <str>  base Sequence of the DNA strand to simulate\n");
	printf("             default: %s\n", DEF_BASE_SEQUENCE);
	printf(" -t <flt>  length of Time steps (in femtoseconds)\n");
	printf("             default: %lf\n", DEF_TIMESTEP);
	printf(" -T <flt><C|K>  initial Temperature (example: 20C or 300K)\n");
	printf("             note: some measurements set their own temperature when sampling!\n");
	printf("             default: %s\n", DEF_INITIAL_TEMPERATURE);
	printf(" -N <flt>  concentration of Na+ in the environment (in mol/m^3)\n");
	printf("             default: %lf\n", DEF_SALT_CONCENTRATION);
	printf(" -m <num>  perform the given number of Monte carlo sweeps and then quit\n");
	printf(" -j <type> base pairing interaction only between bases of <type>:\n");
	printf("             h: matching bases according to hairpin position [default]\n");
	printf("             a: all matching bases\n");
	printf("             x: behave like 'h' for XY pairs, and 'a' for the other bases\n");
	printf(" -Y        only enable base pairing between XY base pairs\n");
	printf(" -o <fact> disable dihedral interaction and modify base pairing interaction by the given factor\n");
	printf(" -n        No attractive base pairing interaction: only repulsive part\n");
	printf(" -y <path> write random number generator seed to this file\n");
	printf(" -z <path> read random number generator seed from this file\n");
	printf(" -r        Render\n");
	printf(" -f <flt>  desired Framerate when rendering.\n");
	printf("             default: %lf)\n", DEF_RENDER_FRAMERATE);
	printf(" -R <flt>  Radius (in Angstrom) of the particles when rendering\n");
	printf("             default: %lf\n", DEF_RENDER_RADIUS);
	printf(" -F        draw Forces on particles when rendering\n");
	printf(" -l <flt>  truncation Length of potentials (in Angstrom).\n");
	printf("             negative value: sets truncation to worldsize/2\n");
	printf("             default: %lf\n", DEF_TRUNCATION_LENGTH);
	printf(" -S <flt>  Size of world (in Angstrom).\n");
	printf("             default: (number of monomers + 2) * %lf\n", DEF_MONOMER_WORLDSIZE_FACTOR);
	printf(" -b <num>  number of Boxes per dimension\n");
	printf("             default: max so that boxsize >= potential truncation length\n");
	printf(" -x <flt>  time interval between reboXing (in femtoseconds)\n");
	printf("             default: %lf\n", DEF_REBOX_INTERVAL);
	printf(" -v <int>  Verbose: dump statistics every <flt> picoseconds\n");
	printf(" -i <type> Integrator to use. Values for <type>:\n");
	printf("             l: Langevin (velocity BBK) [default]\n");
	printf("             v: velocity Verlet with Berendsen thermostat\n");
	printf("             m: use Monte Carlo sampling\n");
	printf("             n: No dynamics (useful for just viewing the world when rendering)\n");
	printf("\n");
	printf("Parameters for Langevin integrator:\n");
	printf(" -g <flt>  Gamma: friction coefficient for Langevin dynamics\n");
	printf("             default: %le\n", DEF_LANGEVIN_GAMMA);
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
	printf(" -P <flt>  measurement Period: maximum time to sample the system (in nanoseconds)\n");
	printf("             default: sample indefinitely\n");
	printf(" -K <flt>  Keep sampling for the given minimum time (in nanoseconds).\n");
	printf("             default: 0\n");
	printf(" -D <path> Data file to Dump measurement output. The directory must exist.\n");
	printf("             default: %s\n", DEF_DATA_PATH);
	printf(" -w <path> data file to dump World state in at the end of the measurement. The directory must exist.\n");
	printf(" -d <path> read initial world Data from file\n");
	printf(" -X <type>  measurement to eXecute. Values for <type>:\n");
	printf("             m: hairpin Melting temperature\n");
	printf("             f: hairpin Formation time\n");
	printf("             s: hairpin State\n");
	printf(" -e        also measure End-to-end distance of the strand\n");
	printf("             output: the data filename (see -D) with suffix: '%s'\n",
							END_TO_END_DIST_FILE_SUFFIX);

	printf(" -u <K>:<Rref>  Umbrella potential on end-to-end distance R: K/2(R - Rref)^2.\n");
	printf("           For Umbrella sampling of end-to-end distance free energy profile.\n");
	printf(" -E <Fx>:<Fy>:<Fz>  constant End-to-end force F.\n");
	printf(" -p        also measure raw hairpin base Pairing state\n");
	printf("             output: the data filename (see -D) with suffix: '%s'\n",
							BASE_PAIRING_FILE_SUFFIX);
	printf("                     the distance file gets suffix: '%s'\n",
							DISTANCE_BASE_PAIRING_FILE_SUFFIX);
	printf("                     the energy file gets suffix: '%s'\n",
							ENERGY_BASE_PAIRING_FILE_SUFFIX);
	printf(" -k        also measure the Kinetic temperature\n");
	printf("             output: the data filename (see -D) with suffix: '%s'\n",
							TEMPERATURE_FILE_SUFFIX);
	printf(" -q        also measure the gyration radius of the strand\n");
	printf("             output: the data filename (see -D) with suffix: '%s'\n",
							GYRATION_FILE_SUFFIX);
	printf(" -J        also measure the force, velocity, and friction for each monomer of the strand\n");
	printf("             output: the data filename (see -D) with suffix: '%s' for the monomers\n",
							FORCE_VEL_FRIC_FILE_SUFFIX);
	printf("                     the data filename (see -D) with suffix: '%s' for the backbone phosphates\n",
							FORCE_VEL_FRIC_P_FILE_SUFFIX);
	printf("                     the data filename (see -D) with suffix: '%s' for the phosphates (relative to backbone)\n",
							FORCE_VEL_P_FILE_SUFFIX);
	printf("                     the data filename (see -D) with suffix: '%s' for the sugars (relative to backbone)\n",
							FORCE_VEL_S_FILE_SUFFIX);
	printf("\n");
	printf("For more info and default values of params below: see code :P\n");
	//TODO parameter names are a mess
	printf("Parameters for hairpin melting temperature measurement:\n");
	printf(" -A <flt><C|K> startTemp\n");
	printf("             default: %lfC\n", TO_CELSIUS(hmtc.Tstart));
	printf(" -B <flt>  stepTemp\n");
	printf("             default: %lf\n", hmtc.Tstep);
	printf(" -C <int>  nSteps\n");
	printf("             default: %d\n", hmtc.numSteps);
	printf(" -G <flt>  measureTime per step (nanoseconds)\n");
	printf("             default: %lf\n", hmtc.measureTime / NANOSECONDS);
	printf(" -L <flt>  reLaxTime per step (nanoseconds)\n");
	printf("             default: %lf\n", hmtc.relaxationTime / NANOSECONDS);
	printf(" -V        Be verbose: also dump hairpin state to file\n");
	printf("Parameters for hairpin formation measurement:\n");
	printf(" -H <int>:<int>  min required bounded base pairs for confirming 'zipped' state,\n");
	printf("           followed by min required bounded base pairs for confirming nucleation\n");
	printf(" -M <int>  max allowed bounded base pairs for confirming 'unzipped' state\n");
	printf("             default: %d\n", hfc.allowedBoundBPs);
	printf(" -O <flt><C|K> zipping temperature\n");
	printf("             default: %lfC\n", TO_CELSIUS(hfc.zippingTemperature));
	printf(" -Q <flt><C|K> unzipping temperature\n");
	printf("             default: %lfC\n", TO_CELSIUS(hfc.unzippingTemperature));
	printf(" -U <flt>  zipped relaxation phase duration (nanoseconds)\n");
	printf("             default: %lf\n", hfc.zippedRelaxationTime / NANOSECONDS);
	printf(" -Z <flt>  minimum duration of Zipping+relaxation phase (nanoseconds).\n");
	printf("           The relaxation phase duration is extended when required.\n");
	printf("             default: %lf\n", hfc.minZippingSamplingTime / NANOSECONDS);
	printf("Parameters for hairpin state measurement:\n");
	printf(" -a <flt><C|K> startTemp\n");
	printf("             default: same as initial temperature\n");
	printf("\n");
}

static void parseArguments(int argc, char **argv)
{
	int c;
	const char *measDescr;

	/* defaults */
	temperature = parseTemperature(DEF_INITIAL_TEMPERATURE);

	/* Unused options:
	 * [All letters have been used >_<] */
	/* TODO rework option system to something sane. Long options? Read 
	 * and parse config file? */
	while ((c = getopt(argc, argv, ":s:t:T:N:m:Yo:j:ny:z:g:c:f:rR:Fl:S:b:x:v:i:W:I:P:K:D:w:d:X:eu:E:pkqJhA:B:C:G:L:VH:M:O:Q:U:Z:a:")) != -1)
	{
		switch (c)
		{
		case 's':
			baseSequence = optarg;
			printf("s: Setting sequence to %s\n", baseSequence);
			break;
		case 't':
			integratorConf.timeStep = atof(optarg) * FEMTOSECONDS;
			if (integratorConf.timeStep <= 0)
				die("Invalid timestep %s\n", optarg);
			printf("t: Setting time step %le\n", integratorConf.timeStep);
			break;
		case 'T':
			temperature = parseTemperature(optarg);
			if (temperature < 0)
				die("Invalid temperature!\n");
			printf("T: Setting temperature to %lf\n", temperature);
			break;
		case 'N':
			interactionSettings.saltConcentration = atof(optarg);
			if (interactionSettings.saltConcentration < 0)
				die("Invalid salt concentration %s\n", optarg);
			printf("N: Setting salt concentration to %lf\n", 
					interactionSettings.saltConcentration);
			break;
		case 'm':
			monteCarloConfig.sweeps = atoi(optarg);
			if (monteCarloConfig.sweeps < 0)
				die("Invalid number of monte carlo sweeps %s\n", optarg);
			printf("m: Stopping after %d monte carlo sweeps\n", 
					monteCarloConfig.sweeps);
			break;
		case 'Y':
			interactionSettings.onlyXYbasePairing = true;
			printf("Y: Only enabling base pairing between XY pairs\n");
			break;
		case 'o':
			interactionSettings.enableDihedral = false;
			interactionSettings.basePairFactor = atof(optarg);
			printf("o: Disabling dihedral interation and modifying base pairing with factor %f\n",
					interactionSettings.basePairFactor);
			break;
		case 'j':
			if (optarg[0] == '\0' || optarg[1] != '\0')
				die("Base pairing: badly formatted base pairing type\n");
			const char *basePairingDescr;
			enum BasePairInteraction bpi;
			switch(optarg[0]) {
			case 'h':
				bpi = BASE_PAIR_HAIRPIN;
				basePairingDescr = "hairpin";
				break;
			case 'a':
				bpi = BASE_PAIR_ALL;
				basePairingDescr = "all";
				break;
			case 'x':
				bpi = BASE_PAIR_XY_HAIRPIN_REST_ALL;
				basePairingDescr = "XY:hairpin, rest:all";
				break;
			default: die("Unknown base pairing type '%s'\n", optarg);
				basePairingDescr = "ERROR"; bpi = 0;
				break;
			}
			printf("j: Setting base pairing interaction to: %s\n", 
						basePairingDescr);
			interactionSettings.basePairInteraction = bpi;
			break;
		case 'n':
			interactionSettings.onlyRepulsiveBasePairing = true;
			printf("n: Disabling base pairing\n");
			break;
		case 'y':
			randomSeedDest = optarg;
			printf("y: writing random number generator seed to %s\n",
							randomSeedDest);
			break;
		case 'z':
			randomSeedSource = optarg;
			printf("z: reading random number generator seed from %s\n",
							randomSeedSource);
			break;
		case 'g':
			langevinSettings.gamma = atof(optarg);
			if (langevinSettings.gamma < 0)
				die("Invalid friction coefficient %s\n", optarg);
			printf("g: Setting gamma to %le\n", 
						langevinSettings.gamma);
			break;
		case 'c':
			verletSettings.tau = atof(optarg) * FEMTOSECONDS;
			if (verletSettings.tau < 0)
				die("Invalid thermostat relaxation time %s\n",
						optarg);
			printf("c: Setting verlet coupling to %le\n", 
						verletSettings.tau);
			break;
		case 'f':
			renderConf.framerate = atof(optarg);
			if (renderConf.framerate < 0)
				die("Invalid framerate %s\n", optarg);
			printf("f: Setting frame rate to %d\n", 
						renderConf.framerate);
			break;
		case 'r':
			render = true;
			printf("r: Enabling rendering\n");
			break;
		case 'R':
			renderConf.radius = atof(optarg) * ANGSTROM;
			if (renderConf.radius <= 0)
				die("Invalid radius %s\n", optarg);
			printf("R: Setting render radius to %le\n", 
						renderConf.radius);
			break;
		case 'F':
			renderConf.drawForces = true;
			printf("F: Enabling rendering of forces\n");
			break;
		case 'l':
			interactionSettings.truncationLen = atof(optarg) * ANGSTROM;
			printf("l: Setting truncationLength to %le\n", 
					interactionSettings.truncationLen);
			break;
		case 'S':
			worldSize = atof(optarg) * ANGSTROM;
			if (worldSize <= 0)
				die("Invalid world size %s\n", optarg);
			printf("S: Setting world size to %le\n", worldSize);
			break;
		case 'b':
			numBoxes = atoi(optarg);
			if (numBoxes <= 0)
				die("Invalid number of boxes %s\n",
						optarg);
			printf("b: Setting number of boxes to %d\n", numBoxes);
			break;
		case 'x':
			integratorConf.reboxInterval = atof(optarg) * FEMTOSECONDS;
			if (integratorConf.reboxInterval <= 0)
				die("Invalid rebox interval %s\n", optarg);
			printf("x: Setting rebox interval to %le\n", 
						integratorConf.reboxInterval);
			break;
		case 'v':
			verboseConf.measureInterval = atof(optarg) * PICOSECONDS;
			if (verboseConf.measureInterval <= 0)
				die("Verbose: invalid verbose interval %s\n",
						optarg);
			printf("v: Setting verbose interval to %le\n", 
						verboseConf.measureInterval);
			break;
		case 'i':
			if (optarg[0] == '\0' || optarg[1] != '\0')
				die("Integrator: badly formatted integrator type\n");
			const char *integratorDescr;
			switch(optarg[0]) {
			case 'l': 
				integratorType = LANGEVIN;
				integratorDescr = "Langevin";
				break;
			case 'v':
				integratorType = VERLET;
				integratorDescr = "Verlet";
				break;
			case 'm':
				useMonteCarlo = true;
				integratorDescr = "Monte Carlo";
				break;
			case 'n':
				noDynamics = true;
				integratorDescr = "No Dynamics";
				break;
			default: die("Unknown integrator type '%s'\n", optarg);
				integratorDescr = "ERROR";
				break;
			}
			printf("i: Integrating with %s integrator\n", 
						integratorDescr);
			break;
		case 'W':
			measurementConf.measureWait = atof(optarg) * NANOSECONDS;
			if (measurementConf.measureWait < 0)
				die("Invalid relaxation time %s\n", optarg);
			printf("W: Setting measurement wait to %le\n", 
						measurementConf.measureWait);
			break;
		case 'I':
			measurementConf.measureInterval = atof(optarg) * PICOSECONDS;
			if (measurementConf.measureInterval<= 0)
				die("Invalid measurement interval %s\n", optarg);
			printf("I: Setting measurement interval to %le\n", 
						measurementConf.measureInterval);
			break;
		case 'P':
			measurementConf.maxMeasureTime = atof(optarg) * NANOSECONDS;
			if (measurementConf.maxMeasureTime <= 0)
				die("Invalid measurement time %s\n", optarg);
			printf("P: Setting maximum measurement time to %le\n", 
						measurementConf.maxMeasureTime);
			break;
		case 'K':
			measurementConf.minMeasureTime = atof(optarg) * NANOSECONDS;
			if (measurementConf.minMeasureTime <= 0)
				die("Invalid measurement time %s\n", optarg);
			printf("P: Setting minimum measurement time to %le\n", 
						measurementConf.minMeasureTime);
			break;
		case 'D':
			measurementConf.measureFile = optarg;
			printf("D: Setting measurement file to %s\n", 
						measurementConf.measureFile);
			break;
		case 'w':
			finalStateFile = optarg;
			printf("w: Setting final state file to %s\n", 
						finalStateFile);
			break;
		case 'd':
			initialStateFile = optarg;
			printf("d: Setting initial state file to %s\n", 
						initialStateFile);
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
			printf("A: Setting melting starting temperature to %lf\n", 
						hmtc.Tstart);
			break;
		case 'B':
			hmtc.Tstep = atof(optarg);
			printf("B: Setting melting temperature step to %lf\n", 
						hmtc.Tstep);
			break;
		case 'C':
			hmtc.numSteps = atoi(optarg);
			printf("C: Setting melting temperature number of steps to %d\n", 
						hmtc.numSteps);
			break;
		case 'G':
			hmtc.measureTime = atof(optarg) * NANOSECONDS;
			printf("G: Setting melting measure time to %le\n", 
						hmtc.measureTime);
			break;
		case 'L':
			hmtc.relaxationTime = atof(optarg) * NANOSECONDS;
			printf("L: Setting melting relaxation time to %le\n", 
						hmtc.relaxationTime);
			break;
		case 'V':
			hmtc.verbose = true;
			printf("v: Being verbose for melting\n");
			break;
		case 'H':
			//TODO error handling
			sscanf(optarg, "%d:%d", &hfc.requiredBoundBPs, &hfc.nucleationBoundBPs);
			printf("H: Setting formation required bound to %d and nucleation bound to %d\n", 
						hfc.requiredBoundBPs, hfc.nucleationBoundBPs);
			break;
		case 'M':
			hfc.allowedBoundBPs = atoi(optarg);
			printf("M: Setting formation allowed bound to %d\n", 
						hfc.allowedBoundBPs);
			break;
		case 'O':
			hfc.zippingTemperature = parseTemperature(optarg);
			if (hfc.zippingTemperature < 0)
				die("Invalid zipping temperature '%s'\n", optarg);
			printf("O: Setting formation zipping temperature to %lf\n", 
						hfc.zippingTemperature);
			break;
		case 'Q':
			hfc.unzippingTemperature = parseTemperature(optarg);
			if (hfc.unzippingTemperature < 0)
				die("Invalid unzipping temperature '%s'\n", optarg);
			printf("Q: Setting formation unzipping temperature to %lf\n", 
						hfc.unzippingTemperature);
			break;
		case 'U':
			hfc.zippedRelaxationTime = atof(optarg) * NANOSECONDS;
			printf("U: Setting formation zipped relaxation time to %le\n", 
						hfc.zippedRelaxationTime);
			break;
		case 'Z':
			hfc.minZippingSamplingTime = atof(optarg) * NANOSECONDS;
			printf("Z: Setting min zipping+relax sampling time to %le\n", 
						hfc.minZippingSamplingTime);
			break;
		case 'a':
			hsc.temperature = parseTemperature(optarg);
			if (hsc.temperature < 0)
				die("Invalid starting temperature '%s'\n", optarg);
			printf("a: Setting state measurement starting temperature to %le\n", 
						hsc.temperature);
			break;
		case 'X':
			switch(optarg[0]) {
			case 'f': 
				hairpinMeasurementType = HAIRPIN_FORMATION_TIME;
				measDescr = "hairpin formation time";
				break;
			case 'm':
				hairpinMeasurementType = HAIRPIN_MELTING_TEMPERATURE;
				measDescr = "hairpin melting temperature";
				break;
			case 's':
				hairpinMeasurementType = HAIRPIN_STATE;
				measDescr = "hairpin state";
				break;
			default:
				die("Unknown hairpin measurement type '%s'\n", optarg);
				measDescr = "ERROR";
				break;
			}
			printf("X: hairpin measurement: %s\n", measDescr);
			break;
		case 'e':
			measureEndToEndDistance = true;
			printf("e: measuring end to end ditance\n");
			break;
		case 'u':
			//TODO error handling
			sscanf(optarg, "%lf:%lf", &hetei.K, &hetei.Rref);
			harmonicEndToEndInteraction = true;
			hetei.K *= ELECTRON_VOLT / SQUARE(ANGSTROM);
			hetei.Rref *= ANGSTROM;
			measureEndToEndDistance = true;
			printf("u: Umbrella potential on end-to-end distance with K=%le and Rref=%le\n",
					hetei.K, hetei.Rref);
			break;
		case 'E':
			//TODO error handling
			sscanf(optarg, "%lf:%lf:%lf", &etefi.F.x, &etefi.F.y, &etefi.F.z);
			endToEndForceInteraction = true;
			printf("E: Constant end-to-end force: ");
			printVectorExp(etefi.F);
			printf("\n");
			break;
		case 'p':
			//TODO: Clean this up asap
			//measureBasePairing = true;
			measureDistanceBasePairing = true;
			measureEnergyBasePairing = true;
			printf("p: measuring base pairing\n");
			break;
		case 'k':
			measureTemperature = true;
			printf("k: measuring temperature\n");
			break;
		case 'q':
			measureGyrationRadius = true;
			printf("q: measuring gyration radius\n");
			break;
		case 'J':
			measureForceVelFric = true;
			printf("J: measuring force, velocity and friction per monomer\n");
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

	if (verletSettings.tau < 0)
		verletSettings.tau = DEF_COUPLING_TIMESTEP_FACTOR
						* integratorConf.timeStep;
}

/* Ideal number of boxes: as many as we can cram into the world, without 
 * getting the box size smaller than the truncation distance. This may not 
 * be beneficial with large world sizes (necessary for long strands), as 
 * the memory complexity scales roughly as N^3 (assuming world size scales 
 * as N). Nonetheless, we only visit occupied boxes, so time complexity is 
 * linear for looping over pairs. */
static void determineIdealNumberOfBoxes(void)
{
	int Nmonomers = strlen(baseSequence);
	double defaultWorldSize = ((Nmonomers  + 2)
				* DEF_MONOMER_WORLDSIZE_FACTOR) * ANGSTROM;
	if (worldSize < 0)
		worldSize = defaultWorldSize;

	if ((harmonicEndToEndInteraction || endToEndForceInteraction || measureEndToEndDistance)
			&&   worldSize < 3 * MAX(defaultWorldSize, hetei.Rref)) {
		worldSize = 3 * MAX(defaultWorldSize, hetei.Rref);
		printf("Need correct end-to-end distance. Increasing world "
				"size to %lf Angstrom\n",
				worldSize / ANGSTROM);
	}

	if (interactionSettings.truncationLen < 0) {
		/* Disable truncation -> no space partitioning */
		printf("Disabling space partitioning, "
				"maximizing truncation length.\n");
		numBoxes = 1;
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
		printf("Truncation (%le) > worldSize/2 (%le)\n   => "
				"Extending worldSize to 2*Truncation (%le).\n",
				interactionSettings.truncationLen, worldSize/2.0,
				2.0 * interactionSettings.truncationLen);
		worldSize = 2.0 * interactionSettings.truncationLen;
		numBoxes = 1; /* No space partitioning */
	} else if (numBoxes == -1) {
		/* Automatically determine ideal number of boxes. */
		numBoxes = worldSize / interactionSettings.truncationLen;
		if (numBoxes < 1)
			numBoxes = 1;
		printf("Automatically determined number of boxes per dimension: %d\n",
				numBoxes);
	}

	if (smallerThanEpsilon(worldSize / numBoxes,
	                       interactionSettings.truncationLen,
			       1e-5))
		die("The boxsize (%le) is smaller than the potential "
			"truncation radius (%le)!\n",
			worldSize / numBoxes / ANGSTROM,
			interactionSettings.truncationLen / ANGSTROM);

	renderConf.numBoxes = numBoxes;
	printf("Worldsize: %e\nNumber of boxes: %d\nBox volume: %e\n",
			worldSize, numBoxes, CUBE(worldSize / numBoxes));
}

int main(int argc, char **argv)
{
	printf("Hairpin, version: %s\n", VERSION);

	parseArguments(argc, argv);

	uint64_t seed = 0; /* initialize for compiler warning */
	if (randomSeedSource == NULL) {
		seed = seedRandom();
	} else {
		if (!readSeed(randomSeedSource, &seed))
			die("Couldn't read random generator seed file: %s\n", 
							randomSeedSource);
		seedRandomWith(seed);
	}

	if (randomSeedDest != NULL)
		if (!writeSeed(randomSeedDest, seed))
			die("Couldn't write random generator seed file: %s\n", 
							randomSeedDest);

	printf("Random number generator seed: %llx\n", (long long unsigned) seed);

	const char *filenameBase = measurementConf.measureFile;

	/* Set temperature before constructing world so the generated 
	 * particles get seeded with appropriate velocities. */
	setHeatBathTemperature(temperature);

	if (initialStateFile != NULL) {
		readWorld(initialStateFile);

		if (worldSize < 0)
			worldSize = world.worldSize; /* from readWorld() */
		else
			resizeWorld(worldSize); /* User input */

		determineIdealNumberOfBoxes();
		/* Resize again in case determineIdealNumberOfBoxes() said 
		 * we need to go bigger */
		resizeWorld(worldSize);
	} else {
		determineIdealNumberOfBoxes();
		allocWorld(1, worldSize);
		fillStrand(&world.strands[0], baseSequence);
		killMomentum();
	}
	initGrid(numBoxes);

	centerStrand(&world.strands[0]);

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
		die("Unknown integrator type!\n");
	}
	integratorConf.integrator = integrator;

	/* End to end interactions */
	if (harmonicEndToEndInteraction) {
		hetei.s = &world.strands[0];
		registerHarmonicEndToEndInt(&hetei);
	}
	if (endToEndForceInteraction) {
		etefi.s = &world.strands[0];
		registerEndToEndForceInt(&etefi);
	}


	/* Reconstruct command line */
	int bufferSize = 0;
	for (int i = 0; i < argc; i++)
		bufferSize += 1 + strlen(argv[i]); /* + 1 for space & 
						      terminating null */
	char *commandLine = calloc(bufferSize, sizeof(*commandLine));
	int cumulativeLength = 0;
	for (int i = 0; i < argc; i++) {
		int n = strlen(argv[i]);
		memcpy(&commandLine[cumulativeLength], argv[i], n);
		cumulativeLength += n;
		commandLine[cumulativeLength] = ' ';
		cumulativeLength += 1;
	}
	/* change last space to a null byte */
	commandLine[cumulativeLength - 1] = '\0';


	/* Measurement header */
	char *measHeaderStrings[2];
	measHeaderStrings[0] = getWorldInfo();
	measHeaderStrings[1] = integratorInfo(&integratorConf);
	char *measHeader = asprintfOrDie("# %s\n# %s\n%s%s"
				"# random number generator seed %llx\n",
				VERSION, commandLine,
				measHeaderStrings[0], measHeaderStrings[1],
				seed, commandLine);
	free(measHeaderStrings[0]); free(measHeaderStrings[1]);
	free(commandLine);
	/* Add extra header if we are using an end-to-end interaction */
	if (harmonicEndToEndInteraction) {
		char *eteHeader = harmonicEndToEndIntHeader(&hetei);
		char *prevMeasHeader = measHeader;
		measHeader = asprintfOrDie("%s%s",
				prevMeasHeader,
				eteHeader);
		free(prevMeasHeader);
	}
	if (endToEndForceInteraction) {
		char *eteHeader = endToEndForceIntHeader(&etefi);
		char *prevMeasHeader = measHeader;
		measHeader = asprintfOrDie("%s%s",
				prevMeasHeader,
				eteHeader);
		free(prevMeasHeader);
	}
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

	/* TODO: do the things below cleaner with less duplication (this is 
	 * just a quick hack to test things atm) */
	/* Distance base pairing task */
	char *distanceBasePairFile = asprintfOrDie("%s%s", filenameBase,
						DISTANCE_BASE_PAIRING_FILE_SUFFIX);
	Measurement distanceBasePairing;
	bpc.mode = BPSM_DISTANCE;
	distanceBasePairing.sampler = basePairingSampler(&bpc);
	distanceBasePairing.measConf = additionalMeasConf; /* struct copy */
	distanceBasePairing.measConf.measureFile = distanceBasePairFile;
	Task distanceBasePairingTask = measurementTask(&distanceBasePairing);

	/* Energy base pairing task */
	char *energyBasePairFile = asprintfOrDie("%s%s", filenameBase,
						ENERGY_BASE_PAIRING_FILE_SUFFIX);
	Measurement energyBasePairing;
	bpc.mode = BPSM_ENERGY;
	energyBasePairing.sampler = basePairingSampler(&bpc);
	energyBasePairing.measConf = additionalMeasConf; /* struct copy */
	energyBasePairing.measConf.measureFile = energyBasePairFile;
	Task energyBasePairingTask = measurementTask(&energyBasePairing);

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

	/* Gyration radius task */
	char *gyradFile = asprintfOrDie("%s%s", filenameBase,
						GYRATION_FILE_SUFFIX);
	Measurement gyradMeas;
	gyradMeas.sampler = gyrationRadiusSampler(&world.strands[0]);
	gyradMeas.measConf = additionalMeasConf; /* struct copy */
	gyradMeas.measConf.measureFile = gyradFile;
	Task gyrationRadiusTask = measurementTask(&gyradMeas);

	/* Force, velocity, friction task */
	char *fvfFile = asprintfOrDie("%s%s", filenameBase,
						FORCE_VEL_FRIC_FILE_SUFFIX);
	Measurement fvfMeas;
	fvfMeas.sampler = forceVelFricSampler(&world.strands[0]);
	fvfMeas.measConf = additionalMeasConf; /* struct copy */
	fvfMeas.measConf.measureFile = fvfFile;
	Task forceVelFricTask = measurementTask(&fvfMeas);
	/* For the backbone phosphate */
	char *fvfPFile = asprintfOrDie("%s%s", filenameBase,
						FORCE_VEL_FRIC_P_FILE_SUFFIX);
	Measurement fvfPMeas;
	fvfPMeas.sampler = forceVelFricPSampler(&world.strands[0]);
	fvfPMeas.measConf = additionalMeasConf; /* struct copy */
	fvfPMeas.measConf.measureFile = fvfPFile;
	Task forceVelFricPTask = measurementTask(&fvfPMeas);
	/* For the phosphates relative to backbone */
	char *fvPFile = asprintfOrDie("%s%s", filenameBase,
						FORCE_VEL_P_FILE_SUFFIX);
	Measurement fvPMeas;
	fvPMeas.sampler = forceVelPSampler(&world.strands[0]);
	fvPMeas.measConf = additionalMeasConf; /* struct copy */
	fvPMeas.measConf.measureFile = fvPFile;
	Task forceVelPTask = measurementTask(&fvPMeas);
	/* For the sugars relative to backbone */
	char *fvSFile = asprintfOrDie("%s%s", filenameBase,
						FORCE_VEL_S_FILE_SUFFIX);
	Measurement fvSMeas;
	fvSMeas.sampler = forceVelSSampler(&world.strands[0]);
	fvSMeas.measConf = additionalMeasConf; /* struct copy */
	fvSMeas.measConf.measureFile = fvSFile;
	Task forceVelSTask = measurementTask(&fvSMeas);

	/* Hairpin config task */
	char *stateFile = asprintfOrDie("%s%s", filenameBase,
						STATE_SUFFIX);
	hcsc.configFilePrefix = stateFile;
	hcsc.s = &world.strands[0];
	Measurement stateMeas;
	stateMeas.sampler = hairpinConfigSampler(&hcsc);
	stateMeas.measConf = additionalMeasConf; /* struct copy */
	stateMeas.measConf.measureFile = stateFile;
	Task stateTask = measurementTask(&stateMeas);

	/* Hairpin task */
	Measurement hairpin;
	char *zippedStateFile = NULL;
	switch (hairpinMeasurementType) {
	case HAIRPIN_MELTING_TEMPERATURE:
		hairpin.sampler = hairpinMeltingTempSampler(&hmtc);
		break;
	case HAIRPIN_FORMATION_TIME:
		if (hfc.requiredBoundBPs < 0 || hfc.nucleationBoundBPs < 0)
			die("You need to give the minimum required number "
					"of (nucleation) bound base pairs "
					"to measure hairpin formation time!\n");
		zippedStateFile = asprintfOrDie("%s%s", filenameBase,
						ZIPPED_STATE_FILE_SUFFIX);
		hfc.zippedStateFile = zippedStateFile;
		hairpin.sampler = hairpinFormationSampler(&hfc);
		break;
	case HAIRPIN_STATE:
		hairpin.sampler = hairpinStateSampler(&hsc);
		break;
	case HAIRPIN_NO_MEASUREMENT:
		hairpin.sampler = trivialSampler();
		break;
	default:
		die("Unknown hairpin measurement type!\n");
	}
	hairpin.measConf = measurementConf;
	Task hairpinTask = measurementTask(&hairpin);

	/* Render task */
	Task renderTask = makeRenderTask(&renderConf);

	//TODO
	if (useMonteCarlo)
		integratorTask = makeMonteCarloTask(&monteCarloConfig);

	/* Combined task */
	Task *tasks[15];
	tasks[0]  = (render ? &renderTask : NULL);
	tasks[1]  = (noDynamics ? NULL : &integratorTask);
	tasks[2]  = &verboseTask;
	tasks[3]  = &hairpinTask;;
	tasks[4]  = (measureEndToEndDistance ? &endToEndTask : NULL);
	tasks[5]  = (measureBasePairing ? &basePairingTask : NULL);
	tasks[6]  = (measureDistanceBasePairing ? &distanceBasePairingTask : NULL);
	tasks[7]  = (measureEnergyBasePairing ? &energyBasePairingTask : NULL);
	tasks[8]  = (measureTemperature ? &temperatureTask : NULL);
	tasks[9]  = (measureGyrationRadius ? &gyrationRadiusTask : NULL);
	tasks[10] = (measureForceVelFric ? &forceVelFricTask : NULL);
	tasks[11] = (measureForceVelFric ? &forceVelFricPTask : NULL);
	tasks[12] = (measureForceVelFric ? &forceVelPTask : NULL);
	tasks[13] = (measureForceVelFric ? &forceVelSTask : NULL);
	tasks[14] = (hairpinMeasurementType == HAIRPIN_STATE ? &stateTask : NULL);
	Task task = sequence(tasks, 15);

	registerInteractions(interactionSettings);

	bool everythingOK = run(&task);

	if (finalStateFile != NULL)
		writeWorld(finalStateFile);

	freeWorld();
	free(basePairFile);
	free(distanceBasePairFile);
	free(energyBasePairFile);
	free(endToEndFile);
	free(temperatureFile);
	free(gyradFile);
	free(fvfFile);
	free(fvfPFile);
	free(fvPFile);
	free(fvSFile);
	free(stateFile);
	free(zippedStateFile);
	free(measHeader);

	if (!everythingOK)
		return 1;
	return 0;
}


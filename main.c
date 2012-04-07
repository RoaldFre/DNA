#define _GNU_SOURCE

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include "main.h"
#include "vmath.h"
#include "physics.h"
#include "render.h"
#include "samplers.h"

#define DATA_FILE_NAME "/tmp/data.txt"

/* Defaults */
#define DEF_TIMESTEP 			1.0
#define DEF_TEMPERATURE 		300.0
#define DEF_LANGEVIN_GAMMA		5e12 //TODO sane?
#define DEF_COUPLING_TIMESTEP_FACTOR 	1000
#define DEF_TRUNCATION_LENGTH		20.0
#define DEF_MONOMER_WORLDSIZE_FACTOR    8.0
#define DEF_MONOMERS_PER_RENDER 	2000
#define DEF_RENDER_RADIUS 		1.5
#define DEF_MEASUREMENT_WAIT 		4e4
#define DEF_RENDER_FRAMERATE 		30.0
#define DEF_INTEGRATOR			LANGEVIN

static void printUsage(void)
{
	printf("Usage: main <number of monomers> [flags]\n");
	printf("\n");
	printf("Flags:\n");
	printf(" -t <flt>  length of Time steps (in femtoseconds)\n");
	printf("             default: %f\n", DEF_TIMESTEP);
	printf(" -T <flt>  Temperature\n");
	printf("             default: %f\n", DEF_TEMPERATURE);
	printf(" -r        Render\n");
	printf(" -f <flt>  desired Framerate when rendering.\n");
	printf("             default: %f)\n", DEF_RENDER_FRAMERATE);
	printf(" -R <flt>  Radius (in Angstrom) of the particles when rendering\n");
	printf("             default: %f\n", DEF_RENDER_RADIUS);
	printf(" -l <flt>  truncation Length of potentials (in Angstrom).\n");
	printf("             negative value: sets truncation to worldsize/2\n");
	printf("             default: %f\n", DEF_TRUNCATION_LENGTH);
	printf(" -S <flt>  Size of world (in Angstrom).\n");
	printf("             default: (number of monomers) * %f\n", DEF_MONOMER_WORLDSIZE_FACTOR);
	printf(" -b <num>  number of Boxes per dimension\n");
	printf("             default: max so that boxsize >= potential truncation length\n");
	printf(" -v <int>  Verbose: dump statistics every <int> iterations\n");
	printf(" -E <flt>  dump Energy statistics every <flt> femtoseconds.\n");
	printf("           Don't forget to set the number of samples with '-s'!\n");
	printf(" -s <int>  accumulate <int> measurement Samples\n");
	printf("             default: Don't sample, loop forever\n");
	printf(" -w <flt>  Wait <flt> femtoseconds before starting the measurements\n");
	printf("             default: %f\n", DEF_MEASUREMENT_WAIT);
	printf(" -i <type> Integrator to use. Values for <type>:\n");
	printf("             l: Langevin (velocity BBK) [default]\n");
	printf("             v: velocity Verlet with Berendsen thermostat\n");
	printf("\n");
	printf("Parameters for velocity Verlet + Berendsen termostat:\n");
	printf(" -g <flt>  Gamma: friction coefficient for Langevin dynamics\n");
	printf("             default: %e\n", DEF_LANGEVIN_GAMMA);
	printf("\n");
	printf("Parameters for velocity Verlet + Berendsen termostat:\n");
	printf(" -c <flt>  thermal bath Coupling: relaxation time (zero to disable)\n");
	printf("             default: %d * timestep\n", DEF_COUPLING_TIMESTEP_FACTOR);
}

static void parseArguments(int argc, char **argv)
{
	int c;

	/* defaults */
	config.verbose          = 0;
	config.measureSamples   = -1; /* loop indefinitely */
	config.measureWait	= DEF_MEASUREMENT_WAIT * TIME_FACTOR;
	config.timeStep 	= DEF_TIMESTEP * TIME_FACTOR;
	config.thermostatTemp	= DEF_TEMPERATURE;
	config.radius		= DEF_RENDER_RADIUS * LENGTH_FACTOR;
	config.truncationLen    = DEF_TRUNCATION_LENGTH;
	config.framerate        = DEF_RENDER_FRAMERATE;
	config.langevinGamma	= DEF_LANGEVIN_GAMMA;
	config.integrator	= DEF_INTEGRATOR;

	/* guards */
	config.worldSize = -1;
	config.thermostatTau = -1;
	config.measureInterval = -1;
	config.numBoxes = -1;

	while ((c = getopt(argc, argv, ":t:T:g:c:f:rR:l:S:b:v:E:s:w:i:h")) != -1)
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
			config.framerate = atof(optarg);
			if (config.framerate < 0)
				die("Invalid framerate %s\n", optarg);
			break;
		case 'r':
			config.render = true;
			break;
		case 'R':
			config.radius = atof(optarg) * LENGTH_FACTOR;
			if (config.radius <= 0)
				die("Invalid radius %s\n", optarg);
			break;
		case 'l':
			config.truncationLen = atof(optarg) * LENGTH_FACTOR;
			break;
		case 'S':
			config.worldSize = atof(optarg) * 1e-10;
			if (config.worldSize <= 0)
				die("Invalid world size %s\n", optarg);
			break;
		case 'b':
			config.numBoxes = atoi(optarg);
			if (config.numBoxes <= 0)
				die("Invalid number of boxes %s\n",
						optarg);
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
		case 'w':
			config.measureWait = atof(optarg) * TIME_FACTOR;
			if (config.measureWait < 0)
				die("Invalid wait time %f\n", config.measureWait);
			break;
		case 'i':
			if (optarg[0] == '\0' || optarg[1] != '\0')
				die("Integrator: badly formatted integrator type\n");
			switch(optarg[0]) {
				case 'l': config.integrator = LANGEVIN; break;
				case 'v': config.integrator = VERLET; break;
				default: die("Unknown integrator type%s\n", optarg);
					 break;
			}
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

	if (config.truncationLen < 0) {
		/* Disable truncation -> no space partitioning */
		config.numBoxes = 1;
		config.truncationLen = config.worldSize / 2.0;
		/* Due to (cubic) periodicity, we still need to truncate at 
		 * worldsize/2 to have correct energy conservation and to 
		 * make sure every particle 'sees' the same (spherical) 
		 * potential, no matter where it is within the cube. */
	} else if (config.truncationLen > config.worldSize / 2.0)
		config.truncationLen = config.worldSize / 2.0; /* same reason */

	if (config.numBoxes == -1) {
		config.numBoxes = config.worldSize / config.truncationLen;
		if (config.numBoxes  < 1)
			config.numBoxes = 1;
	}

	if (config.worldSize / config.numBoxes < config.truncationLen)
		die("The boxsize (%f) is smaller than the potential "
			"truncation radius (%f)!\n",
			config.worldSize / config.numBoxes,
			config.truncationLen);
}

void die(const char *fmt, ...)
{
	va_list args;

	va_start(args, fmt);
	vfprintf(stderr, fmt, args);
	va_end(args);

	exit(1);
}



typedef struct timer {
	struct timeval prev; /* time at last invocation of tickTimer */
	double interval;
	double accum; /* accumulated time since last tick */
} Timer;

static Timer makeTimer(double interval)
{
	Timer timer;
	gettimeofday(&timer.prev, NULL);
	timer.interval = interval;
	timer.accum = 0;
	return timer;
}

static bool tickTimer(Timer *timer)
{
	time_t prev_sec       = timer->prev.tv_sec;
	suseconds_t prev_usec = timer->prev.tv_usec;
	gettimeofday(&timer->prev, NULL);
	time_t sec       = timer->prev.tv_sec;
	suseconds_t usec = timer->prev.tv_usec;

	timer->accum += sec - prev_sec + ((double) (usec - prev_usec)) / 1e6;
	if (timer->accum < timer->interval)
		return false;

	timer->accum = 0;
	//timer->accum -= timer->interval;
	//timer->accum = fmod(timer->accum, timer->interval);
	return true;
}


/* Advance the simulation by one time step. Render and/or dump statistics 
 * if neccesary. Return false if the user wants to quit. */
static bool stepSimulation(Timer *renderTimer) {
	static int  stepsSinceVerbose = 0;

	stepWorld();

	if (config.verbose > 0) {
		stepsSinceVerbose++;
		if (stepsSinceVerbose > config.verbose) {
			stepsSinceVerbose = 0;
			dumpStats();
		}
	}

	if (config.render && tickTimer(renderTimer))
			return stepGraphics();

	return true;
}

int main(int argc, char **argv)
{
	srand(time(NULL)); //seed random generator

	parseArguments(argc, argv);

	allocWorld(1, config.numBoxes, config.worldSize); // TODO
	allocStrand(&world.strands[0], config.numMonomers);
	//allocStrand(&world.strands[1], config.numMonomers);
	fillWorld();


	/* TODO this is just a quick and dirty test of the task and 
	 * measurement framework for now! */
	MeasurementConf measConf;
	measConf.measureSamples = 10;
	measConf.measureInterval = 10e-15;
	measConf.measureWait = 200e-15;

	Measurement meas;
	meas.measConf = measConf;
	meas.sampler = averageTemperatureSampler();

	Task measTask = measurementTask(&meas);

	Task *tasks[3];
	tasks[0] = &renderTask;
	tasks[1] = &integratorTask;
	tasks[2] = &measTask;
	Task task = sequence(tasks, 3);
	run(&task);
	return 0;

#if 0

	Timer renderTimer = makeTimer(1.0 / config.framerate);

	if (config.render)
		initRender();

	if (config.measureSamples < 0) {
		/* Loop forever, or until the user quits the renderer */
		while (stepSimulation(&renderTimer));
	} else {
		printf("Waiting for system to relax.\n");
		for (double t = 0; keepGoing && t < config.measureWait; t += config.timeStep) {
			keepGoing = stepSimulation(&renderTimer);
			if (fmod(t, config.measureWait / 100) < config.timeStep) {
				printf("\rRelax time %13f of %f",
						(t + config.measureWait/100) / TIME_FACTOR, 
						config.measureWait / TIME_FACTOR);
				fflush(stdout);
			}
		}

		/* Perform the measurements */
		printf("\nStarting measurement.\n");
		FILE *outstream = fopen(DATA_FILE_NAME, "w");
		double intervalTime = 0;
		for (long sample = 0; keepGoing && sample < config.measureSamples; sample++) {
			while (keepGoing && intervalTime <= config.measureInterval) {
				keepGoing = stepSimulation(&renderTimer);
				intervalTime += config.timeStep;
			}
			if (!keepGoing)
				break;

			/* Check for numerical drift (or bugs) before 
			 * commiting measurement. */
			if (!physicsCheck())
				die("You broke physics!\n");

			dumpEnergies(outstream);
			printf("\rMeasured sample %ld/%ld", sample + 1, config.measureSamples);
			fflush(stdout);
			intervalTime -= config.measureInterval;
		}
		printf("\n");
		fclose(outstream);
	}

	freeWorld();
	return 0;
#endif
}

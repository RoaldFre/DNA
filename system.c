#include "system.h"
#include "task.h"

World world;
Config config;

static double sim_time = 0;
static long iteration = 0;

double getTime(void)
{
	return sim_time;
}

long getIteration(void)
{
	return iteration;
}

void run(Task *task)
{
	void *state = taskStart(task);
	bool keepGoing = true;
	while (keepGoing) {
		keepGoing = taskTick(task, state);
		sim_time += config.timeStep;
		iteration++;
	}
	taskStop(task, state);
}




#include "system.h"
#include "task.h"

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

/* Returns true if everything went according to plan. False if something 
 * unexpected happened. */
bool run(Task *task)
{
	void *state = taskStart(task);
	TaskSignal taskSig = TASK_OK;
	while (taskSig == TASK_OK) {
		taskSig = taskTick(task, state);
		sim_time += config.timeStep;
		iteration++;
	}
	taskStop(task, state);

	return taskSig != TASK_ERROR;
}




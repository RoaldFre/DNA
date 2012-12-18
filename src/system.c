#include <stdarg.h>
#include "system.h"
#include "task.h"

static long iteration = 0;

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
		iteration++;
	}
	taskStop(task, state);

	return taskSig != TASK_ERROR;
}

void die(const char *fmt, ...)
{
	va_list args;

	va_start(args, fmt);
	vfprintf(stderr, fmt, args);
	va_end(args);

	assert(false); /* Debugging breakpoint */
	exit(1);
}
void dieMem()
{
	die("Eek! Looks like we are out of memory! Bailing out!\n");
}

char *asprintfOrDie(const char *fmt, ...)
{
	va_list args;
	char *ret;

	va_start(args, fmt);
	int err = vasprintf(&ret, fmt, args);
	va_end(args);

	if (err < 0)
		dieMem();

	return ret;
}

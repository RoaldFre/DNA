#include "system.h"
#include "task.h"

#include <stdlib.h>
#include <string.h>

void *seqStart(void *initialData);
bool  seqTick(void *state);
void  seqStop(void *state);


/* Returns the state pointer that gets returned from task.start(), or NULL 
 * if task.start == NULL */
void *taskStart(Task *task)
{
	if (task->start == NULL)
		return NULL;
	return task->start(task->initialData);
}

/* Tick the task and return task.tick(), or do nothing if task.tick == NULL 
 * and return true in that case */
bool taskTick(Task *task, void *state)
{
	if (task->tick == NULL)
		return true;
	return task->tick(state);
}

/* Stop the task, or do nothing if task.stop == NULL */
void taskStop(Task *task, void *state)
{
	if (task->stop == NULL)
		return;
	task->stop(state);
}




typedef struct seqData
{
	int num;
	Task *tasks;
} SeqData;

Task sequence(Task **tasks, int num)
{
	/* We must alloc on the heap because the memory must remain 
	 * resident. Free it in seqStart. */
	SeqData *seqData = malloc(sizeof(*seqData));

	int trueNum = 0;
	/* Find out how many 'actual' (non-NULL) tasks there are */
	for (int i = 0; i < num; i++)
		if(tasks[i] != NULL)
			trueNum++;

	/* We must copy, because the list of pointers we received may point 
	 * to the stack of the caller, and may become invalid by the time 
	 * we actually run the simulation. */
	seqData->num = trueNum;
	seqData->tasks = calloc(trueNum, sizeof(*seqData->tasks));
	int j = 0;
	for (int i = 0; i < num; i++) {
		if(tasks[i] == NULL)
			continue;
		memcpy(&seqData->tasks[j], tasks[i], sizeof(*seqData->tasks));
		j++;
	}
	
	Task seq;
	seq.initialData = seqData;
	seq.start = &seqStart;
	seq.tick  = &seqTick;
	seq.stop  = &seqStop;

	return seq;
}

typedef struct seqState
{
	int num;
	Task *tasks;
	void **states; /* This will hold the states for all tasks */
} SeqState;

void *seqStart(void *initialData)
{
	SeqData *seqData = (SeqData*) initialData;

	SeqState *seqState = calloc(seqData->num, sizeof(*seqState));
	seqState->tasks    = seqData->tasks;
	seqState->states   = calloc(seqData->num, sizeof(seqState->states));
	seqState->num      = seqData->num;

	for (int i = 0; i < seqData->num; i++)
		seqState->states[i] = taskStart(&seqData->tasks[i]);

	/* We can free the initial data struct now */
	free(seqData);
	
	return (void*) seqState;
}

bool seqTick(void *state)
{
	SeqState *seqState = (SeqState*) state;

	for (int i = 0; i < seqState->num; i++) {
		Task *task = &seqState->tasks[i];
		if (!taskTick(task, seqState->states[i]))
			return false;
	}
	return true;
}

void seqStop(void *state)
{
	SeqState *seqState = (SeqState*) state;

	for (int i = 0; i < seqState->num; i++) {
		Task *task = &seqState->tasks[i];
		taskStop(task, seqState->states[i]);
	}

	/* We can free the state now */
	free(seqState);
}


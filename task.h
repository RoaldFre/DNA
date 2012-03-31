#ifndef _TASK_H_
#define _TASK_H_

/* Some functions for manipulating and combining tasks */

#include "system.h"

#define UNUSED(x) ((void) x)

/* Generate a new Task that executes all given tasks in sequence. The list 
 * of tasks will be copied to an internal structure, so it is safe to free 
 * the data after calling this function.
 *
 * args: sequence: a list to task pointers, all of which should be valid! */
Task sequence(Task **tasks, int num);


/* Some helper/wrapper functions */

/* Returns the state pointer that gets returned from task->start(), or NULL 
 * if task->start == NULL */
void *taskStart(Task *task);
/* Tick the task and return task->tick(), or do nothing if task->tick == NULL
 * and return true in that case */
bool taskTick(Task *task, void *state);
/* Stop the task, or do nothing if task->stop == NULL */
void taskStop(Task *task, void *state);

#endif


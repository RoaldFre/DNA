#ifndef _RENDER_H_
#define _RENDER_H_

#include "task.h"

int initRender(void);
bool stepGraphics(void);

extern Task renderTask;

#endif

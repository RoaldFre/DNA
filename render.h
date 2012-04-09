#ifndef _RENDER_H_
#define _RENDER_H_

#include "task.h"

typedef struct
{
	int framerate; /* The desired framerate */
	double radius; /* The radius of the particles to render */
} RenderConf;

Task makeRenderTask(RenderConf *rc);

#endif

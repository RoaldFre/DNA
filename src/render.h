#ifndef _RENDER_H_
#define _RENDER_H_

#include "task.h"

typedef struct
{
	int framerate; /* The desired framerate */
	double radius; /* The radius of the particles to render */
	bool drawForces; /* Draw the forces on the particles? */
} RenderConf;

Task makeRenderTask(RenderConf *rc);

typedef struct {
	const char *string;
	int x, y; /* Pixel-coordinates of the string to render */
} RenderStringConfig;

/* Register the given string to be rendered at the specified position */
void registerString(RenderStringConfig *rsc);

#endif

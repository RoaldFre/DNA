#ifndef _RENDER_H_
#define _RENDER_H_

#include "task.h"

/* NOTE:
 * You can disable rendering by building with a NO_RENDER define (eg by 
 * passing -DNO_RENDER to the compiler), or by uncommenting the line below. */
// #define NO_RENDER 

typedef struct
{
	int framerate; /* The desired framerate */
	double radius; /* The radius of the particles to render */
	bool drawForces; /* Draw the forces on the particles? */
	int numBoxes;  /* Number of boxes (in each dimension) to render 
			  when pressing 'b'. */
} RenderConf;

Task makeRenderTask(RenderConf *rc);

typedef struct {
	const char *string;
	int x, y; /* Pixel-coordinates of the string to render */
} RenderStringConfig;

/* Register the given string to be rendered at the specified position */
void registerString(RenderStringConfig *rsc);

#endif

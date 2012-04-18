#include <stdio.h>
#include <SDL/SDL.h>
#include <GL/gl.h>
#include <math.h>
#include "main.h"
#include "render.h"
#include "physics.h"
#include "task.h"

#define SCREEN_W 1000
#define SCREEN_H 1000

#define SPHERE_SLICES 10
#define CILINDER_FACES 6
#define CILINDER_RADIUS_DIVISOR 10

typedef struct {
	GLfloat x, y, z;
} Vertex3;

const GLfloat light_pos[]  = {2.0, 1.0, 2.0, 0.0};
const GLfloat light_diff[] = {1.0, 1.0, 1.0, 0.0};
const GLfloat light_spec[] = {1.0, 0.0, 0.0, 0.0};
const GLfloat light_ambi[] = {0.8, 0.8, 0.8, 0.0};

const GLfloat red[]   = {1.0, 0.0, 0.0, 0.0};
const GLfloat green[] = {0.0, 1.0, 0.0, 0.0};
const GLfloat blue[]  = {0.0, 0.0, 1.0, 0.0};
const GLfloat gray[]  = {0.2, 0.2, 0.2, 0.0};

const GLfloat S_col[] = {0.1, 0.1, 0.1, 0.0};
const GLfloat P_col[] = {0.1, 0.1, 0.1, 0.0};
const GLfloat A_col[] = {0.1, 0.1, 1.0, 0.0};
const GLfloat T_col[] = {0.0, 0.7, 1.0, 0.0};
const GLfloat C_col[] = {0.0, 0.8, 0.0, 0.0};
const GLfloat G_col[] = {0.7, 1.0, 0.0, 0.0};

static int numVertices;
static int numIndices;
static Vertex3 *sphereVertex;
static GLushort *sphereIndex;

static SDL_Surface *surface;
static GLfloat view_angle;

static void createSphere(int slices, int *numVert, Vertex3 **vertices, int *numInd,
		GLushort **indices);
static void drawLine(Vec3 *p1, Vec3 *p2);
static void drawCilinder(Vec3 *p1, Vec3 *p2, int faces, double radius);
static void renderParticles(int num, Particle *ps, RenderConf *rc);
static void renderBase(Particle *p, RenderConf *rc);
static void renderStrand(Strand *s, RenderConf *rc);
static void renderConnection(Particle *p1, Particle *p2, RenderConf *rc);
static void gluPerspective(GLfloat fovy, GLfloat aspect, GLfloat zNear, 
		GLfloat zFar);

static void initRender(void);
static void calcFps(void);
void renderIfWeMust(RenderConf *rc);
static void render(RenderConf *rc);

void *renderTaskStart(void *initialData);
bool renderTaskTick(void *state);
void renderTaskStop(void *state);

static void calcFps()
{
	static long tock = 0;
	static int frames = 0;
	long tick;
	char string[32];

	frames++;
	tick = SDL_GetTicks(); /* mili seconds */

	if (tick - tock > 1000) {
		tock = tick;
		sprintf(string, "%u FPS", frames);
		SDL_WM_SetCaption(string, string);
		frames = 0;
	}
	return;
}


/* Handle events.
 * Return true if everything went fine, false if user requested to quit. */
static bool handleEvents(void)
{
	SDL_Event event;
	while (SDL_PollEvent(&event)) {
		if (event.type == SDL_QUIT) {
			printf("\nRequested Quit.\n\n");
			return false;
		}
		else if (event.type == SDL_KEYDOWN) {
			switch (event.key.keysym.sym) {
			case SDLK_ESCAPE:
				printf("\nRequested Quit.\n\n");
				return false;
				break;
			case SDLK_LEFT:
				view_angle--;
				break;
			case SDLK_RIGHT:
				view_angle++;
				break;
			case SDLK_UP:
				config.timeStep *= 1.1;
				printf("Time step: %f\n", 
						config.timeStep / TIME_FACTOR);
				break;
			case SDLK_DOWN:
				config.timeStep /= 1.1;
				printf("Time step: %f\n",
						config.timeStep / TIME_FACTOR);
				break;
			case SDLK_SPACE:
				config.thermostatTemp *= 1.1;
				printf("Temperature: %f\n", config.thermostatTemp);
				break;
			case SDLK_BACKSPACE:
				config.thermostatTemp /= 1.1;
				printf("Temperature: %f\n", config.thermostatTemp);
				break;
			case SDLK_RETURN:
				SDL_WM_ToggleFullScreen(surface);
				break;
			default:
				break;
			}
		}
	}
	return true;
}


/* Render the image if we must, based on the requested framerate and the 
 * last invocation of this function. */
void renderIfWeMust(RenderConf *rc)
{
	static long tock = -1000; /* always draw first frame immediately */

	if (rc->framerate <= 0) {
		render(rc);
		return;
	}

	long tick = SDL_GetTicks(); /* mili seconds */

	if (tick - tock > 1000 / rc->framerate) {
		tock = tick;
		render(rc);
	}
}

static void gluPerspective(GLfloat fovy, GLfloat aspect, GLfloat zNear, 
		GLfloat zFar)
{
	GLfloat xMin, xMax, yMin, yMax;

	yMax = zNear * tan(fovy * M_PI / 360.0);
	yMin = -yMax;

	xMin = yMin * aspect;
	xMax = yMax * aspect;

	glFrustum(xMin, xMax, yMin, yMax, zNear, zFar);
}

/* Returns false if we couldn't initialize, true otherwise */
static void initRender(void)
{
	int flags = 0;
	double ws = world.worldSize;
	const SDL_VideoInfo *vidinfo;

	if (SDL_Init(SDL_INIT_VIDEO) < 0)
		die(SDL_GetError());
		//TODO cleaner way
	
	vidinfo = SDL_GetVideoInfo();
	if (vidinfo == NULL)
		die(SDL_GetError());
	
	flags |= SDL_OPENGL;
	flags |= SDL_HWPALETTE;
	flags |= (vidinfo->hw_available ? SDL_HWSURFACE : SDL_SWSURFACE);
	if (vidinfo->blit_hw) flags |= SDL_HWACCEL;

	SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);
	SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 24);
	SDL_GL_SetAttribute(SDL_GL_STENCIL_SIZE, 0);
	SDL_GL_SetAttribute(SDL_GL_SWAP_CONTROL, 1);

	surface = SDL_SetVideoMode(SCREEN_W, SCREEN_H, 24, flags);
	if (surface == NULL)
		die(SDL_GetError());

	SDL_EnableKeyRepeat(1, 
			SDL_DEFAULT_REPEAT_INTERVAL);

	/*SDL_WM_ToggleFullScreen(surface);	*/

	//atexit(SDL_Quit); //Do it when we stop our render task.

	/* OpenGL Init */
	glEnable(GL_DEPTH_TEST);
	glShadeModel(GL_SMOOTH);
	glEnable(GL_NORMALIZE);

	glEnable(GL_LIGHTING);

	glEnable(GL_LIGHT0);
	glLightfv(GL_LIGHT0, GL_POSITION, light_pos);
	glLightfv(GL_LIGHT0, GL_DIFFUSE,  light_diff);
	glLightfv(GL_LIGHT0, GL_SPECULAR, light_spec);
	glLightfv(GL_LIGHT0, GL_AMBIENT,  light_ambi);

	glClearColor(1.0, 1.0, 1.0, 0.0);

	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_NORMAL_ARRAY);

	createSphere(SPHERE_SLICES, &numVertices, &sphereVertex,
			&numIndices, &sphereIndex);
	glVertexPointer(3, GL_FLOAT, sizeof(Vertex3), sphereVertex);
	glNormalPointer(   GL_FLOAT, sizeof(Vertex3), sphereVertex);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(35, SCREEN_W/(double)SCREEN_H, ws/1000, 100*ws);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

static void renderStrand(Strand *s, RenderConf *rc) {
	glLightfv(GL_LIGHT0, GL_DIFFUSE, S_col);
	renderParticles(s->numMonomers, s->Ss, rc);

	glLightfv(GL_LIGHT0, GL_DIFFUSE, P_col);
	renderParticles(s->numMonomers, s->Ps, rc);

	for (int i = 0; i < s->numMonomers; i++)
		renderBase(&s->Bs[i], rc);

	/* Connections */
	glLightfv(GL_LIGHT0, GL_DIFFUSE, gray);
	renderConnection(&s->Ps[0], &s->Ss[0], rc);
	renderConnection(&s->Ss[0], &s->Bs[0], rc);
	for(int i = 1; i < s->numMonomers; i++) {
		renderConnection(&s->Ps[i], &s->Ss[i],   rc);
		renderConnection(&s->Ss[i], &s->Bs[i],   rc);
		renderConnection(&s->Ss[i], &s->Ps[i-1], rc);
	}

	/* Forces */
	if (rc->drawForces) {
		glColor3f(0.8, 0.0, 0.0);
		glBegin(GL_LINES);
			for(int i = 0; i < 3 * s->numMonomers; i++) {
				Vec3 tmp;
				add(&s->all[i].pos, &s->all[i].F, &tmp);
				drawLine(&s->all[i].pos, &tmp);
			}
		glEnd();
	}
}

/* Renders the frame and calls calcFps() */
static void render(RenderConf *rc)
{
	double ws = world.worldSize;

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

/*	view_angle += 0.01;*/

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glTranslatef(0, 0, -ws*2.5);

	glRotatef(view_angle, 0, 1, 0);

	glColor3f(0.0, 1.0, 0.0);
	glBegin(GL_LINE_LOOP);
		glVertex3f(-ws/2, -ws/2, -ws/2);
		glVertex3f(-ws/2, -ws/2, +ws/2);
		glVertex3f(-ws/2, +ws/2, +ws/2);
		glVertex3f(-ws/2, +ws/2, -ws/2);
	glEnd();

	glBegin(GL_LINE_LOOP);
		glVertex3f(+ws/2, -ws/2, -ws/2);
		glVertex3f(+ws/2, -ws/2, +ws/2);
		glVertex3f(+ws/2, +ws/2, +ws/2);
		glVertex3f(+ws/2, +ws/2, -ws/2);
	glEnd();

	glTranslatef(-ws/2, -ws/2, -ws/2);

	for (int s = 0; s < world.numStrands; s++)
		renderStrand(&world.strands[s], rc);

	SDL_GL_SwapBuffers();

	calcFps();
}

static void drawPoint(Vec3 *p)
{
	glVertex3f(p->x, p->y, p->z);
}

static void drawLine(Vec3 *p1, Vec3 *p2)
{
	drawPoint(p1);
	drawPoint(p2);
}

static void drawCilinder(Vec3 *p1, Vec3 *p2, int faces, double radius)
{
	Vec3 d, dnorm;
	sub(p2, p1, &d);
	normalize(&d, &dnorm);
	Vec3 e1 = {1, 0, 0};
	Vec3 e2 = {0, 1, 0};
	/* Make orthogonal basis for cilinders */
	Vec3 b1, b2;
	if (dot(&e1, &dnorm) < 0.8)
		b1 = cross(&e1, &dnorm);
	else
		b1 = cross(&e2, &dnorm);
	normalize(&b1, &b1);
	b2 = cross(&b1, &dnorm);

	assert(fabs(length(&b1) - 1) < 1e-13);
	assert(fabs(length(&b2) - 1) < 1e-13);

	scale(&b1, radius, &b1);
	scale(&b2, radius, &b2);

	glBegin(GL_QUAD_STRIP);
		Vec3 l, r, u, v, dir;
		add(p1, &b1, &l);
		add(p2, &b1, &r);
		drawLine(&l, &r);
		for (int i = 1; i <= faces; i++) {
			double theta = 2*M_PI * i / faces;
			scale(&b1, cos(theta), &u);
			scale(&b2, sin(theta), &v);
			add(&u, &v, &dir);
			add(p1, &dir, &l);
			add(p2, &dir, &r);
			drawLine(&l, &r);
		}
	glEnd();
}

static void renderParticle(Particle *p, RenderConf *rc)
{
	glPushMatrix();
		glTranslatef(p->pos.x, p->pos.y, p->pos.z);
		glScalef(rc->radius, rc->radius, rc->radius);
		glDrawElements(GL_TRIANGLES, numIndices, GL_UNSIGNED_SHORT, sphereIndex);
	glPopMatrix();
}
static void renderParticles(int num, Particle *ps, RenderConf *rc)
{
	for (int i = 0; i < num; i++)
		renderParticle(&ps[i], rc);
}
static void renderBase(Particle *p, RenderConf *rc)
{
	const GLfloat *col;
	switch(p->type) {
		case BASE_A: col = A_col; break;
		case BASE_T: col = T_col; break;
		case BASE_C: col = C_col; break;
		case BASE_G: col = G_col; break;
		default: /* not a base */
		     assert(false);
		     col = gray;
	}
	glLightfv(GL_LIGHT0, GL_DIFFUSE, col);
	renderParticle(p, rc);
}
static void renderConnection(Particle *p1, Particle *p2, RenderConf *rc)
{
	if (distance(&p1->pos, &p2->pos) > world.worldSize / 2)
		return; /* To avoid periodic boundary clutter */
	drawCilinder(&p1->pos, &p2->pos, CILINDER_FACES,
			rc->radius / CILINDER_RADIUS_DIVISOR);
}

static void createSphere(int slices, int *numVert, Vertex3 **vertices, int *numInd,
		GLushort **indices)
{
	int i, j, k;
	double x, y, z;
	double r;
	int stacks;
	Vertex3 *vert;
	GLushort *ind;

	stacks = slices;
	slices *= 2;

	/* Plus two for the poles */
	*numVert = (stacks - 1) * slices + 2;
	*vertices = calloc(*numVert, sizeof(Vertex3));
	vert = *vertices;

	/* All but the top and bottom stack */
	for (i = 1; i < stacks; i++)
	{
		double phi = M_PI * i / (double) stacks - 2*M_PI;
		
		z = cos(phi);
		r = sqrt(1 - z*z);

		for (j = 0; j < slices; j++)
		{
			double theta = 2*M_PI*j/(double) slices;
			x = r * sin(theta);
			y = r * cos(theta);

			vert[(i-1) * slices + j + 1].x = x;
			vert[(i-1) * slices + j + 1].y = y;
			vert[(i-1) * slices + j + 1].z = z;
		}
	}

	/* Top and bottom */
	vert[0].x = 0;
	vert[0].y = 0;
	vert[0].z = 1;

	vert[*numVert-1].x = 0;
	vert[*numVert-1].y = 0;
	vert[*numVert-1].z = -1;

	*numInd = (stacks - 1) * slices * 6;
	*indices = calloc(*numInd, sizeof(GLushort));
	ind = *indices;

	k = 0;

	for (i = 1; i < slices; i++)
	{
		ind[k++] = 0;
		ind[k++] = i;
		ind[k++] = i+1;
	}
	ind[k++] = 0;
	ind[k++] = 1;
	ind[k++] = slices;
	
	for (i = 0; i < slices - 1; i++)
	{
		ind[k++] = *numVert - 1;
		ind[k++] = (*numVert - 1 - slices) + i;
		ind[k++] = (*numVert - 1 - slices) + i + 1;
	}
	ind[k++] = *numVert - 1;
	ind[k++] = *numVert - 1 - 1;
	ind[k++] = *numVert - 1 - slices + 0;

	for (i = 1; i < stacks - 1; i++)
	{
		int base = 1 + (i - 1) * slices;

		for (j = 0; j < slices - 1; j++)
		{
			ind[k++] = base + j;
			ind[k++] = base + slices + j;
			ind[k++] = base + slices + j + 1;

			ind[k++] = base + j;
			ind[k++] = base + j + 1;
			ind[k++] = base + slices + j + 1;
		}

		ind[k++] = base;
		ind[k++] = base + slices - 1;
		ind[k++] = base + slices;

		ind[k++] = base + slices - 1;
		ind[k++] = base + slices;
		ind[k++] = base + slices + slices - 1;
	}

	return;
}


void *renderTaskStart(void *initialData)
{
	assert(initialData != NULL);
	initRender();
	return initialData;
}

bool renderTaskTick(void *state)
{
	RenderConf *rc = (RenderConf*) state;

	bool ret = handleEvents();
	renderIfWeMust(rc);
	return ret;
}

void renderTaskStop(void *state)
{
	SDL_Quit();
	free(state);
}

Task makeRenderTask(RenderConf *rc)
{
	RenderConf *rcCopy = malloc(sizeof(*rcCopy));
	memcpy(rcCopy, rc, sizeof(*rcCopy));

	Task ret = {
		.initialData = rcCopy,
		.start = &renderTaskStart,
		.tick  = &renderTaskTick,
		.stop  = &renderTaskStop,
	};
	return ret;
}

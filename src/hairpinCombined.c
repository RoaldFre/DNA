/* This 'simulates' LTO by putting everything in one single source file. 
 * The result is sometimes faster, even if the usual compilation is done 
 * with LTO! */

#define NO_RENDER

#include "system.c"
#include "math.c"
#include "task.c"
#include "measure.c"
#include "samplers.c"
#include "physics.c"
#include "integrator.c"
#include "world.c"
#include "spgrid.c"
#include "tinymt/tinymt64.c"
#include "render.c"
#include "octave.c"
//#include "font.c"
//#include "mathlib/vector.c"
//#include "mathlib/quaternion.c"
//#include "mathlib/matrix.c"
#include "hairpin.c"

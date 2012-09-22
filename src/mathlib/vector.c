#include <math.h>

#include "mathlib.h"
#include "vector.h"

Vec3 vec3_lerp(Vec3 a, Vec3 b, double t)
{
	Vec3 c;
	/* c = a*t + b*(1-t) */
	a = scale(a, t);
	b = scale(b, 1-t);
	c = add(a, b);

	return c;
}

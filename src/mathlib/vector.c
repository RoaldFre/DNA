#include <math.h>

#include "mathlib.h"
#include "vector.h"

Vec3 vec3_add(Vec3 a, Vec3 b)
{
	Vec3 c;
	c.x = a.x + b.x;
	c.y = a.y + b.y;
	c.z = a.z + b.z;

	return c;
}

Vec3 vec3_sub(Vec3 a, Vec3 b)
{
	Vec3 c;
	c.x = a.x - b.x;
	c.y = a.y - b.y;
	c.z = a.z - b.z;

	return c;
}

double vec3_dot(Vec3 a, Vec3 b)
{
	return a.x*b.x + a.y*b.y + a.z*b.z;
}

double vec3_length(Vec3 a)
{
	return sqrt(vec3_dot(a, a));
}

Vec3 vec3_scale(Vec3 a, double lambda)
{
	Vec3 b;
	b.x = a.x * lambda;
	b.y = a.y * lambda;
	b.z = a.z * lambda;

	return b;
}

Vec3 vec3_normalize(Vec3 a)
{
	return vec3_scale(a, 1/vec3_length(a));
}

Vec3 vec3_cross(Vec3 a, Vec3 b)
{
	Vec3 c;
	c.x = a.y*b.z - a.z*b.y;
	c.y = a.z*b.x - a.x*b.z;
	c.z = a.x*b.y - a.y*b.x;

	return c;
}

Vec3 vec3_lerp(Vec3 a, Vec3 b, double t)
{
	Vec3 c;
	/* c = a*t + b*(1-t) */
	a = vec3_scale(a, t);
	b = vec3_scale(b, 1-t);
	c = vec3_add(a, b);

	return c;
}

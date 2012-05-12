#ifndef KOSMOS_VECTOR_H
#define KOSMOS_VECTOR_H

#include "../math.h" /* For Vec3 definition */

Vec3 vec3_add(Vec3 a, Vec3 b);
Vec3 vec3_sub(Vec3 a, Vec3 b);
Vec3 vec3_scale(Vec3 a, double r);
double vec3_dot(Vec3 a, Vec3 b);
double vec3_length2(Vec3 a);
double vec3_length(Vec3 a);
Vec3 vec3_normalize(Vec3 a);
Vec3 vec3_cross(Vec3 a, Vec3 b);
Vec3 vec3_lerp(Vec3 a, Vec3 b, double t);

#endif

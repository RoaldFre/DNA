#ifndef RENDER_QUATERNION_H
#define RENDER_QUATERNION_H

#include "matrix.h"
#include "vector.h"

typedef struct Quaternion {
	double w, x, y, z;
} Quaternion;

double quat_length2(Quaternion p);
double quat_length(Quaternion p);
double quat_dot(Quaternion a, Quaternion b);
Quaternion quat_add(Quaternion a, Quaternion b);
Quaternion quat_sub(Quaternion a, Quaternion b);
Quaternion quat_scale(Quaternion p, double scale);
Quaternion quat_normalize(Quaternion p);
Quaternion quat_conjugate(Quaternion p);
Quaternion quat_multiply(Quaternion p, Quaternion q);
Quaternion quat_from_angle_axis(double angle, double ax,
		double ay, double az);
Quaternion quat_trackball(int dx, int dy, double radius);
Quaternion quat_from_mat3(RenderMat3 m);
void mat3_from_quat(RenderMat3 m, Quaternion p);
Vec3 quat_transform(Quaternion q, Vec3 v);
Quaternion quat_nlerp(Quaternion a, Quaternion b, double t);
Quaternion quat_slerp(Quaternion a, Quaternion b, double t);
Quaternion quat_euler(double t1, double t2, double t3);

#endif

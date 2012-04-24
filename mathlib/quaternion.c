#include <math.h>

#include "mathlib.h"
#include "quaternion.h"

/* Square of the norm of the quaternion */
double quat_length2(Quaternion p)
{
	return quat_dot(p, p);
}

double quat_length(Quaternion p)
{
	return sqrt(quat_length2(p));
}

double quat_dot(Quaternion a, Quaternion b)
{
	return a.w*b.w + a.x*b.x + a.y*b.y + a.z*b.z;
}

Quaternion quat_add(Quaternion a, Quaternion b)
{
	Quaternion c;
	c.w = a.w + b.w;
	c.x = a.x + b.x;
	c.y = a.y + b.y;
	c.z = a.z + b.z;

	return c;
}

Quaternion quat_sub(Quaternion a, Quaternion b)
{
	Quaternion c;
	c.w = a.w - b.w;
	c.x = a.x - b.x;
	c.y = a.y - b.y;
	c.z = a.z - b.z;

	return c;
}

Quaternion quat_scale(Quaternion p, double lambda)
{
	Quaternion q;

	q.w = lambda * p.w;
	q.x = lambda * p.x;
	q.y = lambda * p.y;
	q.z = lambda * p.z;

	return q;
}

Quaternion quat_normalize(Quaternion p)
{
	return quat_scale(p, 1/quat_length(p));
}

Quaternion quat_conjugate(Quaternion p)
{
	Quaternion q;

	q.w =  p.w;
	q.x = -p.x;
	q.y = -p.y;
	q.z = -p.z;

	return q;
}

/* Multiply two quaternions: r = p * q */
Quaternion quat_multiply(Quaternion p, Quaternion q)
{
	Quaternion r;
	double a = p.w, b = p.x, c = p.y, d = p.z;
	double w = q.w, x = q.x, y = q.y, z = q.z;

	r.w = a*w - b*x - c*y - d*z;
	r.x = b*w + a*x + c*z - d*y;
	r.y = c*w + a*y + d*x - b*z;
	r.z = d*w + a*z + b*y - c*x;

	return r;
}

Quaternion quat_from_angle_axis(double alpha, double ax,
		double ay, double az)
{
	Quaternion q;

	double l2 = ax*ax + ay*ay + az*az;
	double l = sqrt(l2), s = sin(alpha/2), c = cos(alpha/2);

	/* If the axis is invalid, return the identity */
	if (l2 == 0.0)
	{
		return (Quaternion) {1, 0, 0, 0};
	}

	q.w = c;
	q.x = s*ax/l;
	q.y = s*ay/l;
	q.z = s*az/l;
	
	return q;
}

Quaternion quat_trackball(int dx, int dy, double radius)
{
	double dr, sina, cosa, sina2, cosa2;
	Quaternion q = {1, 0, 0, 0};

	if (dx == 0 && dy == 0)
		return q;

	dr = sqrt(dx*dx + dy*dy);

	sina = dr/radius;
	if (sina >= 1)
		sina = 0;
	cosa = sqrt(1 - sina*sina);

	cosa2 = sqrt((1 + cosa)/2);
	sina2 = sina/(2*cosa2);

	q.w = cosa2;
	q.x = -dy/dr * sina2;
	q.y =  dx/dr * sina2;
	q.z = 0;

	return q;
}

Quaternion quat_from_mat3(Mat3 m)
{
	Quaternion q;
	double T;
#define M(i, j) m[3*j + i]
	/* To understand this, research "convert orthogonal matrix to quaternion"
	 * This algorithm comes from "The Matrix and Quaternion FAQ" */
	T = 1 + M(0,0) + M(1,1) + M(2,2);
	if (T > 1e-3)
	{
		q.w = 0.5*sqrt(T);
		q.x = (M(2,1) - M(1,2)) / (4*q.w);
		q.y = (M(0,2) - M(2,0)) / (4*q.w);
		q.z = (M(1,0) - M(0,1)) / (4*q.w);
	} else
	{
		if (M(0,0) > M(1,1) && M(0,0) > M(2,2))
		{
			T = sqrt(1 + M(0,0) - M(1,1) - M(2,2));
			q.w = (M(2,1) - M(1,2)) / (2*T);
			q.x = 0.5*T;
			q.y = (M(0,1) + M(1,0)) / (2*T);
			q.z = (M(0,2) + M(2,0)) / (2*T);
		} else if (M(1,1) > M(2,2))
		{
			T = sqrt(1 - M(0,0) + M(1,1) - M(2,2));
			q.w = (M(0,2) - M(2,0)) / (2*T);
			q.x = (M(0,1) + M(1,0)) / (2*T);
			q.y = 0.5*T;
			q.z = (M(1,2) + M(2,1)) / (2*T);
		} else
		{
			T = sqrt(1 - M(0,0) - M(1,1) + M(2,2));
			q.w = (M(1,0) - M(0,1)) / (2*T);
			q.x = (M(0,2) + M(2,0)) / (2*T);
			q.y = (M(1,2) + M(2,1)) / (2*T);
			q.z = 0.5*T;
		}
	}
#undef M
	return q;
}

/* The Euler-Rodrigues formula */
void mat3_from_quat(Mat3 m, Quaternion p)
{
	double w = p.w, x = p.x, y = p.y, z = p.z;

#define M(i, j) m[3*j + i]
	M(0, 0) = w*w + x*x - y*y - z*z;
	M(0, 1) = 2*x*y - 2*w*z;
	M(0, 2) = 2*x*z + 2*w*y;

	M(1, 0) = 2*x*y + 2*w*z;
	M(1, 1) = w*w - x*x + y*y - z*z;
	M(1, 2) = 2*y*z - 2*w*x;

	M(2, 0) = 2*x*z - 2*w*y;
	M(2, 1) = 2*y*z + 2*w*x;
	M(2, 2) = w*w - x*x - y*y + z*z;
#undef M
}

Vec3 quat_transform(Quaternion q, Vec3 v)
{
	Mat3 m;

	mat3_from_quat(m, q);
	return mat3_transform(m, v);
}

Quaternion quat_nlerp(Quaternion a, Quaternion b, double t)
{
	Quaternion c;

	c = quat_add(quat_scale(a, 1-t), quat_scale(b, t));

	return quat_normalize(c);
}

Quaternion quat_slerp(Quaternion a, Quaternion b, double t)
{
	double cosa, alpha, interAngle;
	Quaternion q0, q1;

	cosa = quat_dot(a, b);
	if (cosa < 0)
	{
		a = quat_scale(a, -1);
		cosa = -cosa;
	}
	if (cosa > 1-1e-6)
		return quat_nlerp(a, b, t);

	alpha = acos(cosa);
	interAngle = alpha*t;

	q0 = a;
	/* q1 = b - (a.b) a */
	q1 = quat_normalize(quat_add(b, quat_scale(a,-cosa)));

	/* return q0 * cos(interAngle) + q1 * sin(interAngle) */
	return quat_add(quat_scale(q0, cos(interAngle)),
	                quat_scale(q1, sin(interAngle)));
}

Quaternion quat_euler(double t1, double t2, double t3)
{
	Mat3 mat;

	mat3_euler(t1, t2, t3, mat);

	return quat_from_mat3(mat);
}

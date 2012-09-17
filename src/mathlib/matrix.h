#ifndef RENDER_MATRIX_H
#define RENDER_MATRIX_H

#include "vector.h"

typedef double RenderMat3[9];

void mat3_mult(RenderMat3 a, RenderMat3 b, RenderMat3 c);
Vec3 mat3_transform(RenderMat3 m, Vec3 v);
void mat3_euler(double t1, double t2, double t3, RenderMat3 mat);

#endif

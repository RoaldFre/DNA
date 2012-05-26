#include <string.h>
#include <stdio.h>
#include <stdarg.h>
#include "octave.h"

void octaveStartComment(void)
{
	printf("## ");
}

void octaveEndComment(void)
{
	printf("\n");
}

void octaveComment(const char *fmt, ...)
{
	va_list args;

	octaveStartComment();
	va_start(args, fmt);
	vprintf(fmt, args);
	va_end(args);
	octaveEndComment();
}

void octaveScalar(const char *name, double value)
{
	printf("\n");
	printf("# name: %s\n", name);
	printf("# type: scalar\n");
	printf("%e\n", value);
}
void octaveString(const char *name, const char *string)
{
	printf("\n");
	printf("# name: %s\n", name);
	printf("# type: string\n");
	printf("# elements: 1\n");
	printf("# length: %d\n", (int) strlen(string));
	printf("%s\n", string);
}

void octaveMatrixHeader(const char *name, int rows, int cols)
{
	printf("\n");
	printf("# name: %s\n", name);
	printf("# type: matrix\n");
	printf("# rows: %d\n", rows);
	printf("# columns: %d\n", cols);
}

void octave3DMatrixHeader(const char *name, int nx, int ny, int nz)
{
	printf("\n");
	printf("# name: %s\n", name);
	printf("# type: matrix\n");
	printf("# ndims: 3\n");
	printf("%d %d %d\n", nx, ny, nz);
}


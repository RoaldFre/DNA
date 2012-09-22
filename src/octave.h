#include "math.h"

/* Don't use line breaks in comments. If you need multiple lines, call me 
 * multiple times.
 * Alternatively, call octaveStartComment(), then print stuff *without 
 * linebreaks* and end the comment by calling octaveEndComment().
 * Also. comments cannot be placed in between matrix data. */
void octaveComment(const char *fmt, ...);
void octaveStartComment(void);
void octaveEndComment(void);

void octaveScalar(const char *name, real value);
void octaveString(const char *name, const char *string);
void octaveMatrixHeader(const char *name, int rows, int cols);

/* After this, you print:
 *  - every column below the previous one in the 2D matrix name(:,:,i)
 *    [ie: dump the matrix in column major form]
 *  - do this for all such matrices i = 1:nz */
void octave3DMatrixHeader(const char *name, int nx, int ny, int nz);


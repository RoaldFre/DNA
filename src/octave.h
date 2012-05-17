/* Don't use line breaks in comments. If you need multiple lines, call me 
 * multiple times.
 * Also. comments cannot be placed in between matrix data. */
void octaveComment(const char *fmt, ...);
void octaveScalar(const char *name, double value);
void octaveString(const char *name, const char *string);
void octaveMatrixHeader(const char *name, int rows, int cols);

/* After this, you print:
 *  - every column below the nex one in the 2D matrix name(:,:,i)
 *  - do this for all such matrices i = 1:nz
 * Put only one value per line */
void octave3DMatrixHeader(const char *name, int nx, int ny, int nz);


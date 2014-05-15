#ifndef READ_DATA_1D_H_INCLUDED
#define READ_DATA_1D_H_INCLUDED

#include <stdio.h>

int read_data_1d(FILE *fp, int *ncols, int *nlines, double **data);
int read_data_1d_cols(FILE *fp, int ncols, int *cols, int *nlines, double **data);

#endif /* READ_DATA_1D_H_INCLUDED */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <line_reader.h>

#define STRTOK(STR,DELIM,SAVE) strtok_r(STR,DELIM,SAVE)

int read_data_1d(FILE *fp, int *ncols, int *nlines, double **data){
	struct line_reader lr;
	size_t len;
	char *line;
	int got_first_line = 0;
	*ncols = 0;
	*nlines = 0;
	lr_init(&lr, fp);
	int nlines_alloc = 0;
	while(NULL != (line = next_line(&lr, &len))){
		if(line[0] == '#'){ continue; }
		line[len] = '\0';
		if(!got_first_line){
			char *saveptr;
			char *field;
			field = STRTOK(line, " \t", &saveptr);
			while(NULL != field){
				(*ncols)++;
				field = STRTOK(NULL, " \t", &saveptr);
			}
			
			nlines_alloc = 256;
			*data = (double*)realloc(*data, sizeof(double) * (*ncols) * nlines_alloc);
			
			got_first_line = 1;
		}
		
		if(*nlines >= nlines_alloc){
			nlines_alloc *= 2;
			*data = (double*)realloc(*data, sizeof(double) * (*ncols) * nlines_alloc);
		}
		
		int nf = 0;
		char *saveptr;
		char *field;
		field = STRTOK(line, " \t", &saveptr);
		while(NULL != field && nf < *ncols){
			(*data)[nf+(*nlines)*(*ncols)] = atof(field);
			++nf;
			field = STRTOK(NULL, " \t", &saveptr);
		}
		(*nlines)++;
	}
	if(!feof(fp)){
		return -1;
	}
	lr_free(&lr);
	return 0;
}

int read_data_1d_cols(FILE *fp, int ncols, int *cols, int *nlines, double **data){
	struct line_reader lr;
	size_t len = 0;
	char *line = NULL;
	int j;
	
	int nlinefields = 64;
	double *linefields = (double*)malloc(sizeof(double) * nlinefields);
	
	*nlines = 0;
	lr_init(&lr, fp);
	int nlines_alloc = 0;
	while(NULL != (line = next_line(&lr, &len))){
		if(line[0] == '#'){ continue; }
		line[len] = '\0';
		
		if(*nlines >= nlines_alloc){
			nlines_alloc *= 2;
			if(0 == nlines_alloc){ nlines_alloc = 64; }
			*data = (double*)realloc(*data, sizeof(double) * (ncols) * nlines_alloc);
		}
		
		int nf = 0;
		char *saveptr = NULL;
		char *field;
		field = STRTOK(line, " \t", &saveptr);
		while(NULL != field){
			if(nf >= nlinefields){
				nlinefields *= 2;
				linefields = (double*)realloc(linefields, sizeof(double) * nlinefields);
			}
			linefields[nf] = atof(field);
			++nf;
			field = STRTOK(NULL, " \t", &saveptr);
		}
		
		for(j = 0; j < ncols; ++j){
			if(cols[j] < nf){
				(*data)[j+(*nlines)*ncols] = linefields[cols[j]];
			}else{
				(*data)[j+(*nlines)*ncols] = 0;
			}
		}
		
		(*nlines)++;
	}
	if(!feof(fp)){
		return -1;
	}
	lr_free(&lr);
	free(linefields);
	return 0;
}

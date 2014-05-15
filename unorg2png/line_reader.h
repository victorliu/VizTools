#ifndef LINE_READER_H_INCLUDED
#define LINE_READER_H_INCLUDED

#include <stddef.h>
#include <stdio.h>

struct line_reader{
	FILE *f;
	char *buf;
	size_t siz;
};

void lr_init(struct line_reader *lr, FILE *f);
char *next_line(struct line_reader *lr, size_t *len);
void lr_free(struct line_reader *lr);

/*
example:

	struct line_reader lr;
	size_t len;
	char *line;
	
	lr_init(&lr, fp);
	
	while(NULL != (line = next_line(&lr, &len))){
		// do something with line
	}
	if(!feof(fp)){
		return -1;
	}
	lr_free(&lr);

*/

#endif /* LINE_READER_H_INCLUDED */

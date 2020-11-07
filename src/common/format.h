#ifndef FORMAT_H
#define FORMAT_H

typedef struct format format;

struct format {
  int dim;
  int* valcount;
};

extern format* 
format_create(int dim);

extern void 
format_read(char* filename, format* fmt);

extern format* 
format_cread(char* filename);

extern void 
format_assign(format* dst, format* src);

extern format* 
format_copy(format* src);

extern void
format_free(format* fmt);

#endif

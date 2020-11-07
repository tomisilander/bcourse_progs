#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "err.h"
#include "format.h"

format* format_create(int dim) {
  format* fmt;
  MEMALL(fmt,1,format);
  fmt->dim = dim;
  MECALL(fmt->valcount, fmt->dim, int);
  return fmt;
}


void format_read(char* filename, format* fmt){
  int c,n,i;
  FILE* fp;

  OPENFILE_OR_DIE(fp,filename,"r");
  n = 1;
  i = -1;
  while(EOF != (c=fgetc(fp))) {
    i += n;
    fmt->valcount[i] += (c == '\t');
    n = (c == '\n');
    
  }
  CLOSEFILE_OR_DIE(fp,filename);
}

format* 
format_cread(char* filename){
  int c,n,dim;
  FILE* fp;
  format* fmt;

  OPENFILE_OR_DIE(fp,filename,"r");

  dim = 0;
  n = 1;
  while(EOF != (c=fgetc(fp))) {
    dim += n;
    n = (c == '\n');
  }

  CLOSEFILE_OR_DIE(fp,filename);

  fmt = format_create(dim);
  format_read(filename, fmt);

  return fmt;
}

void format_assign(format* dst, format* src){
  if (dst == src) return;
  if(dst->dim < src->dim) {
    free(dst->valcount);
    dst->dim = src->dim;
    MEMALL(dst->valcount, dst->dim, int);
  }
  memcpy(dst->valcount, src->valcount, src->dim * sizeof(int));
}

format* 
format_copy(format* src){
  format* fmt = format_create(src->dim);
  format_assign(fmt, src);
  return fmt;
}

void
format_free(format* fmt){
  free(fmt->valcount);
  free(fmt);
}











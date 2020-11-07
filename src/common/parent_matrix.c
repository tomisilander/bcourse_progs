#include <stdlib.h>
#include <string.h>
#include "err.h"
#include "parent_matrix.h"

parent_matrix*
parent_matrix_create(int m){
  parent_matrix* pmx;
  MEMALL(pmx, 1, parent_matrix)
  pmx->m = m;
  pmx->bits_in_int = sizeof(int)*8;
  pmx->one_dim_size = 1 + (m-1) / pmx->bits_in_int;
 
  MECALL(pmx->mx, m * pmx->one_dim_size, unsigned int)
  return pmx; 
}

parent_matrix*
parent_matrix_create_wrap(int m, unsigned int* mx){
  parent_matrix* pmx;
  MEMALL(pmx, 1, parent_matrix)
  pmx->m = m;
  pmx->bits_in_int = sizeof(int)*8;
  pmx->one_dim_size = 1 + (m-1) / pmx->bits_in_int;
 
  pmx->mx = mx;
  return pmx; 
}

void
parent_matrix_assign(parent_matrix* dst, parent_matrix* src){
  memcpy(dst->mx, src->mx, 
	 src->m * src->one_dim_size * sizeof(int));
}

parent_matrix*
parent_matrix_copy(parent_matrix* src){
  parent_matrix* dst = parent_matrix_create(src->m);
  parent_matrix_assign(dst,src);
  return dst;
}

void
parent_matrix_print(parent_matrix* pmx, FILE* fp){
  int c = 0;
  for(c=0; c<pmx->m; ++c) {
    int p;
    for(p=0; p<pmx->m; ++p){
    	fprintf(fp, "%d%c",
    		IS_PARENT(c,p,pmx), (p==pmx->m-1)?'\n':' ');
    }
  }
}

void
parent_matrix_free(parent_matrix* pmx){
  free(pmx->mx);
  free(pmx);
}

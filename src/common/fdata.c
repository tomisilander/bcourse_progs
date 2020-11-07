#include <stdio.h>
#include <stdlib.h>
#include "err.h"
#include "fdata.h"

fdata* fdata_create(int datacount, int attcount){
  fdata* dt;
  int i,j;

  MEMALL(dt, 1, fdata);
  dt->N = datacount;
  dt->om = attcount;
  dt->m = attcount;
  MEMALL(dt->d, dt->N*dt->m, float);
  MEMALL(dt->rowmap, dt->N, int);
  MEMALL(dt->colmap, dt->m, int);

  for(j=0; j<dt->N; ++j) dt->rowmap[j] = j;
  for(i=0; i<dt->m; ++i) dt->colmap[i] = i;

  return dt;
}

fdata* fdata_new(fdata* old_dt){
  int i,j;
  fdata* dt;
  MEMALL(dt, 1, fdata)
  dt->N = old_dt->N;
  dt->om = old_dt->om;
  dt->m = old_dt->m;
  dt->d = old_dt->d;
  MEMALL(dt->rowmap, dt->N, int)
  for(j=0; j<dt->N; ++j){
    dt->rowmap[j] = old_dt->rowmap[j];
  }
  MEMALL(dt->colmap, dt->m, int)
  for(i=0; i<dt->m; ++i){
    dt->colmap[i] = old_dt->colmap[i];
  }

  return dt;
}

fdata* fdata_new_cols(fdata* old_dt, int* sel){
  int i,j,c;
  fdata* dt;
  MEMALL(dt, 1, fdata)
  dt->N = old_dt->N;  
  dt->om = old_dt->om;
  dt->m = old_dt->m;
  dt->d = old_dt->d;
  MEMALL(dt->rowmap, dt->N, int)
  for(j=0; j<dt->N; ++j){
    dt->rowmap[j] = old_dt->rowmap[j];
  }
  MEMALL(dt->colmap, dt->m, int)
  c=0;
  for(i=0; i<dt->m; ++i){
    if(sel[i]) dt->colmap[c++] = old_dt->colmap[i];
  }
  dt->m = c;

  return dt;
}

fdata* fdata_new_rows_by_vals(fdata* old_dt, int* sel){
  int i,j,r;
  fdata* dt;
  MEMALL(dt, 1, fdata)
  dt->N = old_dt->N;  
  dt->om = old_dt->om;
  dt->m = old_dt->m;
  dt->d = old_dt->d;

  MEMALL(dt->colmap, dt->m, int);

  for(i=0; i<dt->m; ++i){
    dt->colmap[i] = old_dt->colmap[i];
  }

  MEMALL(dt->rowmap, dt->N, int);

  r = 0;
  for(j=0; j<dt->N; ++j){
    int match = 1;
    for(i=0; i<dt->m; ++i){
      match &= ((sel[i] == -1) || (FD(old_dt,j,i) == sel[i]));
    }
    if(match) dt->rowmap[r++] = old_dt->rowmap[j];
  }
  
  dt->N = r;

  return dt;
}

void fdata_read(char* filename, fdata* dt){
  FILE* fp;  
  int j;

  if(NULL == (fp = fopen(filename, "r"))){
    fprintf(stderr, "can't open %s\n", filename);
    exit(CANNOT_OPEN_DATAFILE);
  }

  for(j=0; j<dt->N; ++j){
    int i;
    dt->rowmap[j] = j;
    for(i=0; i<dt->m; ++i){
      float tmp;
      if(1 != fscanf(fp, "%f", &tmp)){
	fprintf(stderr, "error while reading %s\n", filename);
	exit(ERROR_IN_DATAREAD);	
      }
      dt->colmap[i] = i;
      FD(dt,j,i) = tmp;
    }
  }

  if(fclose(fp)){
    fprintf(stderr, "can't close %s\n", filename);
    exit(CANNOT_CLOSE_DATAFILE);
  }
}

void fdata_old(fdata* dt){
  free (dt->rowmap);
  free (dt->colmap);
  free (dt);
}

void fdata_free(fdata* dt){
  free (dt->d);
  free (dt->rowmap);
  free (dt->colmap);
  free (dt);
}













#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "err.h"
#include "data.h"

data* data_create(int datacount, int attcount){
  data* dt;

  MEMALL(dt, 1, data);
  dt->N = datacount;
  dt->m = attcount;
  MEMALL(dt->d, dt->N*dt->m, char);

  return dt;
}

data* data_cread(int datacount, int attcount, char* filename){
  data* dt = data_create(datacount, attcount);
  data_read(filename, dt);
  return dt;
}

data* data_copy(data* old_dt){
  data* dt;
  MEMALL(dt, 1, data);
  MEMALL(dt->d, old_dt->N * old_dt->m, char);
  memcpy(dt->d, old_dt->d, dt->N*dt->m*sizeof(char));
  return dt;
}

void data_read(char* filename, data* dt){
  FILE* fp;  
  int j;

  OPENFILE_OR_DIE(fp,filename,"r");

  for(j=0; j<dt->N; ++j){
    int i;
    for(i=0; i<dt->m; ++i){
      int tmp;
      if(1 != fscanf(fp, "%d", &tmp)){
	fprintf(stderr, "error while reading %s\n", filename);
	exit(ERROR_IN_DATAREAD);	
      }
      D(dt,j,i) = (char) tmp;
    }
  }

  CLOSEFILE_OR_DIE(fp,filename);
}

void data_write(char* filename, data* dt){
  FILE* fp;  
  int j;

  OPENFILE_OR_DIE(fp,filename,"w");

  for(j=0; j<dt->N; ++j){
    int i;
    for(i=0; i<dt->m; ++i){
      fprintf(fp,"%d%c", D(dt,j,i), (i==dt->m-1) ? '\n' : '\t');
    }
  }

  CLOSEFILE_OR_DIE(fp,filename);
}

void data_free(data* dt){
  free (dt->d);
  free (dt);
}

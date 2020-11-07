#include <stdio.h>
#include <stdlib.h>
#include "err.h"
#include "format.h"
#include "data.h"
#include "bane.h"


int main(int argc, char* argv[]){
  int n;
  data* dt;
  format* fmt;
  double ess;
  bane* bnf;

  FILE* fp;

  if(argc !=6){
    fprintf(stderr, 
	    "Usage: %s vdfile datafile datacount ess structfile\n", argv[0]);
    exit(-1);
  }

  fmt = format_cread(argv[1]);

  n = atoi(argv[3]);
  dt = data_create(n,fmt->dim);
  data_read(argv[2],dt);
  ess = atof(argv[4]);

  bnf = bane_create_forest(fmt,ess,dt);
  
  OPENFILE_OR_DIE(fp,argv[5],"w");
  bane_write_structure(bnf,fp);
  CLOSEFILE_OR_DIE(fp,argv[5]);

  bane_gather_full_ss(bnf, dt);
  printf("%f\n",bane_get_score(bnf, ess, NULL));

  bane_free(bnf);
  format_free(fmt);
  data_free(dt);
  
  return 0;
}

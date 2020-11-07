#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "err.h"
#include "format.h"
#include "data.h"
#include "bane.h"


int main(int argc, char* argv[]){

  format* fmt;
  data *dt;
  bane* bn;
  double* scoreboard;
  double ess;
  int i;
  FILE* fp;

  if(argc != 6){
    fprintf(stderr, 
	    "Usage: %s vdfile datafile datacount ess struct\n", 
	    argv[0]);
    exit(-1);
  }

  fmt = format_cread(argv[1]);
  dt = data_cread(atoi(argv[3]), fmt->dim, argv[2]);
  bn = bane_create_from_format(fmt);
  OPENFILE_OR_DIE(fp,argv[5],"r"); 
  bane_read_structure(bn,fp);
  CLOSEFILE_OR_DIE(fp,argv[5]);
  ess = atof(argv[4]);

  MECALL(scoreboard, fmt->dim, double);
  BDe_score(bn, dt, ess, scoreboard);

  for(i=0; i<fmt->dim; ++i){
    printf("%g\n", scoreboard[i]);
  }

  format_free(fmt);
  data_free(dt);
  bane_free(bn);
  free(scoreboard);

  return 0;
}











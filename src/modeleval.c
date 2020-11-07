#include <stdio.h>
#include <stdlib.h>
#include "err.h"
#include "format.h"
#include "data.h"
#include "bane.h"
#include "node.h"

extern double lgamma(double);

int main(int argc, char* argv[]){

  format* fmt;
  data* dt;
  bane* bn;
  double ess;
  FILE* fp;

  if(argc != 6){
    fprintf(stderr, 
	    "Usage: %s vdfile datafile datacount structfile ess\n", 
	    argv[0]);
    exit(-1);
  }

  fmt = format_cread(argv[1]);
  dt  = data_cread(atoi(argv[3]), fmt->dim, argv[2]);
  bn  = bane_create_from_format(fmt);
  OPENFILE_OR_DIE(fp,argv[4],"r"); 
  bane_read_structure(bn,fp);
  CLOSEFILE_OR_DIE(fp,argv[4]);
  ess = atof(argv[5]);

  printf("%.6f\n",BDe_score(bn, dt, ess, NULL));

  bane_free(bn);
  data_free(dt);
  format_free(fmt);

  return 0;
}











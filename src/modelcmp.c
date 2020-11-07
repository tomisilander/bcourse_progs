#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "err.h"
#include "format.h"
#include "data.h"
#include "bane.h"
#include "node.h"

int main(int argc, char* argv[]){

  format* fmt;
  data *dt;
  bane* ref_bn;
  double ref_score;
  double ess;
  int b;
  FILE* fp;

  if(argc <= 6){
    fprintf(stderr, 
	    "Usage: %s vdfile datafile datacount "
	    "ref struct ess cmpstruct ..\n", 
	    argv[0]);
    exit(-1);
  }

  fmt = format_cread(argv[1]);
  dt = data_create(atoi(argv[3]), fmt->dim);
  data_read(argv[2], dt);
  ref_bn = bane_create_from_format(fmt);
  OPENFILE_OR_DIE(fp,argv[4],"r"); 
  bane_read_structure(ref_bn,fp);
  CLOSEFILE_OR_DIE(fp,argv[4]);
  ess = atof(argv[5]);

  bane_gather_full_ss_in_order(ref_bn, dt);
  ref_score = bane_get_score(ref_bn,ess,NULL);

  for(b=0; b<argc-6; ++b){
    bane* bn = bane_create_from_format(fmt);
    OPENFILE_OR_DIE(fp,argv[6+b],"r"); 
    bane_read_structure(bn,fp);
    CLOSEFILE_OR_DIE(fp,argv[6+b]);

    bane_gather_full_ss_in_order(bn, dt);
    printf("%g\n", exp(bane_get_score(bn,ess,NULL) - ref_score));
    bane_free(bn);
  }

  format_free(fmt);
  data_free(dt);
  bane_free(ref_bn);

  return 0;
}











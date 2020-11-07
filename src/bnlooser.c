#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include "err.h"
#include "format.h"
#include "data.h"
#include "bane.h"

int main(int argc, char* argv[]){

  format* fmt;
  data* dt;
  bane* bn;
  double ess;


  FILE* fp;
  int j;
  double log_score;

  char** argp;
  int opt;
  char v_level = 'q';
  double* logterms = NULL;

  while (-1 != (opt = getopt(argc, argv, "vV"))) v_level = opt;
  argp = argv + optind;

  if(argc != 6 && argc != 7){
    fprintf(stderr, 
	    "Usage: %s [-v | -V] formatfile datafile datacount strfile ess\n", 
	    argv[0]);
    exit(-1);
  }

  fmt = format_cread(argp[0]);
  dt = data_create(atoi(argp[2]), fmt->dim);
  data_read(argp[1], dt);

  bn = bane_create_from_format(fmt);
  OPENFILE_OR_DIE(fp,argp[3],"r");
  bane_read_structure(bn,fp);
  CLOSEFILE_OR_DIE(fp,argp[3]);                                                 
  ess = atof(argp[4]);

  log_score = 0;
  bane_gather_full_ss(bn,dt);

  if(v_level == 'V') MEMALL(logterms, fmt->dim, double);

  for(j=0; j<dt->N;++j){
    double lp = bane_logprob(bn,dt,j,ess,logterms);
    bane_del_ss(bn,dt,j);
    log_score += lp;
    bane_add_ss(bn,dt,j);

    if(v_level != 'q'){
      fprintf(stdout, "%f", lp);
      if(v_level == 'V'){
	int i;for(i=0; i<fmt->dim; ++i) fprintf(stdout," %f", logterms[i]);
      }
      fprintf(stdout,"\n");
    }
  }

  if(v_level == 'V') free(logterms);

  fprintf(stdout, "%f\n", log_score);

  data_free(dt);
  bane_free(bn);
  format_free(fmt);

  return 0;
}

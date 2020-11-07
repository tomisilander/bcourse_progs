#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include "err.h"
#include "format.h"
#include "data.h"
#include "bane.h"

int main(int argc, char* argv[]){

  format* fmt;
  data* trndt;
  data* tstdt;
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

  if(argc != 8 && argc != 9){
    fprintf(stderr, 
	    "Usage: %s [-v | -V] formatfile strfile trnidt trndc ess tstidt tstdc\n", 
	    argv[0]);
    exit(-1);
  }

  fmt = format_cread(argp[0]);

  bn = bane_create_from_format(fmt);
  OPENFILE_OR_DIE(fp,argp[1],"r");
  bane_read_structure(bn,fp);
  CLOSEFILE_OR_DIE(fp,argp[1]);                                                 

  trndt = data_create(atoi(argp[3]), fmt->dim);
  data_read(argp[2], trndt);

  ess = atof(argp[4]);

  tstdt = data_create(atoi(argp[6]), fmt->dim);
  data_read(argp[5], tstdt);

  log_score = 0;
  bane_gather_full_ss(bn,trndt);

  if(v_level == 'V') MEMALL(logterms, fmt->dim, double);

  for(j=0; j<tstdt->N;++j){
    double lp = bane_logprob(bn,tstdt,j,ess,logterms);
    log_score += lp;

    if(v_level != 'q'){
      fprintf(stdout, "%f", lp);
      if(v_level == 'V'){
	int i;for(i=0; i<fmt->dim; ++i) fprintf(stdout," %f", logterms[i]);
      }
      fprintf(stdout,"\n");
    }
  }

  if(v_level == 'V') free(logterms);

  fprintf(stdout, "%f\n", log_score/tstdt->N);

  data_free(tstdt);
  data_free(trndt);
  bane_free(bn);
  format_free(fmt);

  return 0;
}

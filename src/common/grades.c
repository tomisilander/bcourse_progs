#include <stdio.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "err.h"
#include "grades.h"

grades* grades_create(int N, int K)
{
  grades* gs;
  MEMALL(gs,1,grades);
  gs->N = N;
  gs->K = K;
  MECALL(gs->g, N*K, double);
  return gs;
}

void grades_free(grades* gs)
{
  free(gs->g);
  free(gs);
}

void grades_read(char* filename, grades* gs)
{
  int j,k;
  FILE* fp;

  OPENFILE_OR_DIE(fp,filename,"r");
  for(j=0; j<gs->N; ++j) {
    for(k=0; k<gs->K; ++k){
      fscanf(fp, "%lf", &(GJK(gs,j,k)));
    }
  }
  CLOSEFILE_OR_DIE(fp,filename);
}

void grades_write(char* filename, grades* gs)
{
  int j,k;
  FILE* fp;

  OPENFILE_OR_DIE(fp,filename,"w");
  for(j=0; j<gs->N; ++j) {
    for(k=0; k<gs->K; ++k){
      fprintf(fp, "%lf ", GJK(gs,j,k));
    }
    fprintf(fp, "\n");
  }
  CLOSEFILE_OR_DIE(fp,filename);
}

grades* grades_cread(int N, int K, char* filename)
{
  grades* gs = grades_create(N,K);
  grades_read(filename,gs);
  return gs;
}

void grades_init(gsl_rng* rng, grades* gs)
{
  int j,k;
  double* t_ss;

  MEMALL(t_ss, gs->K, double);
  for(k=0; k<gs->K; ++k) t_ss[k] = 1.0;

  for(j=0; j<gs->N; ++j) {
    gsl_ran_dirichlet(rng, gs->K, t_ss, GJ(gs,j));
  }
  free(t_ss);
}

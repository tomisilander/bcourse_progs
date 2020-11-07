#include <stdio.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "err.h"
#include "format.h"
#include "profiles.h"


profiles* profiles_create(int K, format* fmt)
{
  profiles* profs;
  int i;
  MEMALL(profs,1,profiles);
  profs->K = K;
  MEMALL(profs->var_ofs, fmt->dim, int);
  MEMALL(profs->valcs,   fmt->dim, int);
  profs->var_prof_len = 0;
  for(i=0; i<fmt->dim; ++i){
    profs->var_ofs[i] = profs->var_prof_len;
    profs->var_prof_len += K*(profs->valcs[i] = fmt->valcount[i]);
  }
  MECALL(profs->p, fmt->dim * profs->var_prof_len, double);
  return profs;
}

void profiles_free(profiles* profs)
{
  free(profs->p);
  free(profs->var_ofs);
  free(profs->valcs);
  free(profs);
}

void profiles_read(char* filename, format* fmt, profiles* profs)
{
  FILE* fp;

  OPENFILE_OR_DIE(fp,filename,"r");
  LOOP_KIL(profs->K, fmt){
    fscanf(fp, "%lf", &(PROFSKIL(profs,k,i,l)));
  } END_KIL;
  CLOSEFILE_OR_DIE(fp,filename);
}

void profiles_write(char* filename, format* fmt, profiles* profs)
{
  FILE* fp;

  OPENFILE_OR_DIE(fp,filename,"w");
  LOOP_KIL(profs->K, fmt){
    fprintf(fp, "%lf%c", 
	    PROFSKIL(profs,k,i,l), 
	    (i==fmt->dim-1 && l==fmt->valcount[i]-1) ? '\n':' '); 
    
  } END_KIL;
  CLOSEFILE_OR_DIE(fp,filename);
}

profiles* profiles_cread(char* filename, int K, format* fmt)
{
  profiles* profs = profiles_create(K,fmt);
  profiles_read(filename, fmt, profs);
  return profs;
}

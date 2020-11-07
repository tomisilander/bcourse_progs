#include <stdio.h>
#include <errno.h>
#include <math.h>
#include <sys/types.h>
#include <unistd.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "err.h"
#include "format.h"
#include "data.h"
#include "node.h"
#include "sigcom.h"
#include "grades.h"
#include "profiles.h"


void write_report_n_model(unsigned long t, format* fmt, 
			  grades* gs, profiles* profs,
			  char* report_fn, char* grades_fn, char* profiles_fn){
  FILE* fp;

  grades_write(grades_fn, gs);
  profiles_write(profiles_fn, fmt, profs);

  OPENFILE_OR_DIE(fp,report_fn,"w");
  fprintf(fp,"%lu\n",t);
  CLOSEFILE_OR_DIE(fp,report_fn);
}

#define MULSEL(rng,sel,len,probs,nzer)\
{\
  double x = gsl_rng_uniform(rng) * nzer;\
  double psm = 0.0;\
  int k;\
  for(k=0; k<len; ++k){\
    psm += probs[k];\
    if(psm >= x){\
      sel = k;\
      break;\
    }\
  }\
}

#define SWAP(TYPE,X,Y) {TYPE tmp=(X); (X)=(Y); (Y)=tmp;}

#define INCAVG(C,AVG,NEW) \
{\
  AVG *= (C) - 1;\
  AVG += (NEW);\
  AVG /= (C);\
}

void gibbsgen(gsl_rng* rng,
	      format* fmt, data* d, int K,
	      double a, double b,
	      grades* gs,     profiles* profs,     /* current params */
	      grades* gs_new, profiles* profs_new, /* new     params */
	      double* r)
{
  int k,j,i;
  
  LOOP_KIL(K, fmt) { /* Set the equivalent sample size to b.*/
    PROFSKIL(profs_new,k,i,l) = b;
  } END_KIL;
    
  for(j=0; j<d->N; ++j){
    double* g = GJ(gs, j);
    double* g_new = GJ(gs_new, j);
    for(k=0; k<K; ++k) g_new[k] = a; /* Set the equivalent sample size to a.*/
    for(i=0; i<fmt->dim; ++i){
      int dji = D(d,j,i);
      int z; 
      double nzer = 0.0;
      if (dji == MISSINGVALUE) continue;
      for(k=0; k<K; ++k) nzer += r[k] = g[k] * PROFSKIL(profs,k,i,dji);
      MULSEL(rng,z,K,r,nzer);             /* sample z by r        */
      g_new[z] += 1.0;                    /* and update suffstats */
      PROFSKIL(profs_new,z,i,dji) += 1.0; /* (for profs too)      */
    }
    gsl_ran_dirichlet(rng, K, g_new, g_new); /* sample new grades */
  }

  for(k=0; k<K; ++k){
    for(i=0; i<fmt->dim; ++i){
      double* profki_new = PROFSKI(profs_new,k,i);
      gsl_ran_dirichlet(rng, fmt->valcount[i], profki_new, profki_new);
    }
  }
}


void gibbsample(gsl_rng* rng,
		format* fmt, data* dt, 
		int K, double a, double b, 
		int burnin, int thinnin,
		char* report_fn, char* grades_fn, char* profiles_fn)
{

  /* current, new and resulting grades of membership */
  grades* gs       = grades_create(dt->N, K);
  grades* gs_new   = grades_create(dt->N, K);
  grades* gs_stat  = grades_create(dt->N, K);

  /* current, new and resulting grades of profiles */

  profiles* profs       = profiles_create(K,fmt);
  profiles* profs_new   = profiles_create(K,fmt);
  profiles* profs_stat  = profiles_create(K,fmt);

  /* NB. We collect sufficient statistics to gs_new and profs_new too */

  double* r;           /* work space for GOM distribution per case */
  unsigned long round, c;

  MEMALL(r, K, double);

  /* REAL LOOP */
  
  round = c = 0;
  while(!end_flag_set){
    ++round;
    { /* Call Gibbs and swap old and new */
      gibbsgen(rng,
	       fmt, dt, K, 
	       a, b,
	       gs, profs,
	       gs_new, profs_new,
	       r);
      SWAP(grades*, gs, gs_new);
      SWAP(profiles*, profs, profs_new);
    }

    if((round < burnin) || (round % thinnin != 0)) continue;

    /* UPDATE STATS */
    ++c;
    { /* UPDATE GRADES */
      int j,k;
      for(j=0; j<dt->N; ++j){
	for(k=0; k<K; ++k){
	  INCAVG(c, GJK(gs_stat,j,k), GJK(gs,j,k));
	}
      }
    }
    
    LOOP_KIL(K, fmt) { /* UPDATE PROFILES */
      INCAVG(c, PROFSKIL(profs_stat,k,i,l), PROFSKIL(profs,k,i,l));
    } END_KIL;
  
    
    /* REPORT AND END OF RUN HANDLING */
    if(rep_flag_set || alrm_flag_set) { 
      rep_flag_set = 0;
      write_report_n_model(c, fmt, gs_stat, profs_stat,
			   report_fn, grades_fn, profiles_fn);
      if(alrm_flag_set){
	end_flag_set = 1;
      }
    }
}

  /* FREE STUFF */

  grades_free(gs);
  grades_free(gs_new);
  grades_free(gs_stat);

  profiles_free(profs);
  profiles_free(profs_new);
  profiles_free(profs_stat);

  free(r);
}

int main(int argc, char* argv[]){
  int n, K;
  data* dt;
  format* fmt;
  double a, b;
  FILE* fp;
  int burnin, thinnin;
  gsl_rng * rng;

  if((argc !=13) && (argc!=14)){
    fprintf(stderr, 
	    "Usage: %s vdfile datafile datacount a b K "
	    "reportfile grades profiles "
	    "burnin thinnin searchtime(<0 forever) [pidfile]\n", 
	    argv[0]);
    exit(-1);
  }
   
  set_signal_handlers();

  if(argc==13) {
    fprintf(stdout,"%d\n",getpid());
    fclose(stdout);
  } else {
    fclose(stdin);
    fclose(stdout);
    fclose(stderr);
    OPENFILE_OR_DIE(fp,argv[13],"w");
    fprintf(fp,"%d\n",getpid());
    CLOSEFILE_OR_DIE(fp,argv[10]);
  }
  
  fmt = format_cread(argv[1]);

  n = atoi(argv[3]);
  K = atoi(argv[6]);
  dt = data_cread(n,fmt->dim,argv[2]);

  a = atof(argv[4]);
  b = atof(argv[5]);
  burnin  = atoi(argv[10]);
  thinnin = atoi(argv[11]);

  setalarm(atoi(argv[12]));

  { /* GSL RANDOM NUMBER GENERATOR SETUP */
    const gsl_rng_type * rngT;
    gsl_rng_env_setup();
    rngT = gsl_rng_default;
    rng = gsl_rng_alloc (rngT);
  }

  gibbsample(rng,
	     fmt, dt, 
	     K, a, b, 
	     burnin, thinnin,
	     argv[7], argv[8], argv[9]);

  format_free(fmt);
  data_free(dt);
  
  return 0;
}

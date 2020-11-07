#include <stdio.h>
#include <math.h>
#include <signal.h>
#include <errno.h>
#include <unistd.h>
#include "err.h"
#include "format.h"
#include "data.h"
#include "score_hashtable.h"
#include "nbmodel.h"
#include "bitset.h"
#include "nb_selector.h"
#include "sigcom.h"

typedef struct search_stats search_stats;

struct search_stats {
  unsigned long t;
  unsigned long tn;
  double best_score_zos;
  double best_score_lgs;
  unsigned long tn_when_best;
  double hit_rate_rwa;
  int* sel;
  int n;
};

search_stats* search_stats_create(int* sel, int n) {
  search_stats* me;
  MECALL(me,1,search_stats);
  BITSET_COPY(me->sel, sel, n);
  me->n = n;
  return me;
}

void search_stats_free(search_stats* me) {
  if(me->sel != NULL) {
    BITSET_FREE(me->sel);
  }
  free(me);
}

void search_stats_assign(search_stats* dst, search_stats* src) {
  dst->t = src->t;
  dst->tn = src->tn;
  dst->best_score_zos = src->best_score_zos;
  dst->best_score_lgs = src->best_score_lgs;
  dst->tn_when_best = src->tn_when_best;
  dst->hit_rate_rwa = src->hit_rate_rwa;
  BITSET_ASSIGN(dst->sel, src->sel, src->n);
}

#define BETTERSCORE(Z1,L1,Z2,L2) \
( ((Z1) > (Z2))  || ( ((Z1) == (Z2)) && ((L1) > (L2)) ) )



void write_report_n_structure(search_stats* st, search_stats* stp, int classix,
			      char* reportfilename, char* structfilename){
  FILE* fp;

  /* Structure first since it is bigger */

  nb_selector_write(st->sel, st->n, classix, structfilename);

  /* Report second */

  OPENFILE_OR_DIE(fp,reportfilename,"w");

  fprintf(fp,"%lu %lu %lu %f %f %f %f\n", 
	  st->tn, 
	  st->tn - stp->tn, 
	  st->tn - st->tn_when_best, 
	  stp->best_score_zos,
	  st->best_score_zos,
	  stp->best_score_lgs,
	  st->best_score_lgs
	  );

  CLOSEFILE_OR_DIE(fp,reportfilename);

  search_stats_assign(stp,st);

}

static
void t_times_more_brute(format* fmt, nbmodel* m, double ess, data* dt, 
			int t_limit, 
			search_stats* stl, search_stats* stg){

  int t; 
  int* sel = stl->sel; 

  for(t = 0; 
      (t<t_limit) && !(rep_flag_set || end_flag_set || alrm_flag_set); 
      ++t){
    double new_score_zos;
    double new_score_lgs;

    int i;
    int xor = stl->t ^ (stl->t-1);
    for(i=0; xor >>= 1; ++i);
    if(i >= m->root) ++i;

    if(i >= fmt->dim) {
      sleep(1);
      continue;
    }

    /* Change */
    if(BITSET_ON(sel,i)) BITSET_UNSET(sel,i); else BITSET_SET(sel,i);

    nbmodel_scores(fmt, m, sel, ess, dt, &new_score_zos, &new_score_lgs);
    ++ stl->tn; ++ stg->tn; ++ stl->t; ++ stg->t;

    if(BETTERSCORE(new_score_zos, new_score_lgs, 
		   stl->best_score_zos, stl->best_score_lgs)) { 
      /* fprintf(stderr,"%d %f\n",ch,new_score); */
      stl->best_score_zos = new_score_zos;
      stl->best_score_lgs = new_score_lgs;
      stl->tn_when_best = stl->tn;
      if(BETTERSCORE(new_score_zos, new_score_lgs, 
		     stg->best_score_zos, stg->best_score_lgs)) { 
	stg->best_score_zos = new_score_zos;
	stg->best_score_lgs = new_score_lgs;
	stg->tn_when_best = stg->tn;
	BITSET_ASSIGN(stg->sel, sel, stg->n);
      }
    }
  }
}



void brute_search(format* fmt, nbmodel* m, double ess, data *dt,
		  char* rfn, char* sfn){

  search_stats* stg; /* global search stats */
  search_stats* stl; /* search stats for current round - local*/
  search_stats* stp; /* search stats at previous report */


  int* sel;
  BITSET_CREATE(sel, fmt->dim);
  BITSET_SET(sel,m->root); 
  stg = search_stats_create(sel,fmt->dim);
  stl = search_stats_create(sel,fmt->dim);
  stp = search_stats_create(sel,fmt->dim);
  BITSET_FREE(sel);

  nbmodel_scores(fmt, m, stl->sel, ess, dt,
		 &stl->best_score_zos, &stl->best_score_lgs);
  
  stg->best_score_zos = stl->best_score_zos;
  stg->best_score_lgs = stl->best_score_lgs;

  ++stl->t;
  ++stg->t;
  ++stl->tn;
  ++stg->tn;

  stp->t = -1;
  stp->tn = -1;
  stp->best_score_zos = stg->best_score_zos;
  stp->best_score_lgs = stg->best_score_lgs;

  while(!end_flag_set){
    int tlimit = 100;

    t_times_more_brute(fmt, m, ess, dt, tlimit, stl, stg);

    if(rep_flag_set || alrm_flag_set) {
      rep_flag_set = 0;
      write_report_n_structure(stg, stp, m->root, rfn, sfn);

      if(alrm_flag_set){
	end_flag_set = 1;
      }
    }
  }

  search_stats_free(stl);
  search_stats_free(stg);
  search_stats_free(stp);

}

int main(int argc, char* argv[]){
  int n;
  data* dt;
  format* fmt;
  double ess;
  FILE* fp;
  int classix;
  nbmodel* m;

  if((argc !=9) && (argc!=10)){
    fprintf(stderr, 
	    "Usage: %s vdfile datafile datacount classix ess "
	    "reportfile structfile searchtime(<0 forever) [pidfile]\n", 
	    argv[0]);
    exit(-1);
  }
   
  set_signal_handlers();

  if(argc==9) {
    fprintf(stdout,"%d\n",getpid());
  } else {
    fclose(stdin);
    fclose(stdout);
    fclose(stderr);
    OPENFILE_OR_DIE(fp,argv[9],"w");
    fprintf(fp,"%d\n",getpid());
    CLOSEFILE_OR_DIE(fp,argv[9]);
  }
  
  fmt = format_cread(argv[1]);

  n = atoi(argv[3]);
  dt = data_create(n,fmt->dim);
  data_read(argv[2],dt);

  classix = atoi(argv[4]);
  ess = atof(argv[5]);

  setalarm(atoi(argv[8]));

  m = nbmodel_create(fmt,dt,classix); 
  brute_search(fmt, m, ess, dt, argv[6], argv[7]);

  nbmodel_free(fmt, &m);
  format_free(fmt);
  data_free(dt);

  return 0;
}

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

int better_to_restart(data* dt, search_stats* st){
  int res = 0;
  long vain = st->tn - st->tn_when_best;

  if( (st->tn ==  0) 
      || ((1.0 * vain / st->tn > 0.5)  && (st->t > (unsigned) 10*dt->m))){
    /* fprintf(stderr,"rs1:%lu %ld\n", st->t, vain); */
    res = 1;
  }

  if(((long) st->t > 100 * dt->m) && (st->hit_rate_rwa > 0.95)){
    /* fprintf(stderr,"rs2: %f\n", st->hit_rate_rwa); */
    res = 1;
  }
  return res;
}

static
void remove_bad_arcs(format* fmt, nbmodel* m, double ess, 
		     data* dt,  score_hashtable* sht, search_stats* stg){

  int all_removed = 0; 
  int* sel = stg->sel;

  while(!all_removed){
    int i;
    all_removed = 1;
    for(i=0; i<fmt->dim; ++i){
      if(i == m->root) continue;
      if (BITSET_ON(sel,i)){
	double old_score;
	BITSET_UNSET(sel,i);
	if(!score_hashtable_get(sht, sel, &old_score)) {
	  double zos = 0;
	  double lgs = 0;
	  nbmodel_scores(fmt, m, sel, ess, dt,&zos,&lgs);
	  if(BETTERSCORE(zos,lgs,stg->best_score_zos,stg->best_score_lgs)){
	    stg->tn_when_best = stg->tn;
	    stg->best_score_zos = zos;
	    stg->best_score_lgs = lgs;
	    all_removed = 0;
	  } else {
	    BITSET_SET(sel,i);
	  }
	  ++ stg->tn;
	} else {
	    BITSET_SET(sel,i);
	}
      }
    }
  }
}

static
void t_times_greedy_step(format* fmt, nbmodel* m, double ess, data* dt, 
			 int t_limit,
			 score_hashtable* sht, score_hashtable* sht2,
			 search_stats* stl, search_stats* stg){

  int t; 
  int* sel = stl->sel; 

  for(t = 0; 
      (t<t_limit) && !(rep_flag_set || end_flag_set || alrm_flag_set); 
      ++t){
    double new_score_zos;
    double new_score_lgs;
    int hthit;

    int i;

    i = rand()%(fmt->dim-1);
    if(i >= m->root) ++i;


    /* Change */
    if(BITSET_ON(sel,i)) BITSET_UNSET(sel,i); else BITSET_SET(sel,i);

    hthit = 0;
    if(score_hashtable_get(sht, sel, &new_score_zos)) {
      score_hashtable_get(sht2, sel, &new_score_lgs);
      hthit = 1;
    } else {
      nbmodel_scores(fmt, m, sel, ess, dt, &new_score_zos, &new_score_lgs);
      score_hashtable_put(sht,  sel, new_score_zos);
      score_hashtable_put(sht2, sel, new_score_lgs);
      ++ stl->tn;
      ++ stg->tn;
    }

    stl->hit_rate_rwa += 0.01 * (hthit - stl->hit_rate_rwa);

    ++ stl->t;
    ++ stg->t;

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
    } else { /* change it back */
      if(BITSET_ON(sel,i)) BITSET_UNSET(sel,i); else BITSET_SET(sel,i);
    }
  }
}



void greedy_search(format* fmt, nbmodel* m, double ess, data *dt,
		   char* rfn, char* sfn){

  search_stats* stg; /* global search stats */
  search_stats* stl; /* search stats for current round - local*/
  search_stats* stp; /* search stats at previous report */

  score_hashtable* sht;
  score_hashtable* sht2;

  int* sel;
  BITSET_CREATE(sel, fmt->dim);
  BITSET_SET(sel,m->root); 
  stg = search_stats_create(sel,fmt->dim);
  stl = search_stats_create(sel,fmt->dim);
  stp = search_stats_create(sel,fmt->dim);

  BITSET_FREE(sel);

  { /* Create hashtables */
    int keysize;
    int items_per_pos;
    keysize = BITSET_INTCOUNT(fmt->dim);
    items_per_pos = score_hashtable_items_for_mem(2000000, keysize);
    if(items_per_pos < 1) items_per_pos = 1;
    if(items_per_pos > 10) items_per_pos = 10;
    sht = score_hashtable_create(keysize, items_per_pos);
    sht2 = score_hashtable_create(keysize, items_per_pos);
  }

  nbmodel_scores(fmt, m, stl->sel, ess, dt,
		 &stl->best_score_zos, &stl->best_score_lgs);

  stg->best_score_zos = stl->best_score_zos;
  stg->best_score_lgs = stl->best_score_lgs;
  score_hashtable_put(sht, stl->sel, stg->best_score_zos); 
  score_hashtable_put(sht2,stl->sel, stg->best_score_lgs); 
  stl->hit_rate_rwa = 0;

  ++stl->t;
  ++stg->t;
  ++stl->tn;
  ++stg->tn;

  stp->t = -1;
  stp->tn = -1;
  stp->best_score_zos = stg->best_score_zos;
  stp->best_score_lgs = stg->best_score_lgs;

  while(!end_flag_set){

    t_times_greedy_step(fmt, m, ess, dt, 100, sht, sht2, stl, stg);

    if(rep_flag_set || alrm_flag_set) {
      rep_flag_set = 0;
      if(alrm_flag_set) {
	remove_bad_arcs(fmt, m, ess, dt, sht, stg);
      }
      write_report_n_structure(stg, stp, m->root, rfn, sfn);

      if(alrm_flag_set){
	end_flag_set = 1;
      }
    }

    if(better_to_restart(dt,stl)){
      int n;
      int t;

      /* Take some good one out of hashtable and mutilate it */
      
      n = 1 + rand() % (1 + (int) (sht->keycount * 0.01));

      BITSET_ASSIGN(stl->sel,score_hashtable_get_nth_key(sht,n), stl->n);

      for(t=0; t<fmt->dim / 3; ++t) {
	int i = rand()%(fmt->dim-1);
	if(i >= m->root) ++i;
	if(BITSET_ON(stl->sel,i)) 
	  BITSET_UNSET(stl->sel,i); 
	else 
	  BITSET_SET(stl->sel,i);
      } 

      stl->t = 0;
      stl->tn = 0;
      nbmodel_scores(fmt, m, stl->sel, ess, dt,
		     &stl->best_score_zos, &stl->best_score_lgs);
      stl->tn_when_best = 0;
      stl->hit_rate_rwa = 0;
    }
  }

  search_stats_free(stl);
  search_stats_free(stg);
  search_stats_free(stp);
  score_hashtable_free(sht);
  score_hashtable_free(sht2);

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
  greedy_search(fmt, m, ess, dt, argv[6], argv[7]);

  nbmodel_free(fmt, &m);
  format_free(fmt);
  data_free(dt);

  return 0;
}

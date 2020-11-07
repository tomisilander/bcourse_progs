#include <stdio.h>
#include <math.h>
#include <signal.h>
#include <errno.h>
#include <unistd.h>
#include "err.h"
#include "format.h"
#include "data.h"
#include "parent_matrix.h"
#include "score_hashtable.h"
#include "banesearch.h"
#include "bane.h"
#include "node.h"

int usr1_set;
int usr2_set;
int alrm_set;

void usr2_handler(int signum) {
  signal(12,usr2_handler);
  signum = 0;
  usr2_set = 1;
}

void usr1_handler(int signum) {
  signal(10,usr1_handler);
  signum = 0;
  usr1_set = 1;
}

void alrm_handler(int signum) {
  signal(14,alrm_handler);
  signum = 0;
  alrm_set = 1;
}

void set_signal_handlers(){
  usr1_set = 0;
  usr2_set = 0;
  alrm_set = 0;
  signal(10,usr1_handler);
  signal(12,usr2_handler);
  signal(14,alrm_handler);
}

typedef struct search_stats search_stats;

struct search_stats {
  unsigned long t;
  unsigned long tn;
  double best_score;
  unsigned long tn_when_best;
  double hit_rate_rwa;
  bane* beba;
};

search_stats* search_stats_create(bane* bn) {
  search_stats* me;
  MECALL(me,1,search_stats);
  me->beba = bane_copy(bn);
  return me;
}

void search_stats_free(search_stats* me) {
  if(me->beba != NULL) {
    bane_free(me->beba);
  }
  free(me);
}

void search_stats_assign(search_stats* dst, search_stats* src) {
  dst->t = src->t;
  dst->tn = src->tn;
  dst->best_score = src->best_score;
  dst->tn_when_best = src->tn_when_best;
  dst->hit_rate_rwa = src->hit_rate_rwa;
  bane_assign(dst->beba, src->beba);
}

void write_report_n_structure(search_stats* st, search_stats* stp,
			      char* reportfilename, char* structfilename){
  FILE* fp;
  double better = 0;

  /* Structure first since it is bigger */

  OPENFILE_OR_DIE(fp,structfilename,"w");
  bane_write_structure(st->beba,fp);
  CLOSEFILE_OR_DIE(fp,structfilename);

  /* Report second */

  OPENFILE_OR_DIE(fp,reportfilename,"w");
  

  better = (stp->t == 0) ? 
    0 : exp(st->best_score - stp->best_score);

  fprintf(fp,"%lu %lu %lu %g %f\n", 
	  st->tn, 
	  st->tn - stp->tn, 
	  st->tn - st->tn_when_best, 
	  better, 
	  st->best_score);

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


typedef int (*greedy_change)(bane*, arc*, int, int);
typedef void (*greedy_unchange)(bane*, arc*, int);

static
void remove_bad_arcs(data* dt, double ess, 
		     arc* ar, double* scoreboard, score_hashtable* sht,
		     search_stats* stg){

  int all_removed = 0; 
  bane* bn = stg->beba;

  bane_gather_full_ss(bn, dt);
  bane_get_score(bn, ess, scoreboard);

  while(!all_removed){
    all_removed = 1;
    
    for(ar->to=0; ar->to<bn->nodecount; ++ar->to){
      node* nodi = bn->nodes + ar->to;
      for(ar->from = 0; ar->from<bn->nodecount; ++ar->from){
	if (IS_PARENT(ar->to,ar->from,bn->pmx)){
	  double old_score;
	  double gain = 0;	  
	  bane_del_arc(bn, ar, 1);
	  if(!score_hashtable_get(sht, bn->pmx->mx, &old_score)) {
	    bane_gather_ss_for_i(bn, dt, ar->to);
	    gain = bane_get_score_for_i(nodi,ess) - scoreboard[ar->to];
	    if(gain >= 0) {
	      scoreboard[ar->to] += gain;
	      stg->tn_when_best = stg->tn;
	      stg->best_score += gain;
	      all_removed = 0;	      
	    } else {
	      bane_add_arc(bn, ar, 1);
	    }
	    ++ stg->tn;
	  } else {
	    bane_add_arc(bn, ar, 1);	  
	  }
	}
      }
    }
  }
}

static
void t_times_greedy_step(data* dt, double ess, int t_limit, int maxtblsize,
			 arc* ar, double* scoreboard, score_hashtable* sht,
			 search_stats* stl, search_stats* stg){

  int t; 
  bane* bn = stl->beba;

  greedy_change grch[3];
  greedy_unchange gruch[3];

  grch[0] = bane_add_random_arc;
  grch[1] = bane_del_random_arc;
  grch[2] = bane_rev_random_arc;

  gruch[0] = bane_del_arc;
  gruch[1] = bane_add_arc;
  gruch[2] = bane_rev_arc;

  for(t = 0; (t<t_limit) && (!usr1_set) && (!usr2_set) &&(!alrm_set); ++t){
    int ch = rand()%3;
    double new_score;
    double from_score = 0;
    double to_score = 0;

    int hthit;

    if(!grch[ch](bn,ar,maxtblsize,1)) continue; /* if CHANGE not succeeded */

    hthit = 0;
    if(score_hashtable_get(sht, bn->pmx->mx, &new_score)) {
      hthit = 1;
    } else {
      bane_gather_ss_for_i(bn, dt, ar->to); /* to needs always work */
      to_score = bane_get_score_for_i(bn->nodes + ar->to, ess);	
      from_score = scoreboard[ar->from]; /* cheap default for from */
      if(ch==2) { /* rev - from needs to be calculated */
	bane_gather_ss_for_i(bn, dt, ar->from);
	from_score = bane_get_score_for_i(bn->nodes + ar->from, ess);
      }
      new_score = stl->best_score +  
	from_score + to_score 
      - scoreboard[ar->from] - scoreboard[ar->to];
      score_hashtable_put(sht, bn->pmx->mx, new_score);
      
      ++ stl->tn;
      ++ stg->tn;
    }

    stl->hit_rate_rwa += 0.01 * (hthit - stl->hit_rate_rwa);

    ++ stl->t;
    ++ stg->t;

    if(new_score > stl->best_score) { /* commit */
      /* fprintf(stderr,"%d %f\n",ch,new_score); */
      stl->best_score = new_score;
      stl->tn_when_best = stl->tn;
      if(hthit){ /* We need to update score board */ 
	bane_gather_ss_for_i(bn, dt, ar->to); /* to needs always work */
	to_score = bane_get_score_for_i(bn->nodes + ar->to, ess);	
	from_score = scoreboard[ar->from]; /* cheap default for from */
	if(ch==2) { /* rev - from needs to be calculated */
	  bane_gather_ss_for_i(bn, dt, ar->from);
	  from_score = bane_get_score_for_i(bn->nodes + ar->from, ess);
	}
      }
      scoreboard[ar->from] = from_score;
      scoreboard[ar->to] = to_score;
      if(new_score > stg->best_score) {
	stg->best_score = new_score;
	stg->tn_when_best = stg->tn;
	bane_assign(stg->beba, bn);
      }
    } else { /* retract */
      gruch[ch](bn, ar, 1);
    }
  }
}

void greedy_search(format* fmt, data* dt, double ess, int maxtblsize,
		   char* reportfilename, char* structfilename){

  search_stats* stg; /* global search stats */
  search_stats* stl; /* search stats for current round - local*/
  search_stats* stp; /* search stats at previous report */

  double* scoreboard;
  score_hashtable* sht;
  int keysize;
  int items_per_pos;
  arc* ar;

  bane* bn = NULL;

  /* Load the initial net */
  {  
    FILE* fp;
    bn = bane_create_from_format(fmt);
    OPENFILE_OR_DIE(fp,structfilename,"r");
    bane_read_structure(bn,fp);
    CLOSEFILE_OR_DIE(fp,structfilename);
  }

  stg = search_stats_create(bn);
  stl = search_stats_create(bn);
  stp = search_stats_create(bn);

  bane_free(bn);

  MECALL(scoreboard, stl->beba->nodecount, double);

  keysize = stl->beba->pmx->m * stl->beba->pmx->one_dim_size;
  items_per_pos = score_hashtable_items_for_mem(5000000, keysize);
  if(items_per_pos < 1) items_per_pos = 1;
  if(items_per_pos > 10) items_per_pos = 10;
  sht = score_hashtable_create(keysize, items_per_pos);

  bane_gather_full_ss(stl->beba, dt);
  stg->best_score = stl->best_score = 
    bane_get_score(stl->beba, ess, scoreboard);
  score_hashtable_put(sht,stl->beba->pmx->mx, stg->best_score); 
  stl->hit_rate_rwa = 0;

  ++stl->t;
  ++stg->t;
  ++stl->tn;
  ++stg->tn;

  stp->t = -1;
  stp->tn = -1;
  stp->best_score = stg->best_score;

  MECALL(ar, 1, arc); 

  while(!usr2_set){

    t_times_greedy_step(dt, ess, 100, maxtblsize, ar, scoreboard, sht,stl,stg);

    if(usr1_set || alrm_set) {
      usr1_set = 0;
      if(alrm_set) {
	remove_bad_arcs(dt, ess, ar, scoreboard, sht, stg);
      }
      write_report_n_structure(stg, stp, reportfilename, structfilename);
      if(alrm_set){
	usr2_set = 1;
      }
    }

    if(better_to_restart(dt,stl)){
      int n;
      int t;
      parent_matrix* smx;
      /* Take some good one out of hashtable and mutilate it */
      
      n = 1 + rand() % (1 + (int) (sht->keycount * 0.01));
      smx = parent_matrix_create_wrap(stl->beba->nodecount,
				      score_hashtable_get_nth_key(sht,n));
      bn = bane_create_from_pmx(fmt,smx);
      parent_matrix_free(smx);
      
      for(t=0; t<bn->nodecount / 3; ++t) {
	bane_del_random_arc(bn,ar,maxtblsize,1);
      } 

      bane_assign(stl->beba, bn);
      bane_free(bn);
      
      bane_gather_full_ss(stl->beba, dt);
      stl->t = 0;
      stl->tn = 0;
      stl->best_score = bane_get_score(stl->beba, ess, scoreboard);
      stl->tn_when_best = 0;
      stl->hit_rate_rwa = 0;
    }
  }

  search_stats_free(stl);
  search_stats_free(stg);
  search_stats_free(stp);
  free(scoreboard);
  free(ar);
  score_hashtable_free(sht);

}

int main(int argc, char* argv[]){
  int n;
  data* dt;
  format* fmt;
  double ess;
  int searchtime;
  FILE* fp;
  int maxtblsize;

  if((argc !=9) && (argc!=10)){
    fprintf(stderr, 
	    "Usage: %s vdfile datafile datacount ess "
	    "reportfile structfile searchtime(<0 forever) maxtblsize"
	    "[pidfile]\n", 
	    argv[0]);
    exit(-1);
  }
   
  set_signal_handlers();

  if(argc==9) {
    fprintf(stdout,"%d\n",getpid());
    fclose(stdout);
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

  ess = atof(argv[4]);

  searchtime = atoi(argv[7]);
  if(searchtime >= 0) { /* HEY, 0 WILL NOT DO ANYTHING */
    alarm(searchtime);
  }
  maxtblsize = atoi(argv[8]);

  greedy_search(fmt, dt, ess, maxtblsize, argv[5], argv[6]);

  format_free(fmt);
  data_free(dt);
  
  return 0;
}

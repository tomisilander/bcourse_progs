#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "err.h"
#include "format.h"
#include "data.h"
#include "bane.h"
#include "node.h"

typedef struct arcwstruct arcwstruct;

struct arcwstruct {
  int from;
  int to;
  double weight;
};

int arcsort(const void* a1, const void* a2){
  if(((arcwstruct*) a1)->weight > ((arcwstruct*) a2)->weight)
    return -1;
  else 
    return ((arcwstruct*) a1)->weight < ((arcwstruct*) a2)->weight;
}

int main(int argc, char* argv[]){

  format* fmt;
  data *dt;
  bane* bn;
  double ess;
  FILE* fp;
  double refscore;
  int arccount;
  arcwstruct* arws;
  arc* ar;
  double* scoreboard;

  char* type_mx;
  double* gain_mx;
  int* rank_mx;

  int poscount = 0;
  int zercount = 0;
  int negcount = 0;

  int i;
  
  if(argc != 6){
    fprintf(stderr, 
	    "Usage: %s vdfile datafile datacount structfile ess\n", 
	    argv[0]);
    exit(-1);
  }

  fmt = format_cread(argv[1]);
  dt = data_cread(atoi(argv[3]), fmt->dim, argv[2]);
  bn = bane_create_from_format(fmt);
  OPENFILE_OR_DIE(fp,argv[4],"r"); 
  bane_read_structure(bn,fp);
  CLOSEFILE_OR_DIE(fp,argv[4]);
  ess = atof(argv[5]);

  arccount = bn->nodecount * bn->nodecount;

  MEMALL(arws, arccount, arcwstruct);
  MECALL(scoreboard, bn->nodecount, double);
  MEMALL(ar, 1, arc);
  MECALL(type_mx, arccount, char);
  MECALL(gain_mx, arccount, double);
  MECALL(rank_mx, arccount, int);

  refscore = BDe_score(bn,dt,ess,scoreboard);

  arccount = 0;
  for(i=0; i<bn->nodecount; ++i){
    int p;
    ar->to = i;
    for(p=0; p<bn->nodecount; ++p){
      double posweight = 0;
      char tc = 'X';
      int ix = i * bn->nodecount + p;

      ar->from = p;
      arws[arccount].from = p;
      arws[arccount].to = i;

      if(IS_PARENT(i,p,bn->pmx)){
	tc = '1';
	bane_del_arc(bn, ar, 0);
	posweight = scoreboard[i] - BDe_score_for_v(bn,dt,i,ess);
	bane_add_arc(bn, ar, 0);
      } else if((i!=p) && (bn->nodes[p].path_to_me_count[i]==0)) {
	tc = '0';
	bane_add_arc(bn, ar, 0);
	posweight = scoreboard[i] - - BDe_score_for_v(bn,dt,i,ess);
	bane_del_arc(bn, ar, 0);
      }

      type_mx[ix] = tc;

      if(tc != 'X') { /* Weight is current/change */
	gain_mx[ix] = arws[arccount++].weight = posweight;
      }
    }
  }

  qsort(arws, arccount, sizeof(arcwstruct), arcsort);

  { /* Fill rank matrix entries */

    int a = 0;

    for(a=0; a<arccount; ++a){
      int rnk = -1;
      int ix = arws[a].to * bn->nodecount + arws[a].from;

      if(gain_mx[ix] > 0)      { rnk = poscount++; } 
      else if(gain_mx[ix] < 0) { rnk = negcount++; } 
      else                     { rnk = zercount++; }
      
      rank_mx[ix] = rnk;
    }
  }
  

  /* Print scores */

  printf("%g\n", refscore);
  for(i=0; i<bn->nodecount; ++i){
    printf("%g%c", scoreboard[i], (i+1 == bn->nodecount) ? '\n' : ' ');
  }
  printf("\n");


  /* Print types */

  for(i=0; i<bn->nodecount; ++i){
    int p;
    for(p=0; p<bn->nodecount; ++p){
      int ix = i * bn->nodecount + p;
      printf("%c%c", type_mx[ix], (p+1 == bn->nodecount) ? '\n' : ' ');
    }
  }
  printf("\n");


  /* Print gains */

  for(i=0; i<bn->nodecount; ++i){
    int p;
    for(p=0; p<bn->nodecount; ++p){
      int ix = i * bn->nodecount + p;
      printf("%g%c", 
	     (type_mx[ix] == 'X') ? 0 : gain_mx[ix], 
	     (p+1 == bn->nodecount) ? '\n' : ' ');
    }
  }
  printf("\n");


  /* Print ranks */

  printf("%d %d %d\n", poscount, zercount, negcount);
  for(i=0; i<bn->nodecount; ++i){
    int p;
    for(p=0; p<bn->nodecount; ++p){
      int ix = i * bn->nodecount + p;
      printf("%d%c", 
	     (type_mx[ix] == 'X') ? -1 : rank_mx[ix], 
	     (p+1 == bn->nodecount) ? '\n' : ' ');
    }
  }
  printf("\n");


  /* Free world **/

  free(gain_mx);
  free(type_mx);
  free(rank_mx);

  free(scoreboard);
  free(ar);

  free(arws);
  format_free(fmt);
  data_free(dt);
  bane_free(bn);

  return 0;
}

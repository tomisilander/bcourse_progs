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

  arccount = 0;
  for(i=0; i<bn->nodecount; ++i){
    arccount += bn->nodes[i].parentcount;
  }

  MEMALL(arws, arccount, arcwstruct);
  MECALL(scoreboard, bn->nodecount, double);
  MEMALL(ar, 1, arc);

  refscore = BDe_score(bn,dt,ess,scoreboard);
  
  arccount = 0;
  for(i=0; i<bn->nodecount; ++i){
    Word_t p;
    node* nodi = bn->nodes + i;
    ar->to = i;
    LOOP_OVER_PARENTS(nodi,p){
      ar->from = p;
      bane_del_arc(bn, ar, 0);
      arws[arccount].from = p;
      arws[arccount].to = i;
      arws[arccount].weight = scoreboard[i] - BDe_score_for_v(bn,dt,i,ess);
      ++arccount;
      bane_add_arc(bn, ar, 0);
    } END_PARENT_LOOP;
  }

  free(scoreboard);
  free(ar);

  qsort(arws, arccount, sizeof(arcwstruct), arcsort);

  {
    int a;
    for(a=0; a<arccount; ++a)    
      printf("%d %d %g\n", arws[a].from, arws[a].to, exp(arws[a].weight));
  }
  
  free(arws);
  format_free(fmt);
  data_free(dt);
  bane_free(bn);

  return 0;
}

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "err.h"
#include "format.h"
#include "data.h"
#include "bitset.h"
#include "nb_selector.h"
#include "nbmodel.h"

typedef struct arcwstruct arcwstruct;

struct arcwstruct {
  int from;
  int to;
  double weight;
  double weight2;
};

int arcsort(const void* a1, const void* a2){
  if(((arcwstruct*) a1)->weight > ((arcwstruct*) a2)->weight)
    return -1;
  else if(((arcwstruct*) a1)->weight < ((arcwstruct*) a2)->weight)
     return 1;
  else if(((arcwstruct*) a1)->weight2 > ((arcwstruct*) a2)->weight2)
    return -1;
  else
    return ((arcwstruct*) a1)->weight2 < ((arcwstruct*) a2)->weight2;
}

int main(int argc, char* argv[]){

  format* fmt;
  data *dt;
  double refscore;
  double refscore2;
  int arccount;
  arcwstruct* arws;
  int i, classix;
  double ess;
  nbmodel* m;
  int* sel;

  if(argc != 7){
    fprintf(stderr, 
	    "Usage: %s vdfile datafile datacount ess classix strfile\n", 
	    argv[0]);
    exit(-1);
  }

  fmt = format_cread(argv[1]);
  dt = data_create(atoi(argv[3]), fmt->dim);
  data_read(argv[2], dt);
  ess = atof(argv[4]);
  classix = atoi(argv[5]);
  m = nbmodel_create(fmt, dt, classix);
  sel = nb_selector_cread(argv[6]);

  arccount = 0;
  for(i=0; i<fmt->dim; ++i){
      if ((i!=m->root) && (BITSET_ON(sel,i))) ++ arccount;
  }

  MEMALL(arws, arccount, arcwstruct);

  {
    int i;

    nbmodel_scores(fmt,m,sel, ess,dt, &refscore,&refscore2);
    
    arccount = 0;
    
    for(i=0; i<fmt->dim; ++i){
      double zos = 0;
      double lgs = 0;
      
      if ((i==m->root) || (BITSET_OFF(sel,i))) continue;

      BITSET_UNSET(sel,i);

      nbmodel_scores(fmt, m, sel, ess, dt, &zos,&lgs);
      
      arws[arccount].from = m->root;
      arws[arccount].to   = i;
      arws[arccount].weight  = refscore  - zos;
      arws[arccount].weight2 = refscore2 - lgs;
 
      ++arccount;      

      BITSET_SET(sel,i);
    }
  }

  qsort(arws, arccount, sizeof(arcwstruct), arcsort);

  {
    int a;
    for(a=0; a<arccount; ++a)    
      printf("%d %d %g %g\n", 
	     arws[a].from, arws[a].to, 
	     100*arws[a].weight, exp(arws[a].weight2));
  }
  
  free(sel);
  free(arws);
  nbmodel_free(fmt,&m);
  format_free(fmt);
  data_free(dt);

  return 0;
}

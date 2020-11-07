#ifndef FOREST_H
#define FOREST_H

#include "data.h"
#include "format.h"
#include "parent_matrix.h"

typedef struct forest forest;

struct forest {
  int dim;
  int* valcount;
  int* parent;
  double ess;

  double** th;

  /* private */
  double* dst;
  char* connected;
  double* orphan_lmlh;
  int* orphan_ss_tmp;
  int* ss_tmp;
};

extern forest* 
forest_create(format* fmt, double ess);

extern void forest_free(forest* frt);

extern double 
forest_learn(forest* frt, data* dt);

extern double 
forest_learn_no_arcs(forest* frt, data* dt);

extern void 
forest_update_parameters(forest* frt, data* dt);

extern void 
forest_set_random_parameters(forest* frt);

extern 
double forest_d1_lprob(forest* frt, data* dt, int j);

extern parent_matrix* 
forest_get_parent_matrix(forest* frt);

extern void 
forest_write_structure(forest* frt, FILE* fp);

extern void 
forest_write_params(forest* frt, FILE* fp);

extern void 
forest_write(forest* frt, FILE* fp);

#define FOREST_SS(F,P,J,C,K) ((F)->ss_tmp[(J) * (F)->valcount[C] + (K)])
#define FOREST_TH(F,P,J,C,K) ((F)->th[C][(J) * (F)->valcount[C] + (K)])

#define FOREST_LP(LP,F,DT,J) {\
  int i;\
  (LP) = 0;\
  for(i=0; i<(F)->dim; ++i)\
    (LP) += log(FOREST_TH((F),0,((F)->parent[i] == -1) ? 0 : D((DT),(J),(F)->parent[i]),i,D((DT),(J),i)));\
}\

#endif


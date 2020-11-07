#ifndef BANE_H
#define BANE_H

#include <stdio.h>
#include <Judy.h>
#include "typedef_node.h"
#include "parent_matrix.h"
#include "format.h"
#include "data.h"

typedef struct bane bane;

struct bane {

  int nodecount;
  node* nodes;
  
  parent_matrix* pmx;

};

typedef struct arc arc;  
struct arc {
  int from;
  int to;
};

extern bane* 
bane_cread(char* filename);

extern void 
bane_free(bane* bn);

extern bane* 
bane_create(int nodecount);

extern bane* 
bane_create_from_pmx(format* fmt, parent_matrix* pmx);

extern bane* 
bane_create_from_format(format* fmt);

extern bane* 
bane_copy(bane* src);

extern void
bane_copy_ss(bane* dst, bane* src);

extern void 
bane_assign(bane* dst, bane* src);

extern bane* 
bane_create_forest(format* fmt, double ess, data* dt);

extern void 
bane_add_arc(bane* bn, arc* ar, int commit);

extern void 
bane_del_arc(bane* bn, arc* ar, int commit);

extern void 
bane_rev_arc(bane* bn, arc* ar, int commit);

extern void 
bane_fix_add_arc(bane* bn, arc* ar);

extern void 
bane_fix_del_arc(bane* bn, arc* ar);

extern void 
bane_fix_rev_arc(bane* bn, arc* ar);

extern void bane_clear_ss(bane* bn);

extern void 
bane_add_ss(bane* bn, data* dt, int j);

extern void 
bane_del_ss(bane* bn, data* dt, int j);

extern 
void bane_gather_ss_for_i(bane* bn, data* dt, int i);

extern void 
bane_gather_full_ss(bane* bn, data* dt);

extern void 
bane_gather_full_ss_in_order(bane* bn, data* dt);

extern void 
bane_gather_full_ss_in_order_from_file(bane* bn, char* filename);

extern double 
bane_logprob(bane* bn, data* dt, int j, double ess, double* terms);

extern Pvoid_t 
gather_sparse_ss_for_v(bane* bn, data* dt, int v);

extern void 
free_sparse_ss(Pvoid_t ss);

extern double 
BDe_score_for_v_from_ss(bane* bn, Pvoid_t ss, int v, double ess);

extern double 
BDe_score_for_v(bane* bn, data* dt, int v, double ess);

extern double 
BDe_score(bane* bn, data* dt, double ess, double* scoreboard);

extern double 
bane_get_score_for_i(node* nodi, double ess);

extern double 
bane_get_score(bane* bn, double ess, double* scoreboard);

extern void
bane_write_structure(bane* bn, FILE* fp);

extern bane* 
bane_cread_structure(FILE* fp);

extern void
bane_read_structure(bane* bn, FILE* fp);

extern void
bane_write_with_ss(bane* bn, char* filename);

extern bane*
bane_read_with_ss(char* filename);

#endif





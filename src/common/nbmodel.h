#ifndef NBMODEL_H
#define NBMODEL_H

struct nbmodel_{
  int* h;
  int** hh;
  int*** f;
  int N;
  int root;
};

typedef struct nbmodel_ nbmodel;

#include "format.h"
#include "data.h"

extern nbmodel* 
nbmodel_create(format* fmt, data* dt, int root);

extern void 
nbmodel_free(format* fmt, nbmodel** m);

extern int 
nbmodel_predict(format* fmt, nbmodel* m, int* selector, double ess,
		data* dt, int j, double* dstr);

extern void 
nbmodel_scores(format* fmt, nbmodel* m, int* selector, double ess, data* dt, 
	       double* zosp, double* lgsp);

extern void 
nbmodel_add_ss(format* fmt, nbmodel* m, data* dt, int j);

extern void 
nbmodel_del_ss(format* fmt, nbmodel* m, data* dt, int j);

extern void 
nbmodel_gather_ss_i(format* fmt, nbmodel* m, data* dt, int i);

extern void 
nbmodel_write_ss(format* fmt, nbmodel* m, char* fn);

#endif













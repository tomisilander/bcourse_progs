#ifndef BANE_NB_H
#define BANE_NB_H

#include "bane.h"

extern int 
bane_nb_predict(bane* bn, int classix, data* dt, int j, double ess, 
		double* distr);

void bane_nb_loo_score(bane* bn, data *dt, int classix, double ess,
		       double* zosp, double* lgsp);

#endif

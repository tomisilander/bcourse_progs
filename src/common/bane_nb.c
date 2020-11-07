#include <math.h>
#include "err.h"
#include "node.h"
#include "bane.h"
#include "bane_nb.h"

int bane_nb_predict(bane* bn, int classix, data* dt, int j, double ess, 
		    double* distr){
  int maxk, k;
  int K = bn->nodes[classix].valcount; 

  double nrm = 0;
  for(k=0; k<K; ++k){
    int i;
    distr[k]= (SS(bn->nodes+classix,0,k) + ess) 
      / (bn->nodes[classix].N + K * ess);
    for(i=0; i<bn->nodecount; ++i){
      node* nodi = bn->nodes + i;
      if((i==classix) || (nodi->parentcount == 0) || MISSING(dt,j,i)) 
	continue;
      distr[k] *= (SS(nodi,k,D(dt, j, i)) + ess )  
	/ (nodi->N + nodi->valcount * ess);
    }
    nrm += distr[k];
  }

  maxk = -1;
  for(k=0; k<K; ++k){
    distr[k] /= nrm;
    if(k==0 || distr[k] > distr[maxk]) maxk = k;
  }

  return maxk;
}


void bane_nb_loo_score(bane* bn, data *dt, int classix, double ess,
		       double* zosp, double* lgsp) {
  int j;
  double* distr;

  MEMALL(distr,bn->nodes[classix].valcount,double);

  if(zosp) *zosp = 0;
  if(lgsp) *lgsp = 0;

  for(j=0;j<dt->N;++j) {
    int maxk;
    int cork = D(dt,j,classix);

    bane_del_ss(bn,dt,j);
    maxk = bane_nb_predict(bn,classix, dt, j, ess, distr);
    bane_add_ss(bn,dt,j);
      
    if(zosp && (cork == maxk)) ++ *zosp;
    if(lgsp) *lgsp += log(distr[cork]);
  }

  free(distr);

  if(zosp) *zosp /= dt->N;
  if(lgsp) *lgsp /= dt->N;
}


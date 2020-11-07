#include <math.h>
#include "err.h"
#include "bitset.h"
#include "nbmodel.h"

nbmodel* nbmodel_create(format* fmt, data* dt, int root){
  int i,j,k;
  int K = fmt->valcount[root];
  nbmodel* m;
  MEMALL(m, 1,nbmodel);

  MECALL(m->h, K, int);
  MEMALL(m->hh, K, int*);
  m->root = root;

  MEMALL(m->f, K, int**);
  for(k=0;k<K;++k){    
    MECALL(m->hh[k], fmt->dim, int);
    MEMALL(m->f[k],fmt->dim, int*);
    for(i=0;i<fmt->dim;++i){
      MECALL(m->f[k][i], fmt->valcount[i], int);
    }
  }
  

  /* Gather ss */

  m->N=0;
  for(j=0;j<dt->N;++j) nbmodel_add_ss(fmt,m,dt,j);

  return m;
}

void nbmodel_add_ss(format* fmt, nbmodel* m,
		    data* dt, int j){
  int k, i;

  if(MISSING(dt,j,m->root)) return;
    
  k = D(dt,j,m->root);
  
  ++(m->h[k]);
  for(i=0;i<fmt->dim;++i){
    if(MISSING(dt,j,i)) continue;
    ++(m->f[k][i][D(dt,j,i)]);
    ++(m->hh[k][i]);
  }
  ++m->N;
}

void nbmodel_del_ss(format* fmt, nbmodel* m,
		    data* dt, int j){
  int k, i;

  if(MISSING(dt,j,m->root)) return;
    
  k = D(dt,j,m->root);
  
  --(m->h[k]);
  for(i=0;i<fmt->dim;++i){
    if(MISSING(dt,j,i)) continue;
    --(m->f[k][i][D(dt,j,i)]);
    --(m->hh[k][i]);
  }
  --m->N;
}


void nbmodel_gather_ss_i(format* fmt, nbmodel* m, data* dt, int i){
  int j, k;

  if(i==m->root) return;

  for(k=0; k<fmt->valcount[m->root]; ++k){
    int l;
    m->hh[k][i] = 0;
    for(l=0; l<fmt->valcount[i]; ++l) m->f[k][i][l] = 0;
  }

  for(j=0;j<dt->N;++j) {
    if(MISSING(dt,j,m->root)) continue;
    k = D(dt,j,m->root);
    if(MISSING(dt,j,i)) continue;
    ++(m->f[k][i][D(dt,j,i)]);
    ++(m->hh[k][i]);
  }

}


void nbmodel_free(format* fmt, nbmodel** m){
  int k,i;
  int K=fmt->valcount[(*m)->root];
  free((*m)->h);
  for(k=0;k<K;++k){    
    for(i=0;i<fmt->dim;++i){
      free((*m)->f[k][i]);
    }
    free((*m)->f[k]);
    free((*m)->hh[k]);
  }
  free((*m)->f);
  free((*m)->hh);
  free((*m));
  (*m)=NULL;
}

int nbmodel_predict(format* fmt, nbmodel* m, int* selector, double ess,
		    data* dt, int j, double* distr){
  int maxk, k;
  int K = fmt->valcount[m->root];  

  double nrm = 0;
  for(k=0; k<K; ++k){
    int i;
    distr[k]=(m->h[k] + ess) / (m->N + K*ess);
    for(i=0; i<fmt->dim; ++i){
      if((i==m->root) 
	 || ((selector != NULL) && BITSET_OFF(selector,i))
	 || MISSING(dt,j,i)) continue;
      distr[k] *= (m->f[k][i][D(dt, j, i)] + ess)  
	/ (m->hh[k][i] + fmt->valcount[i] * ess);
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

void 
nbmodel_scores(format* fmt, nbmodel* m, int* selector, double ess, data* dt, 
	       double* zosp, double* lgsp) {
  int j;
  double* distr;

  MEMALL(distr,fmt->valcount[m->root],double);

  if(zosp) *zosp = 0;
  if(lgsp) *lgsp = 0;

  for(j=0;j<dt->N;++j) {
    int maxk;
    int cork = D(dt,j,m->root);
    nbmodel_del_ss(fmt, m, dt, j);
    maxk = nbmodel_predict(fmt, m, selector, ess, dt, j, distr);
    nbmodel_add_ss(fmt, m, dt, j);
      
    if(zosp && (cork == maxk)) ++ *zosp;
    if(lgsp) *lgsp += log(distr[cork]);
  }

  free(distr);

  if(zosp) *zosp /= dt->N;
  if(lgsp) *lgsp /= dt->N;
}


void nbmodel_write_ss(format* fmt, nbmodel* m, char* fn) {
  FILE* fp = fopen(fn,"w");
  int k,i,l;
  for(k=0;k<fmt->valcount[m->root];++k){
    fprintf(fp,"h[%d]=%d\n",k,m->h[k]);
    for(i=0;i <fmt->dim; ++i){
      fprintf(fp,"hh[%d][%i]=%d\n",k,i,m->hh[k][i]);
      for(l=0;l<fmt->valcount[i];++l){
	fprintf(fp,"f[%d][%d][%d]=%d ",k,i,l,m->f[k][i][l]);
      }
      fprintf(fp,"\n");
    }
  }   
  fclose(fp);
}

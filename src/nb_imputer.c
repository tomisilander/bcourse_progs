#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "err.h"
#include "data.h"
#include "format.h"
#include "nbmodel.h"

typedef struct ment ment; 
struct ment {
   int j;
   int i; 
};

int gather_ments(format* fmt, data* dt, ment* mentab){
  int res = 0;
  int j;
  for(j=0;j<dt->N;++j){
    int i;
    for(i=0;i<fmt->dim;++i){
      if(MISSING(dt,j,i)) {
	if(mentab!=NULL){
	  mentab[res].i=i;
	  mentab[res].j=j;
	}
	++res;
      }
    }
  }
  return res;
}

void impute_surest_ment(format* fmt, data* dt, 
			ment* mentab, int mtl,
			nbmodel** modtab, double** dstrtab){
  int sme = -1;
  int smek = -1;
  double sment = 0;

  int me;
  int mi;
  nbmodel* m = NULL;
  int K=-1;
  double* distr = NULL;
  
  /* Find surest me */

  for(me=0;me<mtl; ++me){
    int k;
    int maxk;
    double mentropy;
    /* double nrm; */

    ment* mnt = mentab + me;

    if(mnt->i<0) continue;

    K = fmt->valcount[mnt->i];

    if(modtab[mnt->i]==NULL){
      MEMALL(dstrtab[mnt->i], K, double);
      modtab[mnt->i] = nbmodel_create(fmt,dt,mnt->i);
    }

    m = modtab[mnt->i];
    distr = dstrtab[mnt->i];

    /* Predict me */

    maxk = nbmodel_predict(fmt,m,NULL,1,dt,mnt->j,distr);
#if 0
    nrm = 0;
    for(k=0; k<K; ++k){
      int i;
      distr[k]=(m->h[k]+1.0)/(m->N+K);
      for(i=0; i<fmt->dim; ++i){
	if((i == mnt->i) || MISSING(dt,mnt->j,i)) continue;
	distr[k] *= (m->f[k][i][D(dt, mnt->j, i)]+1.0)  
	  / (m->hh[k][i]+fmt->valcount[i]);
      }
      nrm += distr[k];
    }
#endif
    /* Get max and entropy */
    mentropy = 0;
    for(k=0; k<K; ++k){
      mentropy -= distr[k] * log(distr[k]);
    }

    if((sme==-1) || (mentropy < sment)){
      sme = me;
      sment = mentropy;
      smek = maxk;
    }
  }

  /* Impute surest me */

  D(dt, mentab[sme].j, mentab[sme].i) = smek;

  /* Update models */

  for(mi=0; mi<fmt->dim; ++mi) {
    m = modtab[mi];
    if(NULL==m) continue; /* No new models are created */
    if(MISSING(dt,mentab[sme].j,mi)) continue; /* If root for m is missing */

    if(mi == mentab[sme].i) { /* This model got new vector */
      int i;
      int k = smek;
      ++(m->h[k]);
      for(i=0;i<fmt->dim;++i){
	if(MISSING(dt,mentab[sme].j,i)) continue;
	++(m->f[k][i][D(dt,mentab[sme].j,i)]);
	++(m->hh[k][i]);
      }
      ++m->N;
    } else { /* these models got just a new variable */
      int k = D(dt,mentab[sme].j,mi);
      ++(m->f[k][mentab[sme].i][smek]);
      ++(m->hh[k][mentab[sme].i]);
    }
  }
  /*
  fprintf(stderr, "%d %i imputed with %d (%g)\n", 
	  mentab[sme].j, mentab[sme].i, smek, sment);
  */
  mentab[sme].i = -mentab[sme].j * fmt->dim - mentab[sme].i - 1;

}

void impute(format* fmt, data* dt){
  int me;
  int nof_md;
  ment* mentab;
  nbmodel** modtab;
  double** dstrtab;
  
  nof_md = gather_ments(fmt, dt, NULL);
  MEMALL(mentab,nof_md,ment);
  gather_ments(fmt, dt, mentab);
  MECALL(modtab, fmt->dim, nbmodel*);
  MECALL(dstrtab, fmt->dim, double*);

  /*
  fprintf(stderr, "%d Missing entries\n", nof_md);
  */

  for(me=0; me<nof_md; ++me) {
    impute_surest_ment(fmt,dt,mentab,nof_md,modtab,dstrtab);
  }
  {
    int i;
    for(i=0; i<fmt->dim; ++i) {
      if(NULL!=modtab[i]) nbmodel_free(fmt,modtab+i);
      if(NULL!=dstrtab[i]) free(dstrtab[i]);
    }
    free(modtab);
    free(dstrtab);
  }
  free(mentab);
}

int main(int argc, char **argv){
  format* fmt;
  data* dt;
  FILE* fp;
  int j,n;

  if(argc!=5) {
    fprintf(stderr, 
	    "Usage: %s vdfile datain datacount dataout\n",
	    argv[0]);
    exit(1);
  }

  fmt = format_cread(argv[1]);
  n = atoi(argv[3]);
  dt = data_create(n,fmt->dim);
  data_read(argv[2],dt);

  impute(fmt,dt);

  OPENFILE_OR_DIE(fp,argv[4],"w");
  for(j=0;j<dt->N;++j) {
    int i;
    for(i=0;i<fmt->dim;++i) {
      fprintf(fp,"%d%c", D(dt,j,i), ((i+1)==fmt->dim) ? '\n' : '\t');
    }
  }
  CLOSEFILE_OR_DIE(fp,argv[4]);                                               

  format_free(fmt);
  data_free(dt);
 
  return 0;

}


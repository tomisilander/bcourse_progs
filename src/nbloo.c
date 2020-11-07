#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "err.h"
#include "data.h"
#include "format.h"
#include "nbmodel.h"
#include "nb_selector.h"


int main(int argc, char **argv){
  format* fmt;
  data* dt;
  FILE* fpc;
  FILE* fpd;
  int j,n,classix;
  nbmodel* m;
  int* sel;
  double ess;

  if(argc!=9) {
    fprintf(stderr, 
	    "Usage: %s vdfile datain datacount ess classix structfile outscore outpred\n",
	    argv[0]);
    exit(1);
  }

  fmt = format_cread(argv[1]);
  n = atoi(argv[3]);
  dt = data_create(n,fmt->dim);
  data_read(argv[2],dt);
  ess = atof(argv[4]);
  classix = atoi(argv[5]);

  m = nbmodel_create(fmt,dt,classix);
  sel = nb_selector_cread(argv[6]);

  {

    int ka,kp;
    int K = fmt->valcount[classix];

    double lgs = 0;
    double zos = 0;
    double* distr;
    int** conf;
    double** lonf;

    MEMALL(distr,K,double);
    MEMALL(conf,K,int*);
    MEMALL(lonf,K,double*);
    for(ka=0; ka<K; ++ka) {
      MEMALL(conf[ka],K,int);
      MEMALL(lonf[ka],K,double);
      for(kp=0; kp<K; ++kp) {
	lonf[ka][kp] = conf[ka][kp] = 0;
      }
    }
    

    OPENFILE_OR_DIE(fpd,argv[8],"w");
    zos = lgs = 0;
    for(j=0;j<dt->N;++j) {
      int k, maxk;
      int cork = D(dt,j,classix);

      nbmodel_del_ss(fmt,m,dt,j);
      maxk = nbmodel_predict(fmt,m, sel, ess, dt, j, distr);
      nbmodel_add_ss(fmt,m,dt,j);
    
      if(cork == maxk) ++ zos;
      lgs += log(distr[cork]);
      
      ++ conf[cork][maxk];
      lonf[cork][maxk] += log(distr[cork]);
      
      fprintf(fpd,"%d\t%d",cork,maxk);
      for(k=0; k<K; ++k) {
	fprintf(fpd,"\t%f", distr[k]);
      }
      fprintf(fpd,"\n");
    }
    CLOSEFILE_OR_DIE(fpd,argv[8]);

    OPENFILE_OR_DIE(fpc,argv[7],"w");
    fprintf(fpc,"%f\t%f\n", zos/dt->N, lgs/dt->N);
    for(ka=0; ka<K; ++ka) {
      for(kp=0; kp<K; ++kp) {
	fprintf(fpc,"%d%c", conf[ka][kp],(kp==K-1)?'\n':'\t');
      }
    }

    for(ka=0; ka<K; ++ka) {
      for(kp=0; kp<K; ++kp) {
	fprintf(fpc,"%g%c", 
		(conf[ka][kp]) ? lonf[ka][kp]/conf[ka][kp] : 0,
		(kp==K-1)?'\n':'\t');
      }
    }
    CLOSEFILE_OR_DIE(fpc,argv[7]);

    free(distr);
    for(ka=0; ka<K; ++ka) {
      free(conf[ka]);
      free(lonf[ka]);
    }
    free(conf);
    free(lonf);
  }


  nbmodel_free(fmt,&m);
  free(sel);
  format_free(fmt);
  data_free(dt);

  return 0;

}


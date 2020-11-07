#include <stdlib.h>
#include <math.h>
#include "err.h"
#include "format.h"
#include "data.h"
#include "parent_matrix.h"
#include "forest.h"

void
theta_free(forest* frt) {
  int i;

  if(frt->th != NULL) {
    for(i=0; i<frt->dim; ++i){
      free(frt->th[i]);
    }
    free(frt->th);
    frt->th = NULL;
  }
}

void get_max_valcounts(format* fmt, 
		       int* max_valcount, int* max_valcount_2){
  int i;
  *max_valcount = *max_valcount_2 = 0;
  for(i=0; i<fmt->dim; ++i){
    if((i==0) || (fmt->valcount[i] > *max_valcount)){
      *max_valcount_2 = *max_valcount;
      *max_valcount = fmt->valcount[i];
    } else if (fmt->valcount[i] > *max_valcount_2){
      *max_valcount_2 = fmt->valcount[i];      
    }
  }
}

forest* forest_create(format* fmt, double ess){
  forest* frt;
  int i, max_valcount, max_valcount_2;
  
  MEMALL(frt, 1, forest);

  frt->dim = fmt->dim;

  MEMALL(frt->valcount, fmt->dim, int);
  for(i=0; i<frt->dim; ++i) frt->valcount[i] = fmt->valcount[i];

  MEMALL(frt->parent, fmt->dim, int);

  frt->ess = ess;
  frt->th = NULL;

  get_max_valcounts(fmt, &max_valcount, &max_valcount_2);

  MEMALL(frt->dst, frt->dim, double);
  MEMALL(frt->orphan_lmlh, frt->dim, double);
  MEMALL(frt->connected, frt->dim, char);
  MEMALL(frt->orphan_ss_tmp, max_valcount, int);
  MEMALL(frt->ss_tmp, max_valcount * max_valcount_2, int);

  return frt;
}

void forest_free(forest* frt){
  theta_free(frt);
  free(frt->valcount);
  free(frt->parent);
  free(frt->dst);
  free(frt->orphan_lmlh);
  free(frt->connected);
  free(frt->orphan_ss_tmp);
  free(frt->ss_tmp);

  free(frt);
}

extern double lgamma(double);

void get_orphan_lmlh(forest* frt, data* dt){
  int i;
  double first_term = lgamma(frt->ess) - lgamma(frt->ess + dt->N);

  for(i=0; i<frt->dim; ++i){
    int j,l;
    int vci = frt->valcount[i];

    if(vci == 1){
      frt->orphan_lmlh[i] = 0;
      continue;
    }
    /* clear ss */
    for(l=0; l<vci; ++l) frt->orphan_ss_tmp[l] = 0;
    
    /* collect ss */
    for(j=0; j<dt->N; ++j){
      if(!MISSING(dt,j,i)) ++(frt->orphan_ss_tmp[D(dt,j,i)]);
    }
    
    /* update orphan_lmlh */ 
    frt->orphan_lmlh[i]  = first_term - vci * lgamma(frt->ess / vci);
    for(l=0; l<vci; ++l){
      frt->orphan_lmlh[i] += lgamma(frt->orphan_ss_tmp[l] + frt->ess / vci);
    }
  }
}


double lmlh_term(forest* frt, data* dt, int c, int p) {
  double bde = 0;
  int pcc = frt->valcount[p];
  int vcc = frt->valcount[c];
  double aij         = frt->ess / pcc;
  double aijk        = aij / vcc;
  int j;

  if((pcc == 1) && (vcc ==1)) return 0;

  /* clear ss */
  for(j = 0; j < pcc; ++j){
    int k;
    for(k = 0; k < vcc; ++k){
      FOREST_SS(frt,p,j,c,k) = 0;
    }
  }

  /* collect ss */
  for(j=0; j<dt->N; ++j){
    if(MISSING(dt,j,p) || MISSING(dt,j,c)) continue;
    ++ (FOREST_SS(frt,p,D(dt,j,p),c,D(dt,j,c)));
    /*    fprintf(stderr,"U (%d -> %d) p=%d c=%d ss=%d\n",
	    p,c, D(dt,j,p), D(dt,j,c),
	    FOREST_SS(frt,p,D(dt,j,p),c,D(dt,j,c))); */
  }

  /* calculate bde */

  bde = pcc * lgamma(aij) - pcc * vcc * lgamma(aijk);
  for(j = 0; j < pcc; ++j){
    double nij = 0;
    int k;
    for(k = 0; k < vcc; ++k){
      nij += FOREST_SS(frt,p,j,c,k);
      bde += lgamma(FOREST_SS(frt,p,j,c,k) + aijk);
    }
    bde -= lgamma(nij + aij);
  }

  /*   fprintf(stderr,"c=%d p=%d : %f\n",c,p,bde);  */

  /*
  for(j = 0; j < pcc; ++j){
    int k;
    for(k = 0; k < vcc; ++k){
      fprintf(stderr,"(%d -> %d) p=%d c=%d ss=%d\n",
	      p,c,j,k, FOREST_SS(frt,p,j,c,k));
    }
  }

  */
  return bde;
 
}

void
create_th(forest* frt) {
  int i;

  theta_free(frt);

  MEMALL(frt->th, frt->dim, double*)

  for(i=0; i<frt->dim; ++i){
    int p = frt->parent[i];
    int pcc = (p == -1) ? pcc = 1 : frt->valcount[p];
    MECALL(frt->th[i], pcc * frt->valcount[i], double)
  }
  
}


double forest_learn(forest* frt, data* dt){

  int first_free;
  int new_node;
  double bde;
  int i;

  /* clear data structures */
  for(i=0;i<frt->dim;++i){
    frt->parent[i] = 0;
    frt->connected[i] = 0;
  }

  get_orphan_lmlh(frt,dt);

  /* set distances to node 0 */
  new_node = 0;
  frt->connected[new_node] = 1;
  first_free = 1;
  for(i = 1; i < frt->dim; ++i) {
    frt->dst[i] = lmlh_term(frt, dt, i, new_node) - frt->orphan_lmlh[i];
  }
 
  /* then process other nodes */  
  for(i = 1; i < frt->dim; ++i) {

    /*    printf("Trying to find someone close to tree\n"); */

    /* find closest to tree that is not already connected */
    {
      int c;
      while(frt->connected[first_free]) ++first_free;
      new_node = first_free;
      for(c = first_free+1; c < frt->dim; ++c){
	/* printf("Trying %d - %f (%d)\n",c,frt->dst[c],frt->connected[c]); */
	if(frt->connected[c] || (frt->dst[c] < frt->dst[new_node])) continue; 
	new_node = c;
      }
    }

    /*    printf("Found %d\n",new_node); */

    frt->connected[new_node] = 1;

    /* update distances of the free nodes to the tree */
    {
      int c;
      for(c = first_free; c < frt->dim; ++c){
	double d_to_new;
	if(frt->connected[c]) continue;
	d_to_new = lmlh_term(frt, dt, c, new_node) - frt->orphan_lmlh[c];
	if(d_to_new > frt->dst[c]) {
	  frt->dst[c] = d_to_new;
	  frt->parent[c] = new_node;
	}
	/* printf("Updating distance of %d to tree, now %f\n",c,frt->dst[c]);*/
      }                               
    }
  }

  bde=0;
  for(i=0; i<frt->dim; ++i){
    /* printf("%d - %g %g\n",i,frt->dst[i],frt->orphan_lmlh[i]); */
    if((i==frt->parent[i]) || (frt->dst[i] <= 0.1) /* frt->orphan_lmlh[i]) */){
      frt->parent[i] = -1;
      bde += frt->orphan_lmlh[i];
    } else {
      bde += frt->dst[i] + frt->orphan_lmlh[i];
    }
  }

  create_th(frt);

  return bde;

}

double forest_learn_no_arcs(forest* frt, data* dt){
  double bde;
  int i;

  get_orphan_lmlh(frt,dt);
  
  bde=0;
  for(i=0; i<frt->dim; ++i){
    frt->parent[i] = -1;
    bde += frt->orphan_lmlh[i];
  }
  
  create_th(frt);

  return bde;
}

void forest_update_parameters(forest* frt, data* dt){
  int i;
  int j;

  for(i=0; i<frt->dim; ++i){
    int l;
    int p = frt->parent[i];
    int pcc = (p == -1) ? 1 : frt->valcount[p];
    int ss_size = pcc * frt->valcount[i];
    for(l=0; l<ss_size; ++l) {
      frt->th[i][l] = 0;
    }    
  }

  for(j=0; j<dt->N; ++j){
    for(i=0; i<frt->dim; ++i){
      int p = frt->parent[i];
      ++FOREST_TH(frt, p, ((p == -1) ? 0 : D(dt, j, p)), i, D(dt,j,i));
    }    
  }
  
  for(i=0; i<frt->dim; ++i){
    int pc;
    int p = frt->parent[i];
    int vci = frt->valcount[i];
    int pcc = (p == -1) ? 1 : frt->valcount[frt->parent[i]];

    double ess2 = frt->ess / pcc;
    double ess1 = ess2 / vci;

    for(pc=0; pc<pcc; ++pc){
      int cc;
      int N_p_is_pc = 0; 
      for(cc=0; cc < vci; ++cc){
	N_p_is_pc += FOREST_TH(frt,p,pc,i,cc);
      }
      for(cc=0; cc < vci; ++cc){
	FOREST_TH(frt,p,pc,i,cc) 
	  = (FOREST_TH(frt,p,pc,i,cc) + ess1) / (N_p_is_pc + ess2); 
      }
    }
  }
}

void forest_set_random_parameters(forest* frt){
  int i;
  for(i=0; i<frt->dim; ++i){
    int pc;
    int p = frt->parent[i];
    int vci = frt->valcount[i];
    int pcc = (p == -1) ? 1 : frt->valcount[frt->parent[i]];

    for(pc=0; pc<pcc; ++pc){
      int cc;
      double norm = 0; 
      for(cc=0; cc < vci; ++cc){
	norm += (FOREST_TH(frt,p,pc,i,cc) = (1.0*rand()/(RAND_MAX+1.0)));
      }
      for(cc=0; cc < vci; ++cc){
	FOREST_TH(frt,p,pc,i,cc) /= norm; 
      }
    }
  }
}

double forest_d1_lprob(forest* frt, data* dt, int j){
  int i;
  double lp = 0;
  for(i=0; i<frt->dim; ++i){
    int pv = (frt->parent[i] == -1) ? 0 : D(dt,j,frt->parent[i]);
    lp += log(FOREST_TH(frt,0,pv,i,D(dt,j,i)));
  }
  return lp;
}

void forest_write_structure(forest* frt, FILE* fp){
  int i;
  for(i=0; i<frt->dim; ++i){
    int p = frt->parent[i];
    if(p != -1) fprintf(fp,"%d -> ", p);
    fprintf(fp,"%d\n",i);
  }
}

void forest_write_params(forest* frt, FILE* fp){
  int i;
  for(i=0; i<frt->dim; ++i){
    int pc;
    int p = frt->parent[i];
    int vci = frt->valcount[i];
    int pcc = (p == -1) ? 1 : frt->valcount[frt->parent[i]];
    
    if(p != -1) fprintf(fp,"%d -> ", p);
    fprintf(fp,"%d:\n",i);

    for(pc=0; pc<pcc; ++pc){
      int cc;
      fprintf(fp,"%d: ", pc);
      for(cc=0; cc < vci; ++cc){
	    fprintf(fp," %.1f",100 * FOREST_TH(frt,p,pc,i,cc));
      }
      fprintf(fp,"\n");
    }
  }
}

parent_matrix* forest_get_parent_matrix(forest* frt) {
  int i;
  parent_matrix* mx = parent_matrix_create(frt->dim);
  CLEAR_PARENT_MATRIX(mx)
  for(i=0; i<frt->dim; ++i){
    int p = frt->parent[i];
    if(p != -1) SET_PARENT(i, p, mx);
  }
  return mx;
}

void forest_write(forest* frt, FILE* fp){
  int i;
  parent_matrix* mx;
  int* rootmap;
  int* is_root;
  int rootcount;
  int root_ord;

  /* Setup parent matrix */
  mx = forest_get_parent_matrix(frt);

  /* Find root for each node */
  rootcount = 0;
  MEMALL(rootmap, frt->dim, int)
  MECALL(is_root, frt->dim, int)
  for(i=0; i<frt->dim; ++i){
    int p;
    for(p=i; 1; p=frt->parent[p]){
      if(frt->parent[p] == -1){
	rootcount +=(!is_root[p]);
	is_root[p] = 1;
	rootmap[i] = p;
	break;
      }
    }
  }

  /* Print roots and ord_them */
  fprintf(fp, "%d %d\n", frt->dim, rootcount);
  root_ord = 0;
  for(i=0; i<frt->dim; ++i){
    if(is_root[i]) {
      fprintf(fp, "1 %d\n", i);
      is_root[i] = root_ord++;
    }
  }
  fprintf(fp, "\n");

  /* Print nodes */
  for(i=0; i<frt->dim; ++i){
    int j;
    int pc;
    int p = frt->parent[i];
    int vci = frt->valcount[i];
    int pcc = (p == -1) ? 1 : frt->valcount[frt->parent[i]];
    int childcount;

    fprintf(fp, "%d %d\n", vci, is_root[rootmap[i]]);
    
    if(p == -1) {
      fprintf(fp, "0\n");
    } else {
      fprintf(fp, "1 %d\n", p);
    }

    childcount = 0;
    for(j=0; j<frt->dim; ++j){
      childcount += IS_PARENT(j, i, mx);
    } 
    fprintf(fp, "%d", childcount);
    for(j=0; j<frt->dim; ++j){
      if(IS_PARENT(j, i, mx)) fprintf(fp, " %d", j);
    } 
    fprintf(fp, "\n");

    fprintf(fp, "%d\n", pcc);
    for(pc=0; pc<pcc; ++pc){
      int cc;
      for(cc=0; cc < vci; ++cc){
	fprintf(fp, " %g", FOREST_TH(frt,p,pc,i,cc));
      }
      fprintf(fp,"\n");
    }
    fprintf(fp,"\n");
  }

  parent_matrix_free(mx);

  free(rootmap);
  free(is_root);
}



#include <math.h>
#include <string.h>
#include "err.h"
#include "forest.h"
#include "parent_matrix.h"
#include "node.h"
#include "bane.h"

#define GET_PCI(PCI,BN,NODI,DT,J) \
{\
      Word_t p;\
      int base = 1;\
      (PCI) = 0;\
\
      /* calculate pci */\
      LOOP_OVER_PARENTS((NODI), p){\
        if(MISSING((DT),(J),p)) {\
          (PCI) = -1;\
           break;\
        }\
	(PCI) += base * D((DT),(J),p);\
	base *= (BN)->nodes[p].valcount;\
      } END_PARENT_LOOP\
\
}

#define CLEAR_SS(BN,NODI) \
{\
    memset((NODI)->ss,  0, (NODI)->pcc * (NODI)->valcount * sizeof(int));\
    memset((NODI)->ssp, 0, (NODI)->pcc * sizeof(int));\
}

#define ADD_SS(BN,NODI,DT,J) \
{\
      if(!MISSING((DT),(J),(NODI)->id)) {\
         int pci = 0;\
         GET_PCI(pci,BN,NODI,DT,J); \
         if(pci != -1) {\
           ++ SS((NODI),pci,D((DT),(J),(NODI)->id));\
           ++ (NODI)->ssp[pci];\
           ++ (NODI)->N;\
         }\
      }\
}

#define DEL_SS(BN,NODI,DT,J) \
{\
      if(!MISSING((DT),(J),(NODI)->id)) {\
         int pci = 0;\
         GET_PCI(pci,BN,NODI,DT,J); \
         if(pci != -1) {\
           -- SS((NODI),pci,D((DT),(J),(NODI)->id));\
           -- (NODI)->ssp[pci];\
           -- (NODI)->N;\
         }\
      }\
}

#define IS_ANCESTOR_OF(NODE,AID) ((NODE)->path_to_me_count[AID]>0)

void bane_clear_ss(bane* bn){
  int i;
  for(i=0; i<bn->nodecount; ++i){
    node* nodi = bn->nodes + i;
    CLEAR_SS(bn,nodi);
  }
}

void bane_add_ss(bane* bn, data* dt, int j){
  int i;
  for(i=0; i<bn->nodecount; ++i){
    node* nodi = bn->nodes + i;
    ADD_SS(bn,nodi,dt,j);
  }
}


void bane_del_ss(bane* bn, data* dt, int j){
  int i;
  for(i=0; i<bn->nodecount; ++i){
    node* nodi = bn->nodes + i;
    DEL_SS(bn,nodi,dt,j);
  }
}

double bane_logprob(bane* bn, data* dt, int j, double ess, 
		double* terms){
  int i;
  double lp = 0;
  for(i=0; i<bn->nodecount; ++i){
    node* nodi = bn->nodes + i;
    double lt;
    int pci;
    GET_PCI(pci,bn,nodi,dt,j);
    if((pci == -1) || MISSING(dt,j,i)){
      fprintf(stderr,"logprob asked with missing values at %d %d\n",j, i);
      exit(-1);
    }
    lp += lt = log(SS(nodi,pci,D(dt,j,i)) + ess) 
             - log(nodi->ssp[pci]         + nodi->valcount * ess); 
    if (terms) terms[i] = lt;
  }  
  return lp;
}

/********** NEW SCORING USING TEMPORARY HASHTABLE FOR SS ********************/

#define MISSING_PCFG  ((unsigned int) -1)

Pvoid_t gather_sparse_ss_for_v(bane* bn, data* dt, int v){

  node* nodv = bn->nodes + v;
  int vc_v = nodv->valcount; 

  Pvoid_t ss = NULL;

  /* Collect suff stats for node v */
  
  { /* GET SUFFSTATS */

    unsigned int *pcfgs  = (unsigned int*) calloc(dt->N,sizeof(unsigned int));
    unsigned int* end_pcfg = pcfgs + dt->N;

    { /* CALCULATE PCFGS (i.e. parent configurations) */
      
      Word_t p;
      LOOP_OVER_PARENTS(nodv, p){

	unsigned char* dp_j;
	unsigned int* pcfg_j;
	unsigned int vc_p = bn->nodes[p].valcount;

	for(dp_j = D_I(dt, p), pcfg_j = pcfgs; 
	    pcfg_j < end_pcfg; 
	    ++pcfg_j, ++dp_j){
	  if (*dp_j == MISSINGVALUE || *pcfg_j == MISSING_PCFG)
	    *pcfg_j = MISSING_PCFG;
	  else {
	    *pcfg_j *= vc_p;
	    *pcfg_j += *dp_j;
	  }
	}
      } END_PARENT_LOOP;
    }

    {
	/* ADD SUFFSTATS */

	unsigned char* dp_j;
	unsigned int* pcfg_j;
	Word_t* cnt_ptr;

	for(dp_j = D_I(dt, v), pcfg_j = pcfgs; 
	    pcfg_j < end_pcfg; 
	    ++pcfg_j, ++dp_j){
	  if (*dp_j == MISSINGVALUE || *pcfg_j == MISSING_PCFG) continue;
	  
	  JLG(cnt_ptr, ss, *pcfg_j);
	  if(!cnt_ptr){
	    JLI(cnt_ptr, ss, *pcfg_j);
	    /* printf("Adding table for %d pcfg %d\n", v, *pcfg_j); */
	    *cnt_ptr = (Word_t) calloc(vc_v, sizeof(unsigned int));	
	  }
	  ++ ((*((unsigned int **)cnt_ptr))[*dp_j]);
	}    
    }

    free(pcfgs);
  }
  
  return ss;
}

void free_sparse_ss(Pvoid_t ss){
  Word_t* cnt_ptr;
  int retcode;
  Word_t pcfg = 0;

  JLF(cnt_ptr, ss, pcfg);
  while(cnt_ptr != NULL) {
    free(*(unsigned int **)cnt_ptr);
    JLN(cnt_ptr ,ss, pcfg);
  }
  JLFA(retcode, ss);
}

double BDe_score_for_v_from_ss(bane* bn, Pvoid_t ss, int v, double ess){

  double res; 

  node* nodv = bn->nodes + v;
  int vc_v  = nodv->valcount; 
  int pcc_v = nodv->pcc; 

  Word_t pcfg = 0;
  Word_t* cnt_ptr;

  /* Calculate score freeing suff stats */

  int cc_v  = pcc_v * vc_v;

  double ess_per_pcc_v = ess / pcc_v;
  double ess_per_cc_v  = ess / cc_v;

  unsigned int nof_freqs;    
  JLC(nof_freqs, ss, 0, -1);

  res  = nof_freqs * lgamma(ess_per_pcc_v);
  res -= nof_freqs * vc_v * lgamma(ess_per_cc_v);
        
  pcfg = 0;
  JLF(cnt_ptr, ss, pcfg);
  while(cnt_ptr != NULL) {
    int l;
    unsigned int pcfreq = 0;
    for(l=0; l<vc_v; ++l){
      unsigned int lfreq = (*(unsigned int**) cnt_ptr)[l];
      pcfreq += lfreq;
      res += lgamma(ess_per_cc_v + lfreq);
    }
    res -= lgamma(ess_per_pcc_v + pcfreq);
    JLN(cnt_ptr ,ss, pcfg);
  }
  
  return res;
}

double BDe_score_for_v(bane* bn, data* dt, int v, double ess){
  Pvoid_t ss = gather_sparse_ss_for_v(bn, dt, v);
  double res = BDe_score_for_v_from_ss(bn, ss, v, ess);
  free_sparse_ss(ss);
  return res;
}

double BDe_score(bane* bn, data* dt, double ess, double* scoreboard){
  double res = 0;
  int i;
  for(i=0; i<bn->nodecount; ++i){
    double score_v =  BDe_score_for_v(bn, dt, i, ess);
    if(scoreboard) scoreboard[i] = score_v;
    res += score_v;
  }
  return res;
}

/********************* END OF NEW SCORING ***********************************/

void bane_gather_ss_for_i(bane* bn, data* dt, int i){
    int j;
    node* nodi = bn->nodes + i;

    if(NULL != nodi->ss){
      free(nodi->ss);
      free(nodi->ssp);
    }
    MECALL(nodi->ss, nodi->pcc * nodi->valcount, int);
    MECALL(nodi->ssp, nodi->pcc , int);
    
    nodi->N = 0;
    for(j=0; j<dt->N; ++j){
      ADD_SS(bn,nodi,dt,j);
    }
}

void bane_gather_full_ss(bane* bn, data* dt){
  int i;
  for(i=0; i<bn->nodecount; ++i){
    bane_gather_ss_for_i(bn, dt, i);
  }
}

void bane_gather_full_ss_in_order(bane* bn, data* dt){
  int i;
  for(i=0; i<bn->nodecount; ++i){
    int j;
    node* nodi = bn->nodes + i;

    if(NULL != nodi->ss){
      free(nodi->ss);
      free(nodi->ssp);
    }
    MECALL(nodi->ss, nodi->pcc * nodi->valcount, int);
    MECALL(nodi->ssp, nodi->pcc , int);
    
    nodi->N = 0;
    for(j=0; j<dt->N; ++j){
      int p;
      int base = 1;
      int pci = 0;

      if(MISSING(dt,j,i)) continue;

      /* calculate pci */
      for(p=0; p<bn->nodecount; ++p){
	if(!IS_PARENT(i,p,bn->pmx)) continue;
	if(MISSING(dt,j,p)) {
	  pci = -1;
	  break;
	}
	pci += base * D(dt,j,p);
	base *= bn->nodes[p].valcount;
      }
      
      if (pci == -1) continue;

      ++ SS(nodi,pci,D(dt,j,i));
      ++ nodi->ssp[pci];
      ++ nodi->N;
    }
  }
}

void bane_gather_full_ss_in_order_from_file(bane* bn, char* filename){
  FILE* fp;  
  int i;
  unsigned char* dr;
  int tmp;

  /* create space for SS */
  for(i=0; i<bn->nodecount; ++i){
    node* nodi = bn->nodes + i;

    if(NULL != nodi->ss){
      free(nodi->ss);
      free(nodi->ssp);
    }
    MECALL(nodi->ss, nodi->pcc * nodi->valcount, int);
    MECALL(nodi->ssp, nodi->pcc , int);    
    
    nodi->N = 0;
  }

  MECALL(dr, bn->nodecount, char);
  OPENFILE_OR_DIE(fp,filename,"r");

  while(EOF != fscanf(fp, "%d", &tmp)){
    /* read the line */
    dr[0]=(char) tmp;
    for(i=1; i<bn->nodecount; ++i){
      if(1 != fscanf(fp, "%d", &tmp)){
	fprintf(stderr, "error while reading %s\n", filename);
	exit(ERROR_IN_DATAREAD);	
      }
      dr[i]=(char) tmp;
    }

    /* update SS for each variable */

    for(i=0; i<bn->nodecount; ++i){
      int p;
      int base = 1;
      int pci = 0;
      node* nodi = bn->nodes + i;

      if(dr[i] == MISSINGVALUE) continue;
    
      /* calculate pci */
      for(p=0; p<bn->nodecount; ++p){
	if(!IS_PARENT(i,p,bn->pmx)) continue;
	if(dr[p] == MISSINGVALUE) {
	  pci = -1;
	  break;
	}
	pci += base * dr[p];
	base *= bn->nodes[p].valcount;
      }
    
      if (pci == -1) continue;
    
      ++ SS(nodi,pci,dr[i]);
      ++ nodi->ssp[pci];
      ++ nodi->N;
    }
  }
	
  CLOSEFILE_OR_DIE(fp,filename);
  free(dr);
}

bane* bane_create(int nodecount){
  bane* bn;

  MEMALL(bn,1,bane)

  bn->nodecount = nodecount;

  bn->pmx = parent_matrix_create(nodecount);
  CLEAR_PARENT_MATRIX(bn->pmx)

  MECALL(bn->nodes, nodecount, node)
  {
    int i;
    for(i=0; i<bn->nodecount; ++i){
      node* nodi = bn->nodes + i;
      MEMALL(nodi->path_to_me_count, nodecount, int);
      node_init(nodi, i, nodecount);
    }
  }

  return bn;
}

bane* bane_create_from_format(format* fmt){
  bane* bn = bane_create(fmt->dim);
  int i;
  for(i=0; i<bn->nodecount; ++i){
    node* nodi = bn->nodes + i;
    nodi->valcount = fmt->valcount[i];
  }
  return bn;
}

void bane_assign(bane* dst, bane* src) {
  int i;
  dst->nodecount = src->nodecount;
  parent_matrix_assign(dst->pmx, src->pmx);
  for(i=0; i<dst->nodecount; ++i){
    node_assign(dst->nodes+i, src->nodes + i, dst->nodecount);
  }
}

bane* bane_copy(bane* src) {
  bane* dst = bane_create(src->nodecount);
  bane_assign(dst, src);
  return dst;
}

void bane_copy_ss(bane* dst, bane* src){
  int i;
  for(i=0; i<dst->nodecount; ++i){
    node_copy_ss(dst->nodes + i, src->nodes + i);
  }
}


bane* bane_create_forest(format* fmt, double ess, data* dt){

  bane* bn = bane_create_from_format(fmt);
  forest* frt = forest_create(fmt, ess);
  double score = forest_learn(frt, dt);
  arc* ar;
  MEMALL(ar,1,arc);
  {
    int i;
    for(i=0; i<bn->nodecount; ++i){
      if(frt->parent[i] > -1){
	ar->from = frt->parent[i];
	ar->to = i;
	bane_add_arc(bn, ar, 1);
      }
    }
  }
  free(ar);
  forest_free(frt);
  score = 0;
  return bn;
}

bane* bane_create_from_pmx(format* fmt, parent_matrix* pmx){
  bane* bn = bane_create_from_format(fmt);
  arc* ar;
  MEMALL(ar,1,arc);

  {
    int c;
    for(c=0; c<bn->nodecount; ++c){
      int p;
      for(p=0;p<bn->nodecount; ++p){
	if(!IS_PARENT(c,p,pmx)) continue;
	if((c==p) || (IS_ANCESTOR_OF(bn->nodes+p, c))) {
	  bane_free(bn);
	  return NULL;
	}
	ar->from = p;
	ar->to = c;
	bane_add_arc(bn, ar, 1);
      }
    }
  }
  free(ar);
  return bn;
}

void propagate_path_add(bane* bn, node* me, int* addend, int new_parent_id){
  int i;

  /* first update me */
  ++ me->path_to_me_count[new_parent_id];

  me->ancestorcount = 0;
  for(i=0; i<bn->nodecount; ++i){
    me->path_to_me_count[i] += addend[i];
    me->ancestorcount += IS_ANCESTOR_OF(me,i);
  }

  /* and then propagate to children */

  {
    Word_t c;
    LOOP_OVER_CHILDREN(me, c){
      propagate_path_add(bn, bn->nodes+c, addend, new_parent_id);
    } END_CHILD_LOOP;
  }

}
  
void bane_fix_add_arc(bane* bn, arc* ar){
  node* from = bn->nodes + ar->from;
  node* to   = bn->nodes + ar->to;
  propagate_path_add(bn, to, from->path_to_me_count, from->id);
}

void bane_add_arc(bane* bn, arc* ar, int commit){

  node* from = bn->nodes + ar->from;
  node* to   = bn->nodes + ar->to;
  SET_PARENT(to->id, from->id, bn->pmx);
  ADD_CHILD(from,ar->to);
  ADD_PARENT(to,ar->from,from->valcount);

  if (commit) bane_fix_add_arc(bn, ar);
  /*
  fprintf(stderr,"Added arc %d -> %d\n",from->id, to->id);
  node_write(bn,from,stderr);
  */
}

void propagate_path_del(bane* bn, node* me, int* subend, int old_parent_id){
  int i;

  /* first update me */
  -- me->path_to_me_count[old_parent_id];

  me->ancestorcount = 0;
  for(i=0; i<bn->nodecount; ++i){
    me->path_to_me_count[i] -= subend[i];
    me->ancestorcount += IS_ANCESTOR_OF(me,i);
  }

  /* and then propagate to children */

  {
    Word_t c;
    
    LOOP_OVER_CHILDREN(me, c){
      propagate_path_del(bn, bn->nodes + c, subend, old_parent_id);
    } END_CHILD_LOOP;
  }

}
  
void bane_fix_del_arc(bane* bn, arc* ar){
  node* from = bn->nodes + ar->from;
  node* to   = bn->nodes + ar->to;
  propagate_path_del(bn, to, from->path_to_me_count, from->id);

}

void bane_del_arc(bane* bn, arc* ar, int commit){
  node* from = bn->nodes + ar->from;
  node* to   = bn->nodes + ar->to;
  UNSET_PARENT(to->id, from->id, bn->pmx);
  DEL_CHILD(from,ar->to);
  DEL_PARENT(to,ar->from,from->valcount);

  if (commit) bane_fix_del_arc(bn, ar);
  /*
  fprintf(stderr,"Deleted arc %d -> %d\n",from->id, to->id);
  node_write(bn,from,stderr);
  */
}

void bane_fix_rev_arc(bane* bn, arc* ar){
  bane_rev_arc(bn, ar, 0); /* reverse it back quick */
  bane_rev_arc(bn, ar, 1); /* reverse it again now commiting */
}

void bane_rev_arc(bane* bn, arc* ar, int commit){
  int tmp;
  bane_del_arc(bn,ar,commit);

  tmp = ar->from;      /* save from   */
  ar->from = ar->to;   /* update from */  /* swap ar */
  ar->to = tmp;        /* update to   */

  bane_add_arc(bn,ar,commit);

  /*
  fprintf(stderr,"Reversed arc %d -> %d\n",from->id, to->id);
  node_write(bn,from,stderr);
  */
}

extern double lgamma(double);


/*
double bane_get_score_for_i_if_del_x(node* nodi, double ess, node* nodx){
  # save old ss and ssp and pcc
  # update pcc
  # create and calculate new ss and ssp from old ones
  # scori = bane_get_score_for_i
  # restore old ss and ssp and pcc
  # free new ss and ssp
  # return scori
}
*/

double bane_get_score_for_i(node* nodi, double ess){
  int pci;
  double esp1 = ess / nodi->pcc;
  double esp2 = ess / (nodi->pcc * nodi->valcount);
  double scori = 
    nodi->pcc * lgamma(esp1) - nodi->pcc * nodi->valcount * lgamma(esp2);
  
  /* This can be optimized */
  for(pci = 0; pci < nodi->pcc; ++pci){
    int v;
    for(v=0; v<nodi->valcount; ++v){
      scori += lgamma(esp2 + SS(nodi,pci,v));
    }
    scori -= lgamma(esp1 + nodi->ssp[pci]);
  }
  return scori;
}

double bane_get_score(bane* bn, double ess, double* scoreboard){
  int i;
  double score = 0;
  for(i=0; i<bn->nodecount; ++i){
    double iscore = bane_get_score_for_i(bn->nodes + i, ess);
    score += iscore;
    if(NULL != scoreboard) scoreboard[i] = iscore; 
  }
  return score;
}

void
bane_free(bane* bn){
  int i;
  int retcode;
  for(i=0; i<bn->nodecount; ++i){
    node* nodi = bn->nodes + i;
    J1FA(retcode, nodi->parents);
    J1FA(retcode, nodi->children);
    free(nodi->ss);
    free(nodi->ssp);
    free(nodi->path_to_me_count);
  }
  free(bn->nodes);
  parent_matrix_free(bn->pmx);
  free(bn);
}


void
bane_write_with_ss(bane* bn, char* filename){
  int i;
  FILE* fp;

  if(NULL == (fp = fopen(filename, "w"))){
    fprintf(stderr, "can't open %s\n", filename);
    exit(CANNOT_OPEN_OUTPUTFILE);
  }


  fprintf(fp,"%d\n", bn->nodecount);

  for(i=0; i<bn->nodecount; ++i){
    int p;
    node* nodi = bn->nodes + i;

    fprintf(fp,"%d %d %d %d",
            nodi->valcount, nodi->childcount, nodi->parentcount, nodi->pcc);

    for(p=0; p<bn->nodecount; ++p){
      if(!IS_PARENT(i,p,bn->pmx)) continue;
      fprintf(fp, " %d", p);
    }

    fprintf(fp, " %d", nodi->N);
    for(p=0; p<nodi->pcc * nodi->valcount; ++p){
      fprintf(fp, " %d", nodi->ss[p]);
    }
    fprintf(fp,"\n");
  }

  if(fclose(fp)){
    fprintf(stderr, "can't close %s\n", filename);
    exit(CANNOT_CLOSE_OUTPUTFILE);
  }
  
}

bane* 
bane_read_with_ss(char* filename){
  int i;
  FILE* fp;
  int nodecount;
  bane* bn;

  parent_matrix* pmxtmp;

  if(NULL == (fp = fopen(filename, "r"))){
    fprintf(stderr, "can't open %s\n", filename);
    exit(CANNOT_OPEN_INPUTFILE);
  }

  
  fscanf(fp,"%d\n", &nodecount);
  bn = bane_create(nodecount);
  pmxtmp = parent_matrix_create(bn->nodecount);

  for(i=0; i<bn->nodecount; ++i){
    int p;
    node* nodi = bn->nodes + i;
    
    fscanf(fp,"%d %d %d %d",&(nodi->valcount), 
	    &(nodi->childcount), &(nodi->parentcount), &(nodi->pcc));
    
    for(p=0; p<nodi->parentcount; ++p){
      int pid;
      fscanf(fp, " %d", &pid);
      SET_PARENT(i, pid, pmxtmp);
    }

    MEMALL(nodi->ss, nodi->pcc * nodi->valcount, int);
    MECALL(nodi->ssp, nodi->pcc, int);

    fscanf(fp, " %d", &(nodi->N));
    for(p=0; p<nodi->pcc * nodi->valcount; ++p){
      fscanf(fp, " %d", &(nodi->ss[p]));
    }

    for(p=0; p<nodi->pcc; ++p){
      int v;
      for(v=0; v<nodi->valcount; ++v){
	nodi->ssp[p] += SS(nodi,p,v);
      }
    }

    nodi->childcount = 0;
    nodi->parentcount = 0;
    nodi->pcc = 1;

  }
 
  if(fclose(fp)){
    fprintf(stderr, "can't close %s\n", filename);
    exit(CANNOT_CLOSE_INPUTFILE);
  }
  
  {
    int c;
    arc* ar;
    MEMALL(ar,1,arc);
    for(c=0; c<bn->nodecount; ++c){
      int p;
      for(p=bn->nodecount-1; p>=0; --p){
	if(IS_PARENT(c,p,pmxtmp)){
	  ar->from = p;
	  ar->to = c;
	  bane_add_arc(bn, ar, 1);
	}
      }
    }
    free(ar);
  }

  parent_matrix_free(pmxtmp);

  return bn;
}

void bane_write_structure(bane* bn, FILE* fp) {
  int i;

  fprintf(fp,"%d\n", bn->nodecount);
  
  for(i=0; i<bn->nodecount; ++i){
    int p;
    node* nodi = bn->nodes + i;

    fprintf(fp,"%d %d", nodi->childcount, nodi->parentcount);

    for(p=0; p<bn->nodecount; ++p){
      if(IS_PARENT(i,p,bn->pmx))
	fprintf(fp, " %d", p);
    }
    fprintf(fp,"\n");
  }
}


void bane_read_structure_n_read(bane* bn, FILE* fp, int n) {
  arc* ar;
  MEMALL(ar,1,arc);

  for(ar->to = 0; ar->to < bn->nodecount; ++ ar->to){
    int p;
    int parc, chlc;

    fscanf(fp,"%d %d", &chlc, &parc);
    
    for(p=0; p<parc; ++p){
      fscanf(fp, "%d", &(ar->from));
      bane_add_arc(bn, ar, 1);
    }
  }
  free(ar);
}

void bane_read_structure(bane* bn, FILE* fp) {
  int n;
  fscanf(fp,"%d", &n);  /* read nodecount */
  bane_read_structure_n_read(bn, fp, n);
}

bane* bane_cread_structure(FILE* fp) {
  bane* bn;
  int nodecount;
  
  fscanf(fp,"%d",&nodecount);
  bn = bane_create(nodecount);
  bane_read_structure_n_read(bn,fp,nodecount);

  return bn;

}

#if 0
    for(p=0; p<nodi->pcc; ++p){
      int v;
      int nsum=0;
      for(v=0; v<nodi->valcount; ++v){
	th[v] = SS(nodi,p,v) + esp2;
	nsum += SS(nodi,p,v);
      }
      for(v=0; v<nodi->valcount; ++v){
	fprintf(fp, " %f", th[v] / (esp1+nsum));      
      }
    }
#endif


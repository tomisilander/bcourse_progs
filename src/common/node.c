#include <stdio.h>
#include <string.h>
#include "err.h"
#include "node.h"

void 
node_assign(node* dst, node* src, int nodecount){
  Word_t c, p;
  int retcode;

  dst->id = src->id;
  dst->valcount = src->valcount;

  J1FA(retcode, dst->parents);
  dst->parentcount = 0;
  dst->parents     = (Pvoid_t) NULL;
  LOOP_OVER_PARENTS(src,p){
    ADD_PARENT(dst,p,0);
  } END_PARENT_LOOP;

  J1FA(retcode, dst->children);
  dst->childcount = 0;
  dst->children   = (Pvoid_t) NULL;
  LOOP_OVER_CHILDREN(src, c){
    ADD_CHILD(dst,c);
  } END_CHILD_LOOP;


  dst->pcc = src->pcc;
  dst->N = src->N;
  dst->ancestorcount = src->ancestorcount;
  memcpy(dst->path_to_me_count, src->path_to_me_count, nodecount*sizeof(int));
}

void node_copy_ss(node* dst, node* src){
   MEMALL(dst->ss,  src->pcc * src->valcount, int);
   MEMALL(dst->ssp, src->pcc , int); 
   node_assign_ss(dst,src);
}

void node_assign_ss(node* dst, node* src){
  memcpy(dst->ss,  src->ss,  dst->pcc * dst->valcount * sizeof(int));
  memcpy(dst->ssp, src->ssp, dst->pcc * sizeof(int));
}


void 
node_init(node* nodi, int i, int nodecount){
  int j;
  nodi->id = i;
  nodi->valcount = 1;
  nodi->parentcount = 0;
  nodi->parents = (Pvoid_t) NULL;
  nodi->childcount = 0;
  nodi->children = (Pvoid_t) NULL;
  nodi->pcc = 1;
  nodi->N = 0;
  nodi->ss = NULL;
  nodi->ssp = NULL;
  nodi->ancestorcount = 0;
  for(j=0; j<nodecount; ++j) nodi->path_to_me_count[j] = 0;
}


void node_write(node* nd, int nodecount, FILE* fp){
  int i;
  fprintf(fp, "id=%d vals=%d parents=%d children=%d ancestors=%d\n",
	  nd->id, nd->valcount, nd->parentcount, nd->childcount, 
	  nd->ancestorcount);
  
  /*
  fprintf(fp, "parents: %d ",nd->first_parent);
  for(i=0; i<nodecount;++i){
    fprintf(fp, " %d",nd->parent[i]);
  } 
  fprintf(fp, "\n");
  */

  /*
  fprintf(fp, "children: %d ",nd->first_child);
  for(i=0; i<nodecount;++i){
    fprintf(fp, " %d",nd->child[i]);
  } 
  fprintf(fp, "\n");
  */

  fprintf(fp, "paths from ancestors: ");
  for(i=0; i<nodecount;++i){
    fprintf(fp, " %d",nd->path_to_me_count[i]);
  } 
  fprintf(fp, "\n");
  fprintf(fp, "pcc=%d\n",nd->pcc);
  
}




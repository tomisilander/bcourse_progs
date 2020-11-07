#include "err.h"
#include "bitset.h"
#include "node.h"
#include "bane.h"
#include "nb_selector.h"

int* nb_selector_cread(char* fn){
  int* s;
  int i;
  bane* bn;
  FILE* fp;

  OPENFILE_OR_DIE(fp,fn,"r");
  bn = bane_cread_structure(fp);
  CLOSEFILE_OR_DIE(fp,fn);
  
  BITSET_CREATE(s,bn->nodecount);

  for(i=0; i<bn->nodecount; ++i){
    if(bn->nodes[i].parentcount + bn->nodes[i].childcount > 0){
      BITSET_SET(s,i);
    }
  }

  return s;
}

void nb_selector_write(int* s, int n, int classix, char* fn){
  int i;
  bane* bn;
  arc ar;
  FILE* fp;

  bn = bane_create(n);  

  ar.from = classix;
  for(i=0; i<n; ++i){
    if((i!=classix) && (BITSET_ON(s,i))){
      ar.to = i;
      bane_add_arc(bn, &ar, 0);
    }
  }

  OPENFILE_OR_DIE(fp,fn,"w");
  bane_write_structure(bn,fp);
  CLOSEFILE_OR_DIE(fp,fn);

  bane_free(bn);
}

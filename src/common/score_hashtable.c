#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "err.h"
#include "score_hashtable.h"

extern int
score_hashtable_items_for_mem(int memsize, int key_size){
  memsize -= sizeof(score_hashtable); 
  memsize -= 256 * sizeof(ht_pos);
  memsize /= 256 * (key_size * sizeof(unsigned int) + sizeof(double));
  return (memsize < 1) ?  0 : memsize;
}

score_hashtable*
score_hashtable_create(int key_size, int items_per_pos){
  int p;

  score_hashtable* ht;
  MEMALL(ht, 1, score_hashtable)

  ht->key_size = key_size;
  ht->chars_per_key = key_size * sizeof(unsigned int) / sizeof(unsigned char);
  ht->items_per_pos = items_per_pos;
  ht->keycount = 0;

  MEMALL(ht->tbl,256,ht_pos);
  for(p=0; p<256; ++p){
    MECALL(ht->tbl[p].keys, ht->items_per_pos * ht->key_size, unsigned int);
    MECALL(ht->tbl[p].items, ht->items_per_pos * ht->key_size, double);
    ht->tbl[p].head = 0;
    ht->tbl[p].count = 0;
  }

  return ht;
}

void
score_hashtable_free(score_hashtable* ht){
  int p;
  for(p=0; p<256; ++p){
    free(ht->tbl[p].keys);
    free(ht->tbl[p].items);
  }
  free(ht->tbl);
  free(ht);
}

int keyqual(unsigned int* k1, unsigned int* k2, int key_size){
  int i;
  for(i=0; i<key_size; ++i){
    if(k1[i] != k2[i]) {
      return 0;
    }
  }
  return 1;
}

unsigned char addr(unsigned int* k, int chars_per_key){
  unsigned char* ckp = (unsigned char*) k;
  unsigned char a = 0;
  int i;
  for(i=0; i < chars_per_key; ++i){
    a ^= ckp[i];
  }
  return a;
}

void 
score_hashtable_put(score_hashtable* ht, unsigned int* key , double item){
  ht_pos* htp = ht->tbl + addr(key,ht->chars_per_key);
  unsigned int* htpkey = htp->keys + htp->head * ht->key_size;

  int i;

#ifdef DEBUG
  {
    double testdbl;
    if(score_hashtable_get(ht, key, &testdbl)){
      fprintf(stderr,"key already in table -> %f\n",testdbl);
    }
  }
#endif  

  for(i=0; i<ht->key_size; ++i){
    htpkey[i] = key[i];
  }
  htp->items[htp->head] = item;
  htp->head = (htp->head+1) % ht->items_per_pos;
  if(htp->count < ht->items_per_pos) {
    ++htp->count;
    ++ht->keycount;
  }
}

int
score_hashtable_get(score_hashtable* ht, unsigned int* key, double* item){
  ht_pos* htp = ht->tbl + addr(key,ht->chars_per_key);
  int i;

  for(i=1; i<=htp->count; ++i){
    int sp = ((htp->head-i) >= 0) ? (htp->head-i) : htp->count-1;
    unsigned int* p = htp->keys + sp  * ht->key_size;

    if(keyqual(key,p,ht->key_size)){
      *item = htp->items[sp];
      return 1;
    }
  }  
  return 0;
}

void
score_hashtable_write(FILE* fp, score_hashtable* ht){
  int p;
  fprintf(fp,"ks=%d ch_p_k=%d i_p_p=%d\n",
	  ht->key_size, ht->chars_per_key, ht->items_per_pos);

  for(p=0; p<256; ++p){
    ht_pos* htp = ht->tbl + p;
    int i;
    fprintf(fp,"pos=%d count=%d head=%d\n",p, htp->count, htp->head);
    for(i=0; i<htp->count; ++i){
      int k;
      fprintf(fp,"%f\t",htp->items[i]);
      for(k=0; k<ht->key_size;++k){
	fprintf(fp,"%u ",*(htp->keys+i*ht->key_size+k));
      }
      fprintf(fp,"\n");
    }
  }  
}


typedef struct ht_entry ht_entry;

struct ht_entry {
  unsigned int* key;
  double item;
};


unsigned int* 
score_hashtable_get_nth_key(score_hashtable* ht, int n){ 
/* assumes n <= ht->keycount */
  int p;
  int hte_count = 0;
  ht_entry* best_htes;
  ht_entry* worst_hte = NULL;

  unsigned int* res = NULL;

  MECALL(best_htes,n,ht_entry);

  for(p=0; p<256; ++p){
    int i;
    ht_pos* htp = ht->tbl + p;

    for(i=1; i<=htp->count; ++i){
      int sp = ((htp->head-i) >= 0) ? (htp->head-i) : htp->count-1;

      if(hte_count < n) { 
	best_htes[hte_count].key = htp->keys + sp  * ht->key_size;
	best_htes[hte_count].item = htp->items[sp];
	++ hte_count;
	if(hte_count == n){
	  int h; /* find worst */
	  for(h=0; h<n; ++h){
	    if((h==0) || (best_htes[h].item < worst_hte->item)){
	      worst_hte = best_htes + h;
	    }
	  }
	}
      } else if (htp->items[sp] > worst_hte->item){
	int h;
	worst_hte->key = htp->keys + sp  * ht->key_size;
	worst_hte->item = htp->items[sp];	
	for(h=0; h<n; ++h){
	  if((best_htes[h].item < worst_hte->item)){
	    worst_hte = best_htes + h;
	  }
	}
      }
    }  
  }


  if(worst_hte != NULL) {
    MECALL(res, ht->key_size, unsigned int);
    memcpy(res, worst_hte->key,  ht->key_size * sizeof(unsigned int));
  }
  
  free(best_htes);

  return res;
}



#ifdef SCORE_HASHTABLE_TEST
#include "parent_matrix.h"

void main(int argc, char* argv[]){
  int m;
  parent_matrix* pmx;
  score_hashtable* ht;
  double item;
  int r;

  if(argc !=2){
    fprintf(stderr,"%s attcount\n",argv[0]);
    exit(-1);
  }

  m = atoi(argv[1]);
  pmx = parent_matrix_create(m);

  for(r=0; r<9*m; ++r){
    int c = lrand48()%m;
    int p = lrand48()%m;
    SET_PARENT(c,p,pmx);
  }

  ht = score_hashtable_create(pmx->one_dim_size,4);

  for(r=0; r<m; ++r){
    score_hashtable_put(ht,GET_PARENT_VECTOR(r,pmx),drand48());
  }

  score_hashtable_write(stdout,ht);

  for(r=0; r<m; ++r){
    if(0==score_hashtable_get(ht,GET_PARENT_VECTOR(r,pmx),&item)){
      printf("IIK - not found nothing\n");
      continue;
    }
    printf("%f\n",item);
  }
}

#endif

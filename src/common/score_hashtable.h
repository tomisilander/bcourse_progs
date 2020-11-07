#ifndef SCORE_HASHTABLE_H
#define SCORE_HASHTABLE_H

typedef struct ht_pos ht_pos;

struct ht_pos {
  unsigned int* keys;
  double* items;
  int count;
  int head;
};

typedef struct score_hashtable score_hashtable;

struct score_hashtable {
  int key_size;
  int chars_per_key;
  int items_per_pos;
  ht_pos* tbl;
  int keycount;
};

extern score_hashtable*
score_hashtable_create(int key_size, int items_per_pos);

extern void
score_hashtable_free(score_hashtable* ht);

extern void 
score_hashtable_put(score_hashtable* ht, unsigned int* key , double item);

extern int
score_hashtable_get(score_hashtable* ht, unsigned int* key, double* item);

extern unsigned int* 
score_hashtable_get_nth_key(score_hashtable* ht, int n);

extern int
score_hashtable_items_for_mem(int memsize, int key_size);

#endif


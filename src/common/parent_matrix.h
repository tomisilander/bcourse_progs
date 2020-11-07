#ifndef PARENT_MATRIX_H
#define PARENT_MATRIX_H

typedef struct parent_matrix parent_matrix;

struct parent_matrix {
  int m;
  int one_dim_size;
  int bits_in_int;
  unsigned int* mx;
};

extern parent_matrix*
parent_matrix_create(int m) ;

extern parent_matrix*
parent_matrix_create_wrap(int m, unsigned int* mx) ;

extern void
parent_matrix_assign(parent_matrix* dst, parent_matrix* src);

extern parent_matrix*
parent_matrix_copy(parent_matrix* src);

extern void
parent_matrix_free(parent_matrix* pmx);

extern void
parent_matrix_print(parent_matrix* pmx, FILE* fp);

#define SET_PARENT(C,P,MX) \
  (*((MX)->mx + (C) * (MX)->one_dim_size + (P)/(MX)->bits_in_int) \
   |= (1U << (P)%(MX)->bits_in_int))

#define UNSET_PARENT(C,P,MX) \
  (*((MX)->mx + (C) * (MX)->one_dim_size + (P)/(MX)->bits_in_int) \
   &= ~(1U << (P)%(MX)->bits_in_int))

#define IS_PARENT(C,P,MX) \
  ((*((MX)->mx + (C) * (MX)->one_dim_size + (P)/(MX)->bits_in_int) \
   & (1U << (P)%(MX)->bits_in_int)) != 0)

#define IS_CHILD(P,C,MX)  IS_PARENT((C),(P),(MX))

#define GET_PARENT_VECTOR(C,MX)  ( (MX)->mx + (C) * (MX)->one_dim_size )

#define CLEAR_PARENT_MATRIX(MX)  { \
  int i;\
  int tsize = (MX)->m  * (MX)->one_dim_size;\
  for(i=0; i<tsize; ++i) (MX)->mx[i] = 0U;\
}


#endif

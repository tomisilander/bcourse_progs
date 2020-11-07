#ifndef BITSET_H
#define BITSET_H

#include "err.h"

#define BITSET_INTCOUNT(N) (1 + ((N)-1)/sizeof(int))

#define BITSET_CREATE(S,N) MECALL(S, BITSET_INTCOUNT(N), int);

#define BITSET_CLEAR(S,N) \
{unsigned int i=0;for(; i < BITSET_INTCOUNT(N); ++i) *((S)+i) ^= *((S)+i);}

#define BITSET_ASSIGN(D,S,N) \
{unsigned int i=0;for(; i < BITSET_INTCOUNT(N); ++i) *((D)+i) = *((S)+i);}

#define BITSET_COPY(D,S,N) {BITSET_CREATE(D,N); BITSET_ASSIGN(D,S,N)}

#define BITSET_FREE(S) {free(S);}

#define BITSET_SET(S,I) \
  (*((S) + (I)/(sizeof(int)))  |= (1U << (I)%(sizeof(int))))

#define BITSET_UNSET(S,I) \
  (*((S) + (I)/(sizeof(int)))  &= ~(1U << (I)%(sizeof(int))))

#define BITSET_ON(S,I) \
  (*((S) + (I)/(sizeof(int)))  & (1U << (I)%(sizeof(int))))

#define BITSET_OFF(S,I)  (! BITSET_ON(S,I))


#endif

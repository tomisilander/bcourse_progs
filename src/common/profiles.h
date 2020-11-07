#ifndef PROFILES_H
#define PROFILES_H

typedef struct profiles profiles;

struct profiles {
  int K;
  int var_prof_len;
  int* valcs;
  int* var_ofs;
  /* int* val_ofs; */
  double* p;
};

#define PROFSKI(PFS,K,I) ((PFS)->p+(PFS)->var_ofs[I] + (K)*(PFS)->valcs[I])
#define PROFSKIL(PFS,K,I,L) (PROFSKI(PFS,K,I)[L])

extern profiles* 
profiles_create(int K, format* fmt);

extern void      
profiles_free(profiles* profs);

extern void      
profiles_read(char* filename, format* fmt, profiles* profs);

extern void      
profiles_write(char* filename, format* fmt, profiles* profs);

extern profiles* 
profiles_cread(char* filename, int K, format* fmt);

#define LOOP_KIL(K,FMT) \
{\
  int k,i,l;\
  for(k=0; k<(K); ++k){\
    for(i=0; i<(FMT)->dim; ++i){\
      for(l=0; l<(FMT)->valcount[i]; ++l)

#define END_KIL \
    }\
  }\
}


#endif

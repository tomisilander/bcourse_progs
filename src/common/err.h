#ifndef ERR_H
#define ERR_H

#include <stdio.h>
#include <stdlib.h>

enum error_codes {
  OUT_OF_MEMORY = 1,
  CANNOT_OPEN_DATAFILE,
  ERROR_IN_DATAREAD,
  CANNOT_CLOSE_DATAFILE,
  CANNOT_OPEN_FORMATFILE,
  ERROR_IN_FORMATREAD,
  CANNOT_CLOSE_FORMATFILE,
  CANNOT_OPEN_INPUTFILE,
  CANNOT_CLOSE_INPUTFILE,
  CANNOT_OPEN_OUTPUTFILE,
  CANNOT_CLOSE_OUTPUTFILE,
  CANNOT_OPEN_FILE,
  CANNOT_CLOSE_FILE
};


#define MEMALL(V,T,TYPE) \
{\
  if(NULL == ((V)=(TYPE *) malloc((T)*sizeof(TYPE)))) {\
    fprintf(stderr,"Could not allocate %d times sizeof(%s)\n",T,#TYPE);\
    exit(OUT_OF_MEMORY);\
  }\
}

#define MECALL(V,T,TYPE) \
{\
  if(NULL == ((V)=(TYPE *) calloc((T),sizeof(TYPE)))) {\
    fprintf(stderr,"Could not callocate %d times sizeof(%s)\n",T,#TYPE);\
    exit(OUT_OF_MEMORY);\
  }\
}

#define OPENFILE_OR_DIE(FP,FN,M) \
{\
  if(NULL == ((FP) = fopen((FN), (M)))){\
    fprintf(stderr, "can't open %s\n", (FN));\
    exit(CANNOT_OPEN_FILE);\
  }\
}

#define CLOSEFILE_OR_DIE(FP,FN) \
{\
  if(fclose(FP)){\
    fprintf(stderr, "can't close %s\n", (FN));\
    exit(CANNOT_CLOSE_FILE);\
  }\
}
 

#endif




#ifndef DATA_H
#define DATA_H

typedef struct data data;

struct data {
  unsigned char* d;
  int m;
  int N;
};

extern data* 
data_create(int datacount, int attcount);

extern data* 
data_cread(int datacount, int attcount, char* filename);

extern data*
data_copy(data* old_dt);

extern void
data_read(char* filename, data* dt);

extern void 
data_write(char* filename, data* dt);

extern void 
data_free(data* d);

#define D(DT,J,I) ((DT)->d[(I) * (DT)->N + (J)])
#define D_I(DT,I) ((DT)->d + (I) * (DT)->N)

#define MISSINGVALUE  ((unsigned char) 255)
#define MISSING(DT,J,I) (MISSINGVALUE == D(DT,J,I))

#endif


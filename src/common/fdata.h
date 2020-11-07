#ifndef FDATA_H
#define FDATA_H

typedef struct fdata fdata;

struct fdata {
  int om;
  float* d;
  int refcount;
  int m;
  int N;
  int* rowmap;
  int* colmap;
};


extern fdata* 
fdata_create(int datacount, int attcount);

extern fdata* 
fdata_new(fdata* old_dt);

extern fdata* 
fdata_new_cols(fdata* old_dt, int* sel);

extern fdata* 
fdata_new_rows_by_vals(fdata* old_dt, int* sel);

extern void
fdata_read(char* filename, fdata* dt);

extern void 
fdata_old(fdata* dt);

extern void 
fdata_free(fdata* d);

#define FMISSINGVALUE ((float) -2147483648.0)

#define FD(DT,J,I) ((DT)->d[(DT)->rowmap[J] * (DT)->om + (DT)->colmap[I]])

#define FMISSING(DT,J,I) ( FMISSINGVALUE == FD(DT,J,I))

#define FDATA_DEL_ALL(DT) ((DT)->N = 0)
#define FDATA_ADD(DTD,DTS,J) ((DTD)->rowmap[((DTD)->N)++] = DTS->rowmap[J])

#endif


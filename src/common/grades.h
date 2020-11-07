#ifndef GRADES_H
#define GRADES_h

typedef struct grades grades;

struct grades {
  int N;
  int K;
  double* g;
};

#define GJ(gs,j) ((gs)->g + (j) * (gs)->K)
#define GJK(gs,j,k) (GJ(gs,j)[k])

extern grades* grades_create(int N, int K);
extern void    grades_free(grades* gs);
extern void    grades_read(char* filename, grades* gs);
extern void    grades_write(char* filename, grades* gs);
extern grades* grades_cread(int N, int K, char* filename);
extern void    grades_init(gsl_rng* rng, grades* gs);

#endif

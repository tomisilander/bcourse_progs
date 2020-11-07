#include <stdio.h>
#include "err.h"

int main(int argc, char* argv[]){
  FILE* ifp;
  FILE* ofp;
  int c;
  int ret;

  if(3 != argc) {
    fprintf(stderr, "Usage: %s infile outfile\n", argv[0]);
    exit(1);
  }

  OPENFILE_OR_DIE(ifp, argv[1], "r");
  OPENFILE_OR_DIE(ofp, argv[2], "w");

  ret = 0;
  while(EOF != (c = getc(ifp))){
    if(ret && (c != '\n')) putc('\n',ofp);
    ret = (c == '\r');
    if(!ret) putc(c,ofp);
  }
  if(ret) putc('\n',ofp);

  CLOSEFILE_OR_DIE(ifp, argv[1]);
  CLOSEFILE_OR_DIE(ofp, argv[2]);

  return 0;
}


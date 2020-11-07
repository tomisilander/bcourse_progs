#include <stdio.h>
#include <math.h>
#include <signal.h>
#include <errno.h>
#include <unistd.h>
#include "err.h"
#include "format.h"
#include "data.h"
#include "parent_matrix.h"
#include "bane.h"
#include "node.h"


int main(int argc, char* argv[]){
  int n;
  data* dt;
  format* fmt;
  bane* bn;
  FILE* fp;

  if(argc !=6){
    fprintf(stderr, 
	    "Usage: %s structfile vdfile datafile datacount outfile\n",
	    argv[0]);
    exit(-1);
  }

  fmt = format_cread(argv[2]);
  n = atoi(argv[4]);
  dt = data_cread(n,fmt->dim,argv[3]);

  bn = bane_create_from_format(fmt);

  OPENFILE_OR_DIE(fp,argv[1],"r");
  bane_read_structure(bn,fp);
  CLOSEFILE_OR_DIE(fp,argv[1]);
  

  { /* Write the ss */
    int i;

    OPENFILE_OR_DIE(fp,argv[5],"w");
    
    for(i=0; i<bn->nodecount; ++i){
      Word_t* cnt_ptr;
      node* nodi = bn->nodes + i;
      unsigned int pcfg;
      Pvoid_t ss = gather_sparse_ss_for_v(bn, dt, i);
      fprintf(fp,"%d\n",nodi->pcc);
      for(pcfg=0; pcfg<nodi->pcc; ++pcfg){
	int l;
	JLG(cnt_ptr, ss, pcfg);
	for(l=0; l<nodi->valcount; ++l){
	  unsigned int lval = 0;
	  if(cnt_ptr) lval = ((unsigned int*)(*cnt_ptr))[l];
	  fprintf(fp, "%d%c", lval,(l<nodi->valcount-1) ? ' ' : '\n');
	}
      }
      free_sparse_ss(ss);
    }

    CLOSEFILE_OR_DIE(fp,argv[5]);

  }

  bane_free(bn);
  data_free(dt);
  format_free(fmt);

  return 0;
}











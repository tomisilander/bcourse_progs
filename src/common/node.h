#ifndef NODE_H
#define NODE_H

#include <stdio.h>
#include <Judy.h>
#include "typedef_node.h"
#include "bane.h"

struct node {
 
  int id;
  int valcount;
  
  int parentcount;
  Pvoid_t parents;

  int childcount;
  Pvoid_t children;

  int pcc;
  int N;
  int* ss;
  int* ssp;

  int* path_to_me_count;
  int ancestorcount;

};

#define SS(N,P,C) ((N)->ss[(P) * (N)->valcount + (C)])
#define SSP(N,P) ((N)->ssp[P])


#define FIRST_CHILD(NODE, C, NOT_FINISHED) \
{   (C) = 0;\
    J1F((NOT_FINISHED),(NODE)->children,(C));\
}

#define NEXT_CHILD(NODE,C,NOT_FINISHED) \
    J1N((NOT_FINISHED),(NODE)->children,(C))

#define GET_CHILD_N(NODE,C,N) \
{  Word_t Nth = ((N) + 1);\
   int Rc_int;\
   J1BC(Rc_int,  (NODE)->children, Nth, (C));\
}

#define ADD_CHILD(NODE,C) \
{\
  int Retcode;\
  J1S(Retcode, (NODE)->children, C);\
  ++ (NODE)->childcount;\
}

#define DEL_CHILD(NODE,C) \
{\
  int Retcode;\
  J1U(Retcode, (NODE)->children, C);\
  -- (NODE)->childcount;\
}

#define LOOP_OVER_CHILDREN(NODE,C) \
{  int Not_finished; \
   Word_t* Ip = &(C); \
   node* Nodi = (NODE);\
   FIRST_CHILD(Nodi, *Ip, Not_finished);\
   while(Not_finished){

#define END_CHILD_LOOP \
     NEXT_CHILD(Nodi, *Ip, Not_finished);\
   }\
}



#define FIRST_PARENT(NODE, P, NOT_FINISHED) \
{   (P) = 0;\
    J1F((NOT_FINISHED),(NODE)->parents,(P));\
}

#define NEXT_PARENT(NODE,P,NOT_FINISHED) \
    J1N((NOT_FINISHED),(NODE)->parents,(P))

#define GET_PARENT_N(NODE,P,N) \
{  Word_t Nth = ((N) + 1);\
   int Rc_int;\
   J1BC(Rc_int,  (NODE)->parents, Nth, (P));\
}

#define ADD_PARENT(NODE,P,PVC) \
{\
  int Retcode;\
  J1S(Retcode, (NODE)->parents, P);\
  ++ (NODE)->parentcount;\
  (NODE)->pcc *= PVC;\
}

#define DEL_PARENT(NODE,P,PVC) \
{\
  int Retcode;\
  J1U(Retcode, (NODE)->parents, P);\
  -- (NODE)->parentcount;\
  (NODE)->pcc /= PVC;\
}


#define LOOP_OVER_PARENTS(NODE,P) \
{  int Not_finished; \
   Word_t* Ip = &(P); \
   node* Nodi = (NODE);\
   FIRST_PARENT(Nodi, *Ip, Not_finished);\
   while(Not_finished){

#define END_PARENT_LOOP \
     NEXT_PARENT(Nodi, *Ip, Not_finished);\
   }\
}

#define IS_ANCESTOR_OF(NODE,AID) ((NODE)->path_to_me_count[AID]>0)

extern void 
node_assign(node* dst, node* src, int nodecount);

extern void
node_copy_ss(node* dst, node* src);

extern void
node_assign_ss(node* dst, node* src);

extern void 
node_init(node* nodi, int i, int nodecount);

extern void 
node_write(node* nd, int nodecount, FILE* fp);

extern void 
node_write_prob(node* nd, FILE* fp);

#endif


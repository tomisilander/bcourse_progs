#ifndef BANESEARCH_H
#define BANESEARCH_H

#include <stdio.h>
#include "typedef_node.h"
#include "parent_matrix.h"
#include "format.h"
#include "data.h"
#include "node.h"
#include "bane.h"

extern int 
bane_add_random_arc_from(bane* bn, node* from, arc* ar, 
			 int maxtblsize, int commit);

extern int
bane_del_random_arc_from(bane* bn, node* from, arc* ar, 
			 int maxtblsize, int commit);

extern int 
bane_add_random_arc(bane* bn, arc* ar, int maxtblsize, int commit);

extern int
bane_del_random_arc(bane* bn, arc* ar, int maxtblsize, int commit);

extern int
bane_rev_random_arc(bane* bn, arc* ar, int maxtblsize, int commit);

#endif





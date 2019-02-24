#ifndef COMPACT_HEADER
#define COMPACT_HEADER


#include "global_variables.h"
#include "global_functions.h"

using namespace std;


/******************************************************/
void bin_packing_compact_build(instance_data *inst);
/******************************************************/

/******************************************************/
int bin_packing_compact_solve(instance_data *inst);
/******************************************************/

/******************************************************/
void bin_packing_compact_free(instance_data *inst);
/******************************************************/


#endif

#ifndef BP_BRANCH_HEADER
#define BP_BRANCH_HEADER


#include "global_variables.h"
#include "global_functions.h"

#include "BP_PRICER.h"
#include "BP_MASTER.h"

using namespace std;


/*****************************************************************************/
bool branch_select_couple_branching(instance_data *inst,int *item_i,int *item_j);
/*****************************************************************************/

/*****************************************************************************/
void branch_separate(instance_data *inst,int item_i,int item_j);
/*****************************************************************************/

/*****************************************************************************/
void branch_together(instance_data *inst,int item_i,int item_j);
/*****************************************************************************/

/*****************************************************************************/
void branch_couple_free(instance_data *inst,int item_i,int item_j);
/*****************************************************************************/





#endif

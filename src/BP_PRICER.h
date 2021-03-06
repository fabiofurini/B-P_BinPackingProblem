#ifndef BP_PRICER_HEADER
#define BP_PRICER_HEADER


#include "global_variables.h"
#include "global_functions.h"

using namespace std;


/*****************************************************************************/
void pricer_load_cplex(instance_data *inst);
/*****************************************************************************/

/*****************************************************************************/
double pricer_solve_cplex(instance_data *inst);
/*****************************************************************************/

/*****************************************************************************/
void pricer_free_cplex(instance_data *inst);
/*****************************************************************************/


#endif

#ifndef DP_LEX_HEADER
#define DP_LEX_HEADER


#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <string.h>
#include <errno.h>
#include <sstream>
#include <vector>
#include <algorithm>
#include <set>
#include <iomanip> // for precision



#include "global_variables.h"
#include "global_functions.h"


using namespace std;


/*****************************************************************************/
double DP_kp01_advanced(int n, int C, double* p, int* w,double *sol,bool MAXIMALITY,int* order);
/*****************************************************************************/


#endif

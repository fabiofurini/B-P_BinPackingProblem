#ifndef VARIABLE_HEADER
#define VARIABLE_HEADER


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
#include <iomanip>

using namespace std;

#define EPSILON_BP 0.0001

#define MAX_COLUMNS 10000

//////////////////////////////////////////////////////////////////////////////////////////////
//#include </home/fabio/ILOG/CPLEX_Studio_AcademicResearch127/cplex/include/ilcplex/cplex.h>
#include </Users/fabiofurini/Applications/IBM/ILOG/CPLEX_Studio127/cplex/include/ilcplex/cplex.h>
//////////////////////////////////////////////////////////////////////////////////////////////

typedef struct
{


	///////////////////////////////////////////////////////////////////////////////
	char *INPUT_FILE;
	int CPU_NUMBER;
	int TIME_LIMIT;
	int ALGO;
	int OPTION;
	///////////////////////////////////////////////////////////////////////////////

	///////////////////////////////////////////////////////////////////////////////
	int itemNumber;
	int capacity;
	int *itemWeight;
	///////////////////////////////////////////////////////////////////////////////

	///////////////////////////////////////////////////////////////////////////////
	clock_t time_start,time_finish,time_current;
	double computation_time,current_time;
	///////////////////////////////////////////////////////////////////////////////

	/////////////////////////////////////CPLEX/////////////////////////////////////
	CPXENVptr env_MASTER,env_PRICER,env_COMPACT;
	CPXLPptr  lp_MASTER,lp_PRICER,lp_COMPACT;
	int status,ccnt,rcnt,nzcnt,lpstat,nodecount,cur_numrows,cur_numcols;
	int* rmatbeg,*rmatind,*cmatbeg, *cmatind;
	double *rmatval,*cmatval,*rngval,*x,*pi,*obj, *lb, *ub,*rhs,coef_p,objval,bestobjval;
	char *xctype,*sense;
	char **rowname,**colname;
	///////////////////////////////////////////////////////////////////////////////


	///////////////////////////////////////////////////////////////////
	int **status_item_couple; //free 0 together 1 separated -1
	int best_incumbent;
	int BP_nodes;
	int BP_cols;
	double BP_lp;
	double BP_lp_time;
	double BP_pricer_time;
	double BP_master_time;
	double BP_heur_time;
	bool STOP_TIME_LIMIT;
	int *ind;
	double *column;
	int counter_together_branching;
	int counter_separate_branching;
	///////////////////////////////////////////////////////////////////


	///////////////////////////////////////////////////////////////////
	int **A_COLUMN;
	int n_COLUMNS;
	double *x_master;
	///////////////////////////////////////////////////////////////////


} instance_data;




#endif

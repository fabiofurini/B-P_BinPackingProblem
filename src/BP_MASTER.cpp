

#include "BP_MASTER.h"


/////////////////////////////OPTIONS///////////////////////////////
#define set_covering_version
//#define set_partitioning_version

#define cplex_alg_Primal_Simplex
//#define cplex_alg_Dual_Simplex
//#define cplex_alg_Network_Simplex
//#define cplex_alg_Barrier
//#define cplex_alg_Sifting
//#define cplex_alg_Concurrent
///////////////////////////////////////////////////////////////////

/////////////////////////////PRINTINGS/////////////////////////////
//#define print_column
//#define print_pi
#define print_MASTER_columns
//#define print_MASTER
#define verbose_BP
#define verbose_PRUNING
///////////////////////////////////////////////////////////////////



/*****************************************************************************/
int var_feasibility_insert(instance_data *inst)
/*****************************************************************************/
{

	inst->ccnt=1;
	inst->nzcnt=inst->itemNumber;

	inst->obj=(double*) calloc(inst->ccnt,sizeof(double));

	inst->cmatbeg=(int*) calloc(inst->ccnt,sizeof(int));
	inst->cmatind=(int*) calloc(inst->nzcnt,sizeof(int));
	inst->cmatval=(double*) calloc(inst->nzcnt,sizeof(double));

	inst->obj[0]=inst->itemNumber+1;

	/////////////////////////////////////////////////////
	for(int i = 0; i < inst->itemNumber; i++)
	{

		////////////////////////////////////////
		inst->A_COLUMN[inst->n_COLUMNS][i]=1;
		////////////////////////////////////////

		inst->cmatval[i]=1.0;
		inst->cmatind[i]=i;
	}
	/////////////////////////////////////////////////////


	/////////////////////////////////////////////
	inst->n_COLUMNS++;
	/////////////////////////////////////////////

	inst->cmatbeg[0]=0;
	inst->lb=(double*) calloc(inst->ccnt,sizeof(double));
	inst->ub=(double*) calloc(inst->ccnt,sizeof(double));

	inst->lb[0]=0.0;
	inst->ub[0]=CPX_INFBOUND;

	CPXaddcols(inst->env_MASTER,inst->lp_MASTER,inst->ccnt,inst->nzcnt,inst->obj,inst->cmatbeg,inst->cmatind,inst->cmatval,inst->lb,inst->ub,NULL);
	if(inst->status!=0)
	{
		printf("error in CPXaddcols\n");
		exit(-1);
	}

	free(inst->cmatbeg);
	free(inst->cmatval);
	free(inst->cmatind);
	free(inst->lb);
	free(inst->ub);
	free(inst->obj);


	return 1;

}

/*****************************************************************************/
void master_column_add(instance_data *inst)
/*****************************************************************************/
{

	inst->ccnt=1;
	inst->nzcnt=inst->itemNumber;

	inst->obj=(double*) calloc(inst->ccnt,sizeof(double));
	inst->obj[0]=1.0;

	inst->cmatbeg=(int*) calloc(inst->ccnt,sizeof(int));
	inst->cmatind=(int*) calloc(inst->nzcnt,sizeof(int));
	inst->cmatval=(double*) calloc(inst->nzcnt,sizeof(double));

	for(int i = 0; i < inst->nzcnt; i++)
	{
		////////////////////////////////////////
		inst->A_COLUMN[inst->n_COLUMNS][i]=inst->column[i];
		////////////////////////////////////////

		inst->cmatval[i]=inst->column[i];
		inst->cmatind[i]=i;
	}

	/////////////////////////////////////////////
	inst->n_COLUMNS++;
	/////////////////////////////////////////////

	inst->cmatbeg[0]=0;
	inst->lb=(double*) calloc(inst->ccnt,sizeof(double));
	inst->ub=(double*) calloc(inst->ccnt,sizeof(double));

	inst->lb[0]=0.0;
	inst->ub[0]=CPX_INFBOUND;

	CPXaddcols(inst->env_MASTER,inst->lp_MASTER,inst->ccnt,inst->nzcnt,inst->obj,inst->cmatbeg,inst->cmatind,inst->cmatval,inst->lb,inst->ub,NULL);
	if(inst->status!=0)
	{
		printf("error in CPXaddcols\n");
		exit(-1);
	}

	free(inst->cmatbeg);
	free(inst->cmatval);
	free(inst->cmatind);
	free(inst->lb);
	free(inst->ub);
	free(inst->obj);

}

/*****************************************************************************/
void master_build(instance_data *inst)
/*****************************************************************************/
{


	inst->env_MASTER=CPXopenCPLEX(&inst->status);
	if(inst->status!=0)
	{
		printf("cannot open CPLEX environment\n");
		exit(-1);
	}

	inst->lp_MASTER=CPXcreateprob(inst->env_MASTER,&inst->status,"MASTER");
	if(inst->status!=0)
	{
		printf("cannot create problem\n");
		exit(-1);
	}

	CPXchgobjsen(inst->env_MASTER,inst->lp_MASTER,CPX_MIN);


	//create an empty master
	inst->rcnt=inst->itemNumber;
	inst->nzcnt=0;

	inst->rhs =(double*) calloc(inst->rcnt,sizeof(double));
	inst->sense =(char*) calloc(inst->rcnt,sizeof(double));

	for(int i=0;i<inst->rcnt;i++)
	{
		inst->rhs[i]=1.0;

#ifdef 	set_covering_version
		inst->sense[i]= 'G';
#endif
#ifdef set_partitioning_version
		inst->sense[i]= 'E';
#endif

	}

	inst->status=CPXaddrows(inst->env_MASTER,inst->lp_MASTER,0,inst->rcnt,0,inst->rhs,inst->sense,NULL,NULL,NULL,NULL,NULL);
	if(inst->status!=0)
	{
		printf("error in CPXaddrows\n");
		exit(-1);
	}

	free(inst->rhs);
	free(inst->sense);

	cout << "MASTER CREATED...\n";

}

/*****************************************************************************/
double master_solve_lp(instance_data *inst, int level)
/*****************************************************************************/
{

	bool optimality=false;
	double current_LP;
	double lagrangian_bound;
	double val_KP;


	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	do{

		optimality=true;


		double time_start_master=clock();


		inst->status=CPXlpopt(inst->env_MASTER,inst->lp_MASTER);
		if(inst->status!=0)
		{
			printf("error in CPXlpopt\n");
			exit(-1);
		}

		double time_finish_master=clock();
		inst->BP_master_time+=(double)(time_finish_master-time_start_master)/(double)CLOCKS_PER_SEC;


		CPXgetobjval(inst->env_MASTER,inst->lp_MASTER,&current_LP);

#ifdef	print_MASTER_columns
		cout << "current_LP\t" << current_LP << "\t";
#endif

		lagrangian_bound=current_LP;


		inst->status=CPXgetpi(inst->env_MASTER,inst->lp_MASTER,inst->pi,0,inst->itemNumber-1);
		if(inst->status!=0)
		{
			printf("error in CPXgetpi\n");
			exit(-1);
		}

#ifdef print_pi
		for(int j=0;j<inst->itemNumber;j++)
		{
			printf("%.3f\t",inst->pi[j]);
		}
		cout << endl;
#endif

		//solve the KP

		if(level==0)
		{
			//THE DP ONLY WORKS AT THE ROOT NODE
			val_KP=DP_kp01_advanced(inst->itemNumber,inst->capacity,inst->pi,inst->itemWeight,inst->column,false,inst->ind);
		}
		else
		{
			val_KP=pricer_solve_cplex(inst);
		}

#ifdef print_column
		cout << "val_KP\t" << val_KP << endl;
		for(int j=0;j<inst->itemNumber;j++)
		{
			printf("%f\n",inst->column[j]);
		}
		cout << endl;
#endif

		//positive reduced profit
		if(val_KP-1>EPSILON_BP)
		{

			optimality=false;

			master_column_add(inst);

			inst->BP_cols++;

		}

		lagrangian_bound=(lagrangian_bound/val_KP);

		inst->time_finish=clock();
		inst->current_time=(double)(inst->time_finish-inst->time_start)/(double)CLOCKS_PER_SEC;

#ifdef	print_MASTER_columns
		cout <<  "lagrangian_bound\t" << lagrangian_bound << "\t time \t" <<  inst->current_time <<  endl;
#endif

	}while(optimality==false);
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



	///////////////////////////////////////////////////////////////////////////////////////////
	//if the first variable is in the basis the the current LP is infeasible
	int current_status=1;
	double x_feas;
	inst->status=CPXgetmipx(inst->env_MASTER,inst->lp_MASTER,&x_feas,0,0);
	if(inst->status!=0)
	{
		printf("error in CPXgetmipx MASTER\n");
		exit(-1);
	}
	if(x_feas>EPSILON_BP)
	{
		current_status=-1;
	}
	///////////////////////////////////////////////////////////////////////////////////////////


#ifdef	print_MASTER
	master_display_status(inst);
#endif

	if(current_status==-1)
	{
		cout << "MASTER INFESIBLE\n";
		return inst->itemNumber+1;
	}
	else
	{
		return current_LP;
	}

}

/*****************************************************************************/
void master_display_status(instance_data *inst)
/*****************************************************************************/
{

	int n_var=CPXgetnumcols(inst->env_MASTER,inst->lp_MASTER);


	inst->status=CPXgetmipx(inst->env_MASTER,inst->lp_MASTER,inst->x_master,0,n_var-1);
	if(inst->status!=0)
	{
		printf("error in CPXgetmipx MASTER\n");
		exit(-1);
	}

	double coeff;

	cout << "\nCurrent MASTER STATUS\n";
	for(int j=0;j<n_var;j++)
	{

		//if(x_master[j]<EPSILON_BP){continue;}

		printf("%d\t%.3f\t",j,inst->x_master[j]);

		for(int jj=0;jj<inst->itemNumber;jj++)
		{

			inst->status = CPXgetcoef(inst->env_MASTER,inst->lp_MASTER, jj, j, &coeff);
			if (inst->status != 0)
			{
				cout  << "error in CPXgetcoef constraint trasf\n";
				exit(-1);
			}
			cout << coeff;

		}

		double _ub;
		double _lb;

		CPXgetub(inst->env_MASTER,inst->lp_MASTER,&_ub,j,j);
		CPXgetlb(inst->env_MASTER,inst->lp_MASTER,&_lb,j,j);

		cout << "\tbounds\t"<< _lb << "\t" << _ub << endl;

	}


	cout << "Status Item Couples\n";
	for(int i=0;i<inst->itemNumber;i++)
	{
		for(int j=0;j<inst->itemNumber;j++)
		{
			cout << inst->status_item_couple[i][j] << "";
		}
		cout << endl;
	}
	cout << endl;

	cin.get();


}

/*****************************************************************************/
bool master_check_integrality(instance_data *inst)
/*****************************************************************************/
{

	bool integer=true;

	int n_var=CPXgetnumcols(inst->env_MASTER,inst->lp_MASTER);

	inst->status=CPXgetmipx(inst->env_MASTER,inst->lp_MASTER,inst->x_master,0,n_var-1);
	if(inst->status!=0)
	{
		printf("error in CPXgetmipx check_integrality\n");
		exit(-1);
	}


	for (int i = 0; i < n_var; i++)
	{

		if ((fabs(inst->x_master[i] - floor(inst->x_master[i])) > EPSILON_BP)
				&& (fabs(inst->x_master[i] - floor(inst->x_master[i]) - 1.0)
						> EPSILON_BP))
		{
			integer = 0;
			break;
		}

	}

	return integer;
}

/*****************************************************************************/
void master_initialize(instance_data *inst)
/*****************************************************************************/
{

	///////////////////////////////////////////////////
	inst->A_COLUMN=new int*[MAX_COLUMNS];
	for(int i = 0; i < MAX_COLUMNS; i++)
	{
		inst->A_COLUMN[i]=new int[inst->itemNumber];
	}
	inst->n_COLUMNS=0;
	///////////////////////////////////////////////////

	///////////////////////////////////////////////////
	inst->x_master=new double[MAX_COLUMNS];
	///////////////////////////////////////////////////

	inst->STOP_TIME_LIMIT=false;

	inst->ind=new int[inst->itemNumber];
	inst->column=new double[inst->itemNumber];
	inst->pi=new double[inst->itemNumber];

	for(int i=0;i<inst->itemNumber;i++)
	{
		inst->ind[i]=i;
		inst->column[i]=0.0;
	}

	inst->BP_pricer_time=0;
	inst->BP_master_time=0;
	inst->counter_together_branching=0;
	inst->counter_separate_branching=0;

	inst->time_start=clock();
	inst->best_incumbent=9999999;
	inst->BP_nodes=0;
	inst->BP_cols=0;
	inst->BP_lp=-1;

	inst->status_item_couple=new int*[inst->itemNumber];
	for(int i=0;i<inst->itemNumber;i++)
	{
		inst->status_item_couple[i]=new int[inst->itemNumber];
	}

	for(int i=0;i<inst->itemNumber;i++)
	{
		for(int j=0;j<inst->itemNumber;j++)
		{
			inst->status_item_couple[i][j]=0;
		}
	}

	master_build(inst);


	////////////////////////////////////////////////////////
	//initialize KP01 pricer
	pricer_load_cplex(inst);
	////////////////////////////////////////////////////////


	////////////////////////////////////////////////////////
	var_feasibility_insert(inst);
	////////////////////////////////////////////////////////



	////////////////////////////////////////////////////////
	//ADDING SINGLETONS VARIABLES
	for(int i=0;i< inst->itemNumber;i++)
	{
		for(int j=0;j< inst->itemNumber;j++){

			if(j!=i)
			{
				inst->column[j]=0;
			}
			else{
				inst->column[j]=1;
			}
		}

		master_column_add(inst);
	}
	////////////////////////////////////////////////////////




#ifdef cplex_alg_Primal_Simplex
	CPXsetintparam(inst->env_MASTER, CPX_PARAM_LPMETHOD, 1);
#endif

#ifdef cplex_alg_Dual_Simplex
	CPXsetintparam(inst->env_MASTER, CPX_PARAM_LPMETHOD, 2);
#endif

#ifdef cplex_alg_Network_Simplex
	CPXsetintparam(inst->env_MASTER, CPX_PARAM_LPMETHOD, 3);
#endif

#ifdef cplex_alg_Barrier
	CPXsetintparam(inst->env_MASTER, CPX_PARAM_LPMETHOD, 4);
#endif

#ifdef cplex_alg_Sifting
	CPXsetintparam(inst->env_MASTER, CPX_PARAM_LPMETHOD, 5);
#endif

#ifdef cplex_alg_Concurrent
	CPXsetintparam(inst->env_MASTER, CPX_PARAM_LPMETHOD, 6);
#endif

	CPXsetintparam (inst->env_MASTER, CPX_PARAM_THREADS, 1);


}

/*****************************************************************************/
void master_free(instance_data *inst)
/*****************************************************************************/
{


	//////////////////////////////////////////////
	delete []inst->x_master;
	//////////////////////////////////////////////

	//////////////////////////////////////////////
	for(int i = 0; i < MAX_COLUMNS; i++)
	{
		delete []inst->A_COLUMN[i];
	}
	delete []inst->A_COLUMN;
	//////////////////////////////////////////////

	inst->time_finish=clock();
	inst->computation_time=(double)(inst->time_finish-inst->time_start)/(double)CLOCKS_PER_SEC;	/////////////////////////////////////////////////////////////

	int numcols=CPXgetnumcols(inst->env_MASTER,inst->lp_MASTER);
	int	numrows=CPXgetnumrows(inst->env_MASTER,inst->lp_MASTER);


	cout << "\n********************\n";
	cout << "OPT VAL:\t" << inst->best_incumbent<< endl;
	cout << "LP ROOT VAL:\t" << inst->BP_lp<< endl;
	cout << "numcols\t" << numcols << endl;
	cout << "numrows\t" << numrows<< endl;
	cout << "total_time\t" << inst->computation_time << endl;
	cout << "ROOT time\t" << inst->BP_lp_time << endl;
	cout << "PRICING_time (CPLEX)\t" << inst->BP_pricer_time << endl;
	cout << "MASTER_time\t" << inst->BP_master_time << endl;
	cout << "really_generated_cols\t" << inst->BP_cols<< endl;
	cout << "\n********************\n";

	delete[] inst->ind;
	delete[] inst->column;
	delete[] inst->pi;

	///////////////////////////
	pricer_free_cplex(inst);
	///////////////////////////

	inst->status = CPXfreeprob(inst->env_MASTER, &inst->lp_MASTER);
	if (inst->status != 0)
	{
		cout  << "error in CPXfreeprob\n";
		exit(-1);
	}

	inst->status = CPXcloseCPLEX(&inst->env_MASTER);
	if (inst->status != 0)
	{
		cout << "cannot close CPLEX environment\n";
		exit(-1);
	}

	bool status_couple_OK=true;
	for(int i=0;i<inst->itemNumber;i++)
	{
		for(int j=0;j<inst->itemNumber;j++)
		{
			if(inst->status_item_couple[i][j]!=0)
			{
				status_couple_OK=false;
				break;
			}
		}
	}
	if(!status_couple_OK)
	{
		cout << "*************WARNING status_item_couple not clean!!!!*************\n";
	}

	for(int i=0;i<inst->itemNumber;i++)
	{
		delete[] inst->status_item_couple[i];
	}
	delete[] inst->status_item_couple;

}




/*****************************************************************************/
int branch_and_price(instance_data *inst,int level)
/*****************************************************************************/
{

#ifdef	verbose_BP
	cout << "NODE\t" << inst->BP_nodes << "\t" << "\t level \t" << level << endl;
#endif

	///////////////////////////////////////////////////////////////////////////////////////////////////
	if (!inst->STOP_TIME_LIMIT)
	{
		inst->time_finish=clock();
		inst->current_time=(double)(inst->time_finish-inst->time_start)/(double)CLOCKS_PER_SEC;
		if(inst->current_time>inst->TIME_LIMIT)
		{
			inst->STOP_TIME_LIMIT=true;
		}

		if (inst->STOP_TIME_LIMIT)
		{
			return -1;
		}
	}
	else
	{
		return -1;
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////



	///////////////////////////////////////////////////////////////////////////////////////////////////
	//solve the current LP
	inst->BP_nodes++;

	//if infeasible it returns +infinite
	double current_LP=master_solve_lp(inst,level);

	if(level==0)
	{
		inst->BP_lp=current_LP;
		inst->time_finish=clock();
		inst->current_time=(double)(inst->time_finish-inst->time_start)/(double)CLOCKS_PER_SEC;
		inst->BP_lp_time=inst->current_time;

	}

	int dual_bound=my_ceil(current_LP);

#ifdef	verbose_BP
	inst->time_finish=clock();
	inst->current_time=(double)(inst->time_finish-inst->time_start)/(double)CLOCKS_PER_SEC;
	cout << "Level\t" << level << "\tbest_incumbent\t" << inst->best_incumbent << "\tcurrent bound\t"<< current_LP << "\ttime\t"<< inst->current_time << endl;
#endif

	///////////////////////////////////////////////////////////////////////////////////////////////////


	///////////////////////////////////////////////////////////////////////////////////////////////////
	//try to prune and run heuristics
	if(dual_bound>=inst->best_incumbent){

#ifdef	verbose_PRUNING
		cout << "PRUNE!\t" << "level\t" << level << "\tcols\t" << inst->BP_cols << "\tbound\t" << current_LP << "\tbest_incumbent\t" << inst->best_incumbent <<  endl;
#endif

		return 1;

	}
	///////////////////////////////////////////////////////////////////////////////////////////////////


	///////////////////////////////////////////////////////////////////////////////////////////////////
	//check if the solution is integer then return
	bool integer=master_check_integrality(inst);
	if(integer){

		cout << "integer!!\t" << integer << "\tvalue\t" << dual_bound <<  endl;

		if(inst->best_incumbent>dual_bound)
		{
			inst->best_incumbent=dual_bound;
		}

#ifdef	verbose_BP
		cout << "return..." << endl;
#endif

		return 1;
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////


	///////////////////////////////////////////////////////////////////////////////////////////////////
	int item_i;
	int item_j;

	bool _found=branch_select_couple_branching(inst,&item_i,&item_j);

#ifdef	verbose_BP
	cout << "Couple Branching\t" << item_i << "\t" << item_j << endl;
#endif
	///////////////////////////////////////////////////////////////////////////////////////////////////


	///////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef	verbose_BP
	cout << "together branching" << endl;
#endif

	//together branching
	//set master
	branch_together(inst,item_i,item_j);
	branch_and_price(inst,level+1);

	//reset master
	branch_couple_free(inst,item_i,item_j);
	///////////////////////////////////////////////////////////////////////////////////////////////////


	///////////////////////////////////////////////////////////////////////////////////////////////////
	//check bound then return
	if(dual_bound>=inst->best_incumbent)
	{

#ifdef	verbose_BP
		cout << "return..." << endl;
#endif
		return 1;
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////


	///////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef	verbose_BP
	cout << "separate branching\t"  << endl;
#endif

	//separate branching
	//set master
	branch_separate(inst,item_i,item_j);
	branch_and_price(inst,level+1);

	//reset master
	branch_couple_free(inst,item_i,item_j);
	///////////////////////////////////////////////////////////////////////////////////////////////////


#ifdef	verbose_BP
	cout << "return..." << endl;
#endif

	return 1;



}

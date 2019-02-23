

#include "BP_PRICER.h"


//#define CPLEX_OUTPUT
//#define writer_PRICER



/*****************************************************************************/
void pricer_load_cplex(instance_data *inst)
/*****************************************************************************/
{


	inst->env_PRICER=CPXopenCPLEX(&inst->status);
	if(inst->status!=0)
	{
		printf("cannot open CPLEX environment\n");
		exit(-1);
	}

	inst->lp_PRICER=CPXcreateprob(inst->env_PRICER,&inst->status,"PRICER");
	if(inst->status!=0)
	{
		printf("cannot create problem\n");
		exit(-1);
	}


	CPXchgobjsen(inst->env_PRICER,inst->lp_PRICER,CPX_MAX);

	CPXsetintparam (inst->env_PRICER, CPX_PARAM_THREADS, 1);

	//////////////////////////////////////////////////////////////////////////////////////////
	//variables
	inst->ccnt=inst->itemNumber;
	inst->obj=(double*) calloc(inst->ccnt,sizeof(double));
	inst->lb=(double*) calloc(inst->ccnt,sizeof(double));
	inst->ub=(double*) calloc(inst->ccnt,sizeof(double));
	inst->xctype=(char*) calloc(inst->ccnt,sizeof(char));

	for(int i=0; i<inst->ccnt; i++)
	{
		inst->obj[i]=0.0;
		inst->lb[i]=0.0;
		inst->ub[i]=1.0;
		inst->xctype[i]='B';
	}

	inst->status=CPXnewcols(inst->env_PRICER,inst->lp_PRICER,inst->ccnt,inst->obj,inst->lb,inst->ub,inst->xctype,NULL);
	if(inst->status!=0)
	{
		printf("error in CPXnewcols\n");
		exit(-1);
	}


	free(inst->obj);
	free(inst->lb);
	free(inst->ub);
	free(inst->xctype);

	//////////////////////////////////////////////////////////////////////////////////////////


	//////////////////////////////////////////////////////////////////////////////////////////

	//KP constraint
	inst->rcnt=1;
	inst->nzcnt=inst->itemNumber;
	inst->rhs=(double*) calloc(inst->rcnt,sizeof(double));
	inst->sense=(char*) calloc(inst->rcnt,sizeof(double));

	inst->rhs[0]=inst->capacity;
	inst->sense[0]='L';

	inst->rmatbeg=(int*) calloc(inst->rcnt,sizeof(int));
	inst->rmatind=(int*) calloc(inst->nzcnt,sizeof(int));
	inst->rmatval=(double*) calloc(inst->nzcnt,sizeof(double));

	for(int i=0; i<inst->itemNumber; i++)
	{
		inst->rmatval[i]=inst->itemWeight[i];
		inst->rmatind[i]=i;
	}

	inst->rmatbeg[0]=0;

	inst->status=CPXaddrows(inst->env_PRICER,inst->lp_PRICER,0,inst->rcnt,inst->nzcnt,inst->rhs,inst->sense,inst->rmatbeg,inst->rmatind,inst->rmatval,NULL,NULL);
	if(inst->status!=0)
	{
		printf("error in CPXaddrows\n");
		exit(-1);
	}

	free(inst->rmatbeg);
	free(inst->rmatval);
	free(inst->rmatind);
	free(inst->rhs);
	free(inst->sense);

	//////////////////////////////////////////////////////////////////////////////////////////


	cout << "PRICER INITIALIZED...\n";

}

/*****************************************************************************/
double pricer_solve_cplex(instance_data *inst)
/*****************************************************************************/
{

#ifdef CPLEX_OUTPUT
	CPXsetintparam (inst->env_PRICER, CPX_PARAM_SCRIND, CPX_ON);
#endif


	inst->status = CPXchgobj(inst->env_PRICER, inst->lp_PRICER,inst->itemNumber, inst->ind, inst->pi);
	if (inst->status != 0)
	{
		printf("error in CPXchgobj\n");
		exit(-1);
	}


#ifdef writer_PRICER
	inst->status=CPXwriteprob(inst->env_PRICER,inst->lp_PRICER,"pricer.lp",NULL);
	if(inst->status!=0)
	{
		printf("error in CPXwriteprob\n");
		exit(-1);
	}
#endif

	double time_start_pricer=clock();

	inst->status=CPXmipopt(inst->env_PRICER,inst->lp_PRICER);
	if(inst->status!=0)
	{
		printf("error in CPXmipopt\n");
		exit(-1);
	}

	double time_finish_pricer=clock();
	inst->BP_pricer_time+=(double)(time_finish_pricer-time_start_pricer)/(double)CLOCKS_PER_SEC;


	inst->status=CPXgetmipx(inst->env_PRICER,inst->lp_PRICER,inst->column,0,inst->itemNumber-1);
	if(inst->status!=0)
	{
		printf("error in CPXgetmipx\n");
		exit(-1);
	}

	double _objval_p=0;

	inst->status=CPXgetmipobjval(inst->env_PRICER,inst->lp_PRICER,&_objval_p);
	if(inst->status!=0)
	{
		printf("error in CPXgetmipobjval\n");
		exit(-1);
	}


	return _objval_p;

}

/*****************************************************************************/
void pricer_free_cplex(instance_data *inst)
/*****************************************************************************/
{


	inst->status = CPXfreeprob(inst->env_PRICER, &inst->lp_PRICER);
	if (inst->status != 0)
	{
		cout  << "error in CPXfreeprob\n";
		exit(-1);
	}

	inst->status = CPXcloseCPLEX(&inst->env_PRICER);
	if (inst->status != 0)
	{
		cout << "cannot close CPLEX environment\n";
		exit(-1);
	}

}

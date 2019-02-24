

#include "COMPACT.h"

//#define write_COMPACT

#define PRINT_SOL

/******************************************************/
int position_var_x(int item,int bin,int itemNumber, int binNumber)
/******************************************************/
{

	return binNumber*item+bin;

}

/******************************************************/
int position_var_y(int bin,int itemNumber, int binNumber)
/******************************************************/
{

	return binNumber*itemNumber+bin;

}

/******************************************************/
void bin_packing_compact_build(instance_data *inst)
/******************************************************/
{

	inst->env_COMPACT=CPXopenCPLEX(&inst->status);
	if(inst->status!=0)
	{
		printf("cannot open CPLEX environment\n");
		exit(-1);
	}

	inst->lp_COMPACT=CPXcreateprob(inst->env_COMPACT,&inst->status,"COMPACT");
	if(inst->status!=0)
	{
		printf("cannot create problem\n");
		exit(-1);
	}


	int binNumber=inst->itemNumber;

	cout << "fixing the binNumber to the item number\t" << binNumber << endl;


	CPXchgobjsen(inst->env_COMPACT,	inst->lp_COMPACT, CPX_MIN);

	inst->ccnt=inst->itemNumber * binNumber + binNumber;

	inst->obj=(double*) calloc(inst->ccnt,sizeof(double));
	inst->lb=(double*) calloc(inst->ccnt,sizeof(double));
	inst->ub=(double*) calloc(inst->ccnt,sizeof(double));
	inst->xctype=(char*) calloc(inst->ccnt,sizeof(char));

	inst->colname=new char*[inst->ccnt];
	for(int i=0;i<inst->ccnt;i++){inst->colname[i]=new char[100];}

	int dummy=0;
	for(int i=0; i<inst->itemNumber; i++)
	{
		for(int j=0; j<binNumber; j++)
		{
			inst->obj[dummy]=0.0;
			inst->lb[dummy]=0.0;
			inst->ub[dummy]=1.0;
			inst->xctype[dummy]='B';
			sprintf(inst->colname[dummy], "x(%d.%d)",i,j);
			dummy++;
		}
	}

	for(int j=0; j<binNumber; j++)
	{
		inst->obj[dummy]=1.0;
		inst->lb[dummy]=0.0;
		inst->ub[dummy]=1.0;
		inst->xctype[dummy]='B';
		sprintf(inst->colname[dummy], "y(%d)",j);
		dummy++;
	}

	inst->status=CPXnewcols(inst->env_COMPACT,	inst->lp_COMPACT,inst->ccnt,inst->obj,inst->lb,inst->ub,inst->xctype,inst->colname);
	if(inst->status!=0)
	{
		printf("error in CPXnewcols\n");
		exit(-1);
	}

	for(int i=0;i<inst->ccnt;i++)
	{
		delete []inst->colname[i];
	}
	delete [] inst->colname;

	free(inst->obj);
	free(inst->lb);
	free(inst->ub);
	free(inst->xctype);


	for(int j=0; j<binNumber; j++)
	{

		inst->rcnt=1;
		inst->nzcnt=inst->itemNumber+1;

		inst->rhs=(double*) calloc(inst->rcnt,sizeof(double));
		inst->sense=(char*) calloc(inst->rcnt,sizeof(double));

		inst->rhs[0]=0.0;
		inst->sense[0]='L';

		inst->rmatbeg=(int*) calloc(inst->rcnt,sizeof(int));
		inst->rmatind=(int*) calloc(inst->nzcnt,sizeof(int));
		inst->rmatval=(double*) calloc(inst->nzcnt,sizeof(double));

		for(int i=0; i<inst->itemNumber; i++)
		{
			inst->rmatval[i]=inst->itemWeight[i];
			inst->rmatind[i]=position_var_x(i,j,inst->itemNumber,binNumber);

		}

		inst->rmatval[inst->nzcnt-1]=-inst->capacity;
		inst->rmatind[inst->nzcnt-1]=position_var_y(j,inst->itemNumber,binNumber);

		inst->rmatbeg[0]=0;

		inst->status=CPXaddrows(inst->env_COMPACT,	inst->lp_COMPACT,0,inst->rcnt,inst->nzcnt,inst->rhs,inst->sense,inst->rmatbeg,inst->rmatind,inst->rmatval,NULL,NULL);
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

	}

	for(int i=0; i<inst->itemNumber; i++)
	{
		inst->rcnt=1;
		inst->nzcnt=binNumber;

		inst->rhs=(double*) calloc(inst->rcnt,sizeof(double));
		inst->sense=(char*) calloc(inst->rcnt,sizeof(double));

		inst->rhs[0]=1.0;

		inst->sense[0]= 'G';

		inst->rmatbeg=(int*) calloc(inst->rcnt,sizeof(int));
		inst->rmatind=(int*) calloc(inst->nzcnt,sizeof(int));
		inst->rmatval=(double*) calloc(inst->nzcnt,sizeof(double));

		for(int j=0; j<binNumber; j++)
		{
			inst->rmatval[j]=1.0;
			inst->rmatind[j]=position_var_x(i,j,inst->itemNumber,binNumber);
		}

		inst->rmatbeg[0]=0;

		inst->status=CPXaddrows(inst->env_COMPACT,	inst->lp_COMPACT,0,inst->rcnt,inst->nzcnt,inst->rhs,inst->sense,inst->rmatbeg,inst->rmatind,inst->rmatval,NULL,NULL);
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
	}


#ifdef write_COMPACT
	inst->status=CPXwriteprob( inst->env_COMPACT, inst->lp_COMPACT,"COMPACT.lp", NULL ) ;
	if(inst->status!=0)
	{
		printf("errore in CPXwriteprob");
		exit(-1);
	}
#endif


}

/******************************************************/
int bin_packing_compact_solve(instance_data *inst)
/******************************************************/
{

	CPXsetintparam (inst->env_COMPACT, CPX_PARAM_SCRIND, CPX_ON);

	CPXsetintparam (inst->env_COMPACT, CPX_PARAM_THREADS, 1);
	CPXsetdblparam (inst->env_COMPACT, CPX_PARAM_TILIM,inst->TIME_LIMIT);


	inst->status=CPXmipopt(inst->env_COMPACT,inst->lp_COMPACT);
	if(inst->status!=0)
	{
		printf("error in CPXmipopt\n");
	}


	inst->status=CPXgetstat(inst->env_COMPACT,inst->lp_COMPACT);

	cout << "CPXgetstat ->\t" << inst->status << endl;


	double OBJ_VAL=0;
	inst->status=CPXgetobjval(inst->env_COMPACT,inst->lp_COMPACT,&OBJ_VAL);
	if(inst->status!=0) {
		printf("error in CPXgetmipobjval\n");
	}



	cout << "OBJ_VAL ->\t" << OBJ_VAL << endl;


#ifdef PRINT_SOL

	int nn=CPXgetnumcols(inst->env_COMPACT,inst->lp_COMPACT);

	double *x=(double*) calloc(nn,sizeof(double));

	inst->status=CPXgetmipx(inst->env_COMPACT,inst->lp_COMPACT,x,0,nn-1);
	if(inst->status!=0)
	{
		printf("error in CPXgetmipx\n");
		exit(-1);
	}

	int binNumber=inst->itemNumber;

	cout << "fixing the binNumber to the item number\t" << binNumber << endl;

	for(int j=0; j<binNumber; j++)
	{
		if(x[position_var_y(j,inst->itemNumber,binNumber)]>0.5)
		{
			cout << "BIN\t" << j << "\t" << x[position_var_y(j,inst->itemNumber,binNumber)] << "\t\t";

			for(int i=0; i<inst->itemNumber; i++)
			{

				if((int) (x[position_var_x(i,j,inst->itemNumber,binNumber)]+0.5) == 1)
				{
					printf("1");
				}
				else
				{
					printf("0");
				}
			}
			cout << endl;
		}
	}

	free(x);

#endif


	return (int)(OBJ_VAL+0.999);

}

/******************************************************/
void bin_packing_compact_free(instance_data *inst)
/******************************************************/
{
	inst->status = CPXfreeprob(inst->env_COMPACT, &inst->lp_COMPACT);
	if (inst->status != 0)
	{
		cout  << "error in CPXfreeprob\n";
		exit(-1);
	}

	inst->status = CPXcloseCPLEX(&inst->env_COMPACT);
	if (inst->status != 0)
	{
		cout << "cannot close CPLEX environment\n";
		exit(-1);
	}
}



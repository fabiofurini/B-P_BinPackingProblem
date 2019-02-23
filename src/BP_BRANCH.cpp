

#include "BP_BRANCH.h"


/*****************************************************************************/
bool branch_select_couple_branching(instance_data *inst,int *item_i,int *item_j)
/*****************************************************************************/
{

	//SELECT FIRST FRACTIONAL COUPLE!


	int n_var=CPXgetnumcols(inst->env_MASTER,inst->lp_MASTER);

	inst->status=CPXgetmipx(inst->env_MASTER,inst->lp_MASTER,inst->x_master,0,n_var-1);
	if(inst->status!=0)
	{
		printf("error in CPXgetmipx branch_select_couple_branching\n");
		exit(-1);
	}

	for (int ii = 0; ii < inst->itemNumber; ii++)
	{

		for (int jj = ii+1; jj < inst->itemNumber; jj++)
		{

			if(inst->status_item_couple[ii][jj]!=0){continue;}

			double dummy=0;
			double coeff_i;
			double coeff_j;

			for (int j = 0; j < n_var; j++)
			{

//				inst->status = CPXgetcoef(inst->env_MASTER,inst->lp_MASTER, ii, j, &coeff_i);
//				if (inst->status != 0)
//				{
//					cout  << "error in CPXgetcoef constraint trasf\n";
//					exit(-1);
//				}
//
//				inst->status = CPXgetcoef(inst->env_MASTER,inst->lp_MASTER, jj, j, &coeff_j);
//				if (inst->status != 0)
//				{
//					cout  << "error in CPXgetcoef constraint trasf\n";
//					exit(-1);
//				}

				coeff_i=inst->A_COLUMN[j][ii];
				coeff_j=inst->A_COLUMN[j][jj];


				if(coeff_i>0.5 && coeff_j>0.5)
				{
					dummy+=inst->x_master[j];
				}

			}


			if ((fabs(dummy - floor(dummy)) > EPSILON_BP)
					&& (fabs(dummy - floor(dummy) - 1.0)
							> EPSILON_BP))
			{

				*item_i=ii;
				*item_j=jj;

				return true;
			}
			//			}
		}
	}


	return false;


}

/*****************************************************************************/
void branch_separate(instance_data *inst,int item_i,int item_j)
/*****************************************************************************/
{

	inst->counter_separate_branching++;

	inst->status_item_couple[item_i][item_j]=-1;

	int n_var=CPXgetnumcols(inst->env_MASTER,inst->lp_MASTER);


	////////////////////////////////////////////////////////////////////////////////////////

	double coeff_i;
	double coeff_j;

	int ind;
	double d;
	char lu;

	//the first is a feasibility variable
	for (int j = 1; j < n_var; j++)
	{

//		inst->status = CPXgetcoef(inst->env_MASTER,inst->lp_MASTER, item_i, j, &coeff_i);
//		if (inst->status != 0)
//		{
//			cout  << "error in CPXgetcoef constraint\n";
//			exit(-1);
//		}
//
//		inst->status = CPXgetcoef(inst->env_MASTER,inst->lp_MASTER, item_j, j, &coeff_j);
//		if (inst->status != 0)
//		{
//			cout  << "error in CPXgetcoef constraint\n";
//			exit(-1);
//		}

		coeff_i=inst->A_COLUMN[j][item_i];
		coeff_j=inst->A_COLUMN[j][item_j];

		if(coeff_i>0.5 && coeff_j>0.5){

			ind = j;
			d = 0.0;
			lu = 'U';

			inst->status = CPXchgbds(inst->env_MASTER,inst->lp_MASTER, 1, &ind,&lu, &d);
			if (inst->status != 0)
			{
				cout  << "error in CPXchgobj \n";
				exit(-1);
			}

		}


	}

	////////////////////////////////////////////////////////////////////////////////////////



	////////////////////////////////////////////////////////////////////////////////////////

	inst->rcnt=1;
	inst->nzcnt=2;
	inst->rhs=(double*) calloc(inst->rcnt,sizeof(double));
	inst->sense=(char*) calloc(inst->rcnt,sizeof(double));

	inst->rhs[0]=1.0;
	inst->sense[0]='L';

	inst->rmatbeg=(int*) calloc(inst->rcnt,sizeof(int));
	inst->rmatind=(int*) calloc(inst->nzcnt,sizeof(int));
	inst->rmatval=(double*) calloc(inst->nzcnt,sizeof(double));

	inst->rmatval[0]=1.0;
	inst->rmatind[0]=item_i;

	inst->rmatval[1]=1.0;
	inst->rmatind[1]=item_j;

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

	////////////////////////////////////////////////////////////////////////////////////////


}

/*****************************************************************************/
void branch_together(instance_data *inst,int item_i,int item_j)
/*****************************************************************************/
{
	inst->counter_together_branching++;

	inst->status_item_couple[item_i][item_j]=1;

	int n_var=CPXgetnumcols(inst->env_MASTER, inst->lp_MASTER);

	////////////////////////////////////////////////////////////////////////////////////////

	double coeff_i;
	double coeff_j;

	int ind;
	double d;
	char lu;

	//the first is a feasibility variable
	for (int j = 1; j < n_var; j++)
	{

//		inst->status = CPXgetcoef(inst->env_MASTER, inst->lp_MASTER, item_i, j, &coeff_i);
//		if (inst->status != 0)
//		{
//			cout  << "error in CPXgetcoef constraint trasf\n";
//			exit(-1);
//		}
//
//		inst->status = CPXgetcoef(inst->env_MASTER, inst->lp_MASTER, item_j, j, &coeff_j);
//		if (inst->status != 0)
//		{
//			cout  << "error in CPXgetcoef constraint trasf\n";
//			exit(-1);
//		}

		coeff_i=inst->A_COLUMN[j][item_i];
		coeff_j=inst->A_COLUMN[j][item_j];

		if((coeff_i>0.5 && coeff_j<0.5) || (coeff_i<0.5 && coeff_j>0.5)){

			ind = j;
			d = 0.0;
			lu = 'U';

			inst->status = CPXchgbds(inst->env_MASTER, inst->lp_MASTER, 1, &ind,&lu, &d);
			if (inst->status != 0)
			{
				cout  << "error in CPXchgobj \n";
				exit(-1);
			}

		}


	}

	////////////////////////////////////////////////////////////////////////////////////////



	////////////////////////////////////////////////////////////////////////////////////////

	inst->rcnt=1;
	inst->nzcnt=2;
	inst->rhs=(double*) calloc(inst->rcnt,sizeof(double));
	inst->sense=(char*) calloc(inst->rcnt,sizeof(double));

	inst->rhs[0]=0.0;
	inst->sense[0]='E';

	inst->rmatbeg=(int*) calloc(inst->rcnt,sizeof(int));
	inst->rmatind=(int*) calloc(inst->nzcnt,sizeof(int));
	inst->rmatval=(double*) calloc(inst->nzcnt,sizeof(double));

	inst->rmatval[0]=1.0;
	inst->rmatind[0]=item_i;

	inst->rmatval[1]=-1.0;
	inst->rmatind[1]=item_j;

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

	////////////////////////////////////////////////////////////////////////////////////////



}

/*****************************************************************************/
void branch_couple_free(instance_data *inst,int item_i,int item_j)
/*****************************************************************************/
{

	inst->status_item_couple[item_i][item_j]=0;

	int n_var=CPXgetnumcols(inst->env_MASTER,inst->lp_MASTER);

	////////////////////////////////////////////////////////////////////////////////////////

	double coeff_i;
	double coeff_j;
	double dummi;

	int ind;
	double d;
	char lu;

	//the first is a feasibility variable
	for (int j = 1; j < n_var; j++)
	{

		bool freedom=true;

		CPXgetub(inst->env_MASTER,inst->lp_MASTER, &dummi, j,j);

		if(dummi > 0.5){continue;}


		for (int ii = 0; ii < inst->itemNumber && freedom; ii++)
		{

			for (int jj = ii+1; jj < inst->itemNumber && freedom; jj++)
			{

//				inst->status = CPXgetcoef(inst->env_MASTER,inst->lp_MASTER, ii, j, &coeff_i);
//				if (inst->status != 0)
//				{
//					cout  << "error in CPXgetcoef constraint trasf\n";
//					exit(-1);
//				}
//
//				inst->status = CPXgetcoef(inst->env_MASTER,inst->lp_MASTER, jj, j, &coeff_j);
//				if (inst->status != 0)
//				{
//					cout  << "error in CPXgetcoef constraint trasf\n";
//					exit(-1);
//				}

				coeff_i=inst->A_COLUMN[j][ii];
				coeff_j=inst->A_COLUMN[j][jj];

				if( (  coeff_i>0.5 && coeff_j>0.5) && inst->status_item_couple[ii][jj]==-1)
				{

					freedom=false;

				}

				if( ( (coeff_i>0.5 && coeff_j<0.5) || (coeff_i<0.5 && coeff_j>0.5) ) && (inst->status_item_couple[ii][jj]==1) )
				{

					freedom=false;

				}

			}
		}
		if(freedom)
		{
			ind = j;
			d = CPX_INFBOUND;
			lu = 'U';

			inst->status = CPXchgbds(inst->env_MASTER,inst->lp_MASTER, 1, &ind,&lu, &d);
			if (inst->status != 0)
			{
				cout  << "error in CPXchgobj \n";
				exit(-1);
			}

		}

	}

	////////////////////////////////////////////////////////////////////////////////////////

	int n_row=CPXgetnumrows(inst->env_PRICER,inst->lp_PRICER);

	inst->status = CPXdelrows(inst->env_PRICER,inst->lp_PRICER,n_row-1,n_row-1);

	if (inst->status != 0)
	{
		cout  << "error in CPXchgobj \n";
		exit(-1);
	}

	////////////////////////////////////////////////////////////////////////////////////////

}

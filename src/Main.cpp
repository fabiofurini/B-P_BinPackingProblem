/*
 *		Created on: 20/02/2019
 *		Author: Fabio Furini
 */

#include <stdlib.h>
#include <stdio.h>

using namespace std;

#include "BP_MASTER.h"
#include "BP_PRICER.h"
#include "BP_BRANCH.h"
#include "COMPACT.h"
#include "global_functions.h"
#include "global_variables.h"

/*****************************************************************************/
int main(int argc, char** argv)
/*****************************************************************************/
{
	instance_data inst;

	inst.INPUT_FILE = (char *) calloc(1000, sizeof(char));

	inst.CPU_NUMBER=1;

	if (argc == 5)
	{
		strcpy(inst.INPUT_FILE, argv[1]);
		inst.TIME_LIMIT=atoi(argv[2]);
		inst.ALGO=atoi(argv[3]);
		inst.OPTION=atoi(argv[4]);

	}
	else
	{
		cout << "ERROR NUMBER STANDARD PARAMETERS -- 1) instance 2) time limit 3) option" << endl;
		exit(2);
	}

	cout <<         "---------------------------------------------" << endl ;
	cout << "INPUT_FILE\t->" << inst.INPUT_FILE << endl;
	cout << "TIME_LIMIT\t->" << inst.TIME_LIMIT << endl;
	cout << "ALGO\t->" << inst.ALGO << endl;
	cout << "OPTION\t->" << inst.OPTION << endl;
	cout <<         "---------------------------------------------" << endl ;

	/////////////////////////////////////
	bin_packing_instance_read(&inst);
	/////////////////////////////////////

	//////////////////////////////////////////////////////////////////////////
	if(inst.ALGO==1)
	{

		cout << "\n\n**************BRANCH-AND-PRICE...\n\n";

		master_initialize(&inst);

		cout << "MASTER INITIALIZED...\n";

		branch_and_price(&inst,0);

		master_free(&inst);

		cout << "MASTER FREE...\n";
	}
	//////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////////////////////////////////
	if(inst.ALGO==2)
	{

		cout << "\n\n**************COMPACT FORMULATION...\n\n";

		bin_packing_compact_build(&inst);

		cout << "FORMULATION INITIALIZED...\n";

		bin_packing_compact_solve(&inst);


		bin_packing_compact_free(&inst);


		cout << "FORMULATION FREE...\n";
	}
	//////////////////////////////////////////////////////////////////////////



	/////////////////////////////////////
	bin_packing_instance_free(&inst);
	/////////////////////////////////////

	free(inst.INPUT_FILE);

	return 1;
}


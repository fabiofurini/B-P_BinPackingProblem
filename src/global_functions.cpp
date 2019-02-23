

#include "global_functions.h"

//#define print_instance


/*****************************************************************************/
int my_floor(double a)
/*****************************************************************************/
{
	int q;
	q = (int) (a + EPSILON_BP);
	return (q);
}

/*****************************************************************************/
int my_ceil(double a)
/*****************************************************************************/
{
	int q;
	double b;
	b = a + 1.0 - EPSILON_BP;
	q = (int) (b);
	return (q);
}

/*****************************************************************************/
void bin_packing_instance_free(instance_data *inst)
/*****************************************************************************/
{

	delete [] inst->itemWeight;

}


/*****************************************************************************/
void bin_packing_instance_read(instance_data *inst)
/*****************************************************************************/
{


	ifstream in(inst->INPUT_FILE);

	if(in.fail())
	{
		cout << "Instance reading ERROR" << endl ;
		exit(1);
	}

	in >> inst->itemNumber;
	in >> inst->capacity;


	inst->itemWeight=new int[inst->itemNumber];


	for (int i = 0; i < inst->itemNumber; i++)
	{
		in >> inst->itemWeight[i];

	}

	in.close();

	cout << "ItemTypeNumber\t\t" << inst->itemNumber << endl;
	cout << "Capacity\t\t" << inst->capacity << endl;

#ifdef print_instance
	for (int i = 0; i < inst->itemNumber; i++)
	{
		cout << "w\t" << inst->itemWeight[i] << "\n";
	}
#endif

	cout << "END READING...\n";


}






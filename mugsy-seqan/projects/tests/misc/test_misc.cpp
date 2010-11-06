#include <iostream>
#include <fstream>
#include <sstream>
#include <typeinfo>
#include <ctime>
#define SEQAN_DEBUG
#define SEQAN_TEST

#include <seqan/basic.h>

#include <seqan/misc/misc_random.h>
#include <seqan/map.h>

#include <seqan/misc/misc_set.h>
//#include <seqan/misc/misc_map.h> geht nicht
#include <seqan/misc/priority_type_base.h>
#include <seqan/misc/priority_type_heap.h>

using namespace std;
using namespace seqan;

//////////////////////////////////////////////////////////////////////////////

void testRandom()
{
	mtRandInit();

	for (unsigned int i=0; i<100; ++i)
	{
		cout << mtRand() << ", ";
	}
	cout << "\n\n";

	for (unsigned int i=0; i<100; ++i)
	{
		cout << geomRand<int>() << ", ";
	}

}

//////////////////////////////////////////////////////////////////////////////

int main() 
{
	SEQAN_TREPORT("TEST BEGIN")
//____________________________________________________________________________

	testRandom();

//____________________________________________________________________________

	SEQAN_TREPORT("TEST END")

	return 0;
}

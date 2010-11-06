#define SEQAN_DEBUG
#define SEQAN_TEST
#define SEQAN_NOSRAN

#include <seqan/basic.h>

using namespace seqan;

//////////////////////////////////////////////////////////////////////////////

void Main_Test_Map();
void Main_TestSumlist();

//////////////////////////////////////////////////////////////////////////////

int main() 
{
	SEQAN_TREPORT("TEST BEGIN")
//____________________________________________________________________________

    Main_Test_Map();
    Main_TestSumlist();

//____________________________________________________________________________

	SEQAN_TREPORT("TEST END")

	return 0;
}

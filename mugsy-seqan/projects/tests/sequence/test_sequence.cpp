#include <iostream>
#include <fstream>
#include <typeinfo>
#include <time.h>

#define SEQAN_DEBUG
#define SEQAN_TEST

#include <seqan/basic.h>

using namespace std;
using namespace seqan;

int mainTestString();
int mainTestStringSet();
int mainTestSegment();

//////////////////////////////////////////////////////////////////////////////

int main()  
{
	SEQAN_TREPORT("TEST BEGIN")

	mainTestString();
	mainTestStringSet();
	mainTestSegment();

	SEQAN_TREPORT("TEST END")

	return 0;
}

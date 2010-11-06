#define SEQAN_DEBUG
#define SEQAN_TEST
//#define SEQAN_VVERBOSE

#define TEST_PATH "projects/tests/blast/"
#define LIB_PATH "projects/library/seqan/blast/"
//#define TEST_PATH "C:\\seqan\\projects\\tests\\blast\\"
//#define LIB_PATH "C:\\seqan\\projects\\library\\seqan\\blast\\"



#include <seqan/file.h>
#include <seqan/blast.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <typeinfo>
#include <time.h>
#include <cstdio>
#include <vector>
#include <time.h>
#include <string>
#include <ctime>


//#include "test_blast_library.h"
#include "test_blast_parsing.h"

using namespace std;
using namespace seqan;


//////////////////////////////////////////////////////////////////////////////

int main() 
{
	SEQAN_TREPORT("TEST BEGIN")

	Test_BlastStoreReport<int>();
	Test_BlastStoreReportBasic<int>();
	Test_BlastParsing<int>(BlastN());
	Test_BlastParsing<int>(BlastP());
	Test_BlastParsingBasic<int>(BlastN());
	Test_BlastParsingBasic<int>(BlastP());
	
	debug::verifyCheckpoints(LIB_PATH "blast_parsing.h");
	debug::verifyCheckpoints(LIB_PATH "blast_base.h");
	debug::verifyCheckpoints(LIB_PATH "blast_report.h");
	debug::verifyCheckpoints(LIB_PATH "blast_hit.h");
	debug::verifyCheckpoints(LIB_PATH "blast_hsp.h");
	debug::verifyCheckpoints(LIB_PATH "blast_stream_report.h");
	debug::verifyCheckpoints(LIB_PATH "blast_stream_hit.h");

	debug::verifyCheckpoints(LIB_PATH "blast_iterator.h");
	debug::verifyCheckpoints(LIB_PATH "blast_hit_iterator.h");
	debug::verifyCheckpoints(LIB_PATH "blast_hsp_iterator.h");
	debug::verifyCheckpoints(LIB_PATH "blast_stream_hit_iterator.h");
	debug::verifyCheckpoints(LIB_PATH "blast_stream_hsp_iterator.h");


	SEQAN_TREPORT("TEST END")

	return 0;
}

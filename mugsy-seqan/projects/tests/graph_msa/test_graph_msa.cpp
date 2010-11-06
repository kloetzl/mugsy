#define SEQAN_DEBUG
#define SEQAN_TEST
#define SEQAN_VERBOSE


// Test path
#define TEST_PATH "projects/tests/graph_msa/"

// External / STL
#include <iostream>
#include <fstream>
#include <string>

//SeqAn
#include <seqan/graph_msa.h>
#include <seqan/misc/misc_random.h>

// Test files
#include "test_graph_tcoffee.h"


using namespace seqan;

//////////////////////////////////////////////////////////////////////////////

int main() 
{
	SEQAN_TREPORT("TEST BEGIN")

	// Redirect std::cout
    std::ofstream file(TEST_PATH "redirect.txt");
    std::streambuf* strm_puffer = std::cout.rdbuf();
    std::cout.rdbuf(file.rdbuf());

	Test_GraphTCoffee();		// Test T-Coffee

	// Restore std::cout
	std::cout.rdbuf(strm_puffer);

	SEQAN_TREPORT("TEST END")

	return 0;
}

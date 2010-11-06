#define SEQAN_DEBUG
#define SEQAN_TEST
#define SEQAN_VERBOSE

// Test path
#define TEST_PATH "projects/tests/consensus/"

// External / STL
#include <iostream>
#include <fstream>
#include <string>

//SeqAn
#include <seqan/consensus.h>

// Test files
#include "test_consensus.h"


using namespace seqan;

//////////////////////////////////////////////////////////////////////////////

int main() 
{
	SEQAN_TREPORT("TEST BEGIN")

	// Redirect std::cout
    std::ofstream file(TEST_PATH "redirect.txt");
    std::streambuf* strm_puffer = std::cout.rdbuf();
    std::cout.rdbuf(file.rdbuf());

	Test_Consensus();		// Test Graph Consensus

	// Restore std::cout
	std::cout.rdbuf(strm_puffer);

	SEQAN_TREPORT("TEST END")

	return 0;
}

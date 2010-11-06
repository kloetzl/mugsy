#define SEQAN_DEBUG
#define SEQAN_TEST
#define SEQAN_VERBOSE

// Test path
#define TEST_PATH "projects/tests/graph_align/"

// External / STL
#include <iostream>
#include <fstream>
#include <string>

#include <seqan/map.h>


// SeqAn Includes
#include <seqan/graph_align.h>


// Test files
#include "test_graph_align.h"


using namespace seqan;


//////////////////////////////////////////////////////////////////////////////

int main() 
{
	SEQAN_TREPORT("TEST BEGIN")
	//simulateSequences();

	// Redirect std::cout
    std::ofstream file(TEST_PATH "redirect.txt");
    std::streambuf* strm_puffer = std::cout.rdbuf();
    std::cout.rdbuf(file.rdbuf());
	
	Test_GraphAlignment();		// Test Graph Alignment

	// Restore std::cout
	std::cout.rdbuf(strm_puffer);

	SEQAN_TREPORT("TEST END")

	return 0;
}

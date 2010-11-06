#define SEQAN_DEBUG
#define SEQAN_TEST
#define SEQAN_VVERBOSE

// Test path
#define TEST_PATH "projects/tests/graph_types/"


// External STL
#include <iostream>
#include <fstream>
#include <string>


// Seqan
#include <seqan/graph_types.h>

// Test files
#include "test_graph_basic.h"
#include "test_graph_types.h"
#include "test_graph_iterators.h"
#include "test_graph_properties.h"
#include "test_graph_derived.h"
#include "test_graph_utils.h"


using namespace seqan;

//////////////////////////////////////////////////////////////////////////////

int main() 
{
	SEQAN_TREPORT("TEST BEGIN")

	// Redirect std::cout
    std::ofstream file(TEST_PATH "redirect.txt");
    std::streambuf* strm_puffer = std::cout.rdbuf();
    std::cout.rdbuf(file.rdbuf());

	// Execute Tests
	Test_GraphBasics();			// Test Graph Basic
	Test_GraphTypes();			// Test Graph Types
	Test_GraphIterators();		// Test Graph Iterators
	Test_GraphProperties();		// Test internal and external property maps
	Test_GraphDerivedTypes();	// Test Additional graph types, e.g., oracle, trie,...
	Test_GraphUtils();			// Test the graph drawing, etc.

	// Restore std::cout
	std::cout.rdbuf(strm_puffer);

	debug::verifyCheckpoints("projects/library/seqan/graph_types/graph_interface.h");

	
	SEQAN_TREPORT("TEST END")

	return 0;
}

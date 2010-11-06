#define SEQAN_DEBUG
#define SEQAN_TEST
#define SEQAN_VERBOSE

// Test path
#define TEST_PATH "projects/tests/refinement/"
#define LIB_PATH "projects/library/seqan/refinement/"

//#define LIB_PATH "/home/bude2/emde/seqan/projects/library/seqan/refinement/"
//#define TEST_PATH "/home/bude2/emde/seqan/projects/tests/refinement/"

// External / STL
#include <iostream>
#include <fstream>
#include <string>

// SeqAn Includes
#include <seqan/refinement.h>

// Test files
#include "test_graph_impl_align.h"
#include "test_graph_match_refinement.h"
#include "test_graph_interval_tree.h"


using namespace seqan;


//////////////////////////////////////////////////////////////////////////////

int main() 
{
	SEQAN_TREPORT("TEST BEGIN")

	// Redirect std::cout
    std::ofstream file(TEST_PATH "redirect.txt");
    std::streambuf* strm_puffer = std::cout.rdbuf();
    std::cout.rdbuf(file.rdbuf());

	Test_AlignmentGraph();
	Test_GraphMatchRefinement();// Test Match Refinement
	Test_GraphIntervalTree();	// Test Interval Tree

	// Restore std::cout
	std::cout.rdbuf(strm_puffer);


	debug::verifyCheckpoints(LIB_PATH "graph_impl_interval_types.h");
	debug::verifyCheckpoints(LIB_PATH "graph_impl_interval_tree.h");
	debug::verifyCheckpoints(LIB_PATH "graph_algorithm_refine_scoring.h");
	debug::verifyCheckpoints(LIB_PATH "graph_algorithm_refine_fragment.h");
	debug::verifyCheckpoints(LIB_PATH "graph_algorithm_refine_exact.h");
//	debug::verifyCheckpoints(LIB_PATH "graph_algorithm_refine_aligngraph.h");
//	debug::verifyCheckpoints(LIB_PATH "graph_algorithm_refine_align.h");
//	debug::verifyCheckpoints(LIB_PATH "graph_algorithm_refine_annotation.h");
//	debug::verifyCheckpoints(LIB_PATH "graph_algorithm_refine_inexact.h");

	SEQAN_TREPORT("TEST END")

	return 0;
}

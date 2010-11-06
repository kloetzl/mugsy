#define SEQAN_DEBUG
#define SEQAN_TEST
#define SEQAN_VERBOSE

// SeqAn
#include <seqan/graph_algorithms.h>
#include "test_graph_algorithms.h"

// SeqAn Namespace
using namespace seqan;

//////////////////////////////////////////////////////////////////////////////

int main() 
{
	SEQAN_TREPORT("TEST BEGIN")

	Test_GraphAlgorithms();		// Test Graph Algorithms

	SEQAN_TREPORT("TEST END")

	return 0;
}

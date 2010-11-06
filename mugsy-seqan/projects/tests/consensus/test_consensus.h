#ifndef SEQAN_HEADER_TEST_CONSENSUS_H
#define SEQAN_HEADER_TEST_CONSENSUS_H

namespace SEQAN_NAMESPACE_MAIN
{



//////////////////////////////////////////////////////////////////////////////


void Test_Consensus() {



	debug::verifyCheckpoints("projects/library/seqan/consensus/graph_consensus_base.h");
	debug::verifyCheckpoints("projects/library/seqan/consensus/graph_consensus_library.h");
	debug::verifyCheckpoints("projects/library/seqan/consensus/graph_consensus_ctgstore.h");
	debug::verifyCheckpoints("projects/library/seqan/consensus/graph_consensus_frgstore.h");
	debug::verifyCheckpoints("projects/library/seqan/consensus/graph_consensus_libstore.h");
	debug::verifyCheckpoints("projects/library/seqan/consensus/graph_consensus_readstore.h");
}


}

#endif


// Projekt, mit dem die Demos getestet werden koennen

#define TEST_BASICS
#define TEST_FIND
#define TEST_ALIGNMENT
#define TEST_MODIFIERS
#define TEST_INDEX
#define TEST_GRAPH
#define TEST_MISCELLANEOUS


#ifdef TEST_BASICS

	#define main runAllocator
	#include "../../library/demos/allocator.cpp"
	#undef main

	#define main runAlphabet
	#include "../../library/demos/alphabet.cpp"
	#undef main


	#define main runIterator
	#include "../../library/demos/iterator.cpp"
	#undef main


	#define main runRootedIterator
	#include "../../library/demos/rooted_iterator.cpp"
	#undef main

	#define main runFileFormat
	#include "../../library/demos/file_format.cpp"
	#undef main

#endif

#ifdef TEST_FIND

	#define main runFindExact
	#include "../../library/demos/find_exact.cpp"
	#undef main

	#define main runFindApprox
	#include "../../library/demos/find_approx.cpp"
	#undef main

	#define main runFindWild
	#include "../../library/demos/find_wild.cpp"
	#undef main

	#define main runFindMotif
	#include "../../library/demos/find_motif.cpp"
	#undef main

#endif

#ifdef TEST_ALIGNMENT

	#define main runAlignment
	#include "../../library/demos/alignment.cpp"
	#undef main

	#define main runMsaAlignment
	#include "../../library/demos/alignment_msa.cpp"
	#undef main

	#define main runAlignmentLocal
	#include "../../library/demos/alignment_local.cpp"
	#undef main

#endif

#ifdef TEST_MODIFIERS

	#define main runModifierModReverse
	#include "../../library/demos/modifier_modreverse.cpp"
	#undef main

	#define main runModifierModView
	#include "../../library/demos/modifier_modview.cpp"
	#undef main

	#define main runModifierNested
	#include "../../library/demos/modifier_nested.cpp"
	#undef main

#endif

#ifdef TEST_INDEX

	#define main runIndexSA
	#include "../../library/demos/index_sufarray.cpp"
	#undef main

	#define main runIndexFind
	#include "../../library/demos/index_find.cpp"
	#undef main

	#define main runIndexFindStringSet
	#include "../../library/demos/index_find_stringset.cpp"
	#undef main

	#define main runIndexMUMs
	#include "../../library/demos/index_mums.cpp"
	#undef main

	#define main runIndexSuperMaxRepeats
	#include "../../library/demos/index_supermaxrepeats.cpp"
	#undef main

	#define main runIndexMaxRepeats
	#include "../../library/demos/index_maxrepeats.cpp"
	#undef main

	#define main runIndexMummy
	#include "../../library/demos/index_mummy.cpp"
	#undef main

	#define main runIndexNodePredicate
	#include "../../library/demos/index_node_predicate.cpp"
	#undef main

#endif

#ifdef TEST_GRAPH

	#define main runGraphBfs
	#include "../../library/demos/graph_algo_bfs.cpp"
	#undef main

	#define main runGraphDfs
	#include "../../library/demos/graph_algo_dfs.cpp"
	#undef main

	#define main runGraphTopsort
	#include "../../library/demos/graph_algo_topsort.cpp"
	#undef main

	#define main runGraphScc
	#include "../../library/demos/graph_algo_scc.cpp"
	#undef main

	#define main runGraphPrim
	#include "../../library/demos/graph_algo_tree_prim.cpp"
	#undef main

	#define main runGraphKruskal
	#include "../../library/demos/graph_algo_tree_kruskal.cpp"
	#undef main

	#define main runGraphDag
	#include "../../library/demos/graph_algo_path_dag.cpp"
	#undef main

	#define main runGraphBellman
	#include "../../library/demos/graph_algo_path_bellmanford.cpp"
	#undef main

	#define main runGraphDijkstra
	#include "../../library/demos/graph_algo_path_dijkstra.cpp"
	#undef main

	#define main runGraphAllpairs
	#include "../../library/demos/graph_algo_path_allpairs.cpp"
	#undef main

	#define main runGraphFloyd
	#include "../../library/demos/graph_algo_path_floydwarshall.cpp"
	#undef main

	#define main runGraphTransitive
	#include "../../library/demos/graph_algo_path_transitive.cpp"
	#undef main

	#define main runGraphFord
	#include "../../library/demos/graph_algo_flow_fordfulkerson.cpp"
	#undef main

	#define main runGraphLis
	#include "../../library/demos/graph_algo_lis.cpp"
	#undef main

	#define main runGraphHis
	#include "../../library/demos/graph_algo_his.cpp"
	#undef main

	#define main runGraphLcs
	#include "../../library/demos/graph_algo_lcs.cpp"
	#undef main

	#define main runGraphHmm
	#include "../../library/demos/graph_hmm.cpp"
	#undef main

	#define main runGraphHmmSilent
	#include "../../library/demos/graph_hmm_silent.cpp"
	#undef main

#endif

#ifdef TEST_MISCELLANEOUS

	#define main runBlastReport
	#include "../../library/demos/blast_report.cpp"
	#undef main

	#define main runSeeds
	#include "../../library/demos/seeds.cpp"
	#undef main

	#define main runLagan
	#include "../../library/demos/lagan.cpp"
	#undef main

#endif

int main(int argc, const char *argv[]) 
{
#ifdef TEST_BASICS
	runAllocator();
	runAlphabet();
	runIterator();
	runRootedIterator();
	runFileFormat();
#endif

#ifdef TEST_FIND
	runFindExact();
	runFindApprox();
	runFindWild();
	runFindMotif();
#endif

#ifdef TEST_ALIGNMENT
	runAlignment();
	runMsaAlignment();
	runAlignmentLocal();
#endif

#ifdef TEST_MODIFIERS
	runModifierModReverse();
	runModifierModView();
	runModifierNested();
#endif

#ifdef TEST_INDEX
	runIndexSA();
	runIndexFind();
	runIndexFindStringSet();
	runIndexMUMs();
	runIndexSuperMaxRepeats();
	runIndexMaxRepeats();
	runIndexMummy(argc, argv);
	runIndexNodePredicate();
#endif

#ifdef TEST_GRAPH
	runGraphBfs();
	runGraphDfs();
	runGraphTopsort();
	runGraphScc();
	runGraphPrim();
	runGraphKruskal();
	runGraphDag();
	runGraphBellman();
	runGraphDijkstra();
	runGraphAllpairs();
	runGraphFloyd();
	runGraphTransitive();
	runGraphFord();
	runGraphLis();
	runGraphHis();
	runGraphLcs();
	runGraphHmm();
	runGraphHmmSilent();
#endif

#ifdef TEST_MISCELLANEOUS
	runBlastReport();
	runSeeds();
	runLagan(argc, argv);
#endif

	return 0;
}

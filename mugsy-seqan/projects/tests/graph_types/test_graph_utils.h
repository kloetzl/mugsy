#ifndef SEQAN_HEADER_TEST_GRAPH_UTILS_H
#define SEQAN_HEADER_TEST_GRAPH_UTILS_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////

template <typename TGraph>
void  Test_GraphDrawing_Tmp() {
	// Create a dummy graph
	TGraph g;
	addVertex(g);addVertex(g);addEdge(g,0,1);
	addVertex(g);addEdge(g,0,2);addVertex(g); 
	addVertex(g);addVertex(g);addEdge(g,3,4);

	// Dot Drawing
	write(std::cout,g,DotDrawing());
}

//////////////////////////////////////////////////////////////////////////////

void  Test_GraphDrawing() {
	Test_GraphDrawing_Tmp<Graph<Directed<> > >();
	Test_GraphDrawing_Tmp<Graph<Undirected<> > >();
	Test_GraphDrawing_Tmp<Graph<Tree<> > >();

	// Automat
	Graph<Automaton<Dna> > automat;
	createRoot(automat);addVertex(automat);addEdge(automat,0,1,'a');
	addVertex(automat);addEdge(automat,0,2,'g');
	// Dot Drawing
	write(std::cout,automat,DotDrawing());

	// Trie
	Graph<Automaton<char> > trie;
	String<String<unsigned int> > pos;
	String<String<char> > keywords;
	appendValue(keywords, String<char>("announce"));
	appendValue(keywords, String<char>("annual"));
	appendValue(keywords, String<char>("annually"));
	createTrie(trie,pos,keywords);
	// Dot Drawing
	String<String<char> > nodeMap;
	_createTrieNodeAttributes(trie, pos, nodeMap);
	String<String<char> > edgeMap;
	_createEdgeAttributes(trie,edgeMap);
	write(std::cout,trie,nodeMap, edgeMap, DotDrawing());

	// WordGraph
	typedef Graph<Automaton<Dna, String<Dna>, WordGraph<> > > TWordGraph;
	TWordGraph wordGr;
	addVertex(wordGr);addVertex(wordGr);addVertex(wordGr);addVertex(wordGr);
	assignRoot(wordGr,3);root(wordGr) = 2;addEdge(wordGr,0,3,"ag");
	addVertex(wordGr);addVertex(wordGr);addEdge(wordGr,0,5,"g");
	addVertex(wordGr);addVertex(wordGr);addEdge(wordGr,3,1,"aggg");
	addEdge(wordGr,3,4,"gg");addEdge(wordGr,5,2,"aggg");addEdge(wordGr,5,7,"g");
	addEdge(wordGr,7,6,"g");assignRoot(wordGr,0);
	write(std::cout,wordGr,DotDrawing());
}

	
//////////////////////////////////////////////////////////////////////////////


void Test_GraphUtils() {
	Test_GraphDrawing();

	debug::verifyCheckpoints("projects/library/seqan/graph_types/graph_drawing.h");
	debug::verifyCheckpoints("projects/library/seqan/graph_types/graph_utility_parsing.h");
}


}

#endif


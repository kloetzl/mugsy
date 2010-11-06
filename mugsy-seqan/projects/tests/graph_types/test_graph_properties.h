#ifndef SEQAN_HEADER_TEST_GRAPH_PROPERTIES_H
#define SEQAN_HEADER_TEST_GRAPH_PROPERTIES_H


namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////

template <typename TGraphType>
void Test_ExternalProperty() {
//____________________________________________________________________________
// Graph external property maps
	typedef Graph<TGraphType> TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename Size<TGraph>::Type TSize;
	
	TGraph g;
	TVertexDescriptor v0 = addVertex(g);
	TVertexDescriptor v1 = addVertex(g);
	addVertex(g);
	TEdgeDescriptor e2 =addEdge(g,0,2,'t');
	TEdgeDescriptor e1 =addEdge(g,v0,v1,'a');

	
	// Test external property maps
	String<int> dMap;
	resizeVertexMap(g,dMap);

	String<char> eMap;
	resizeEdgeMap(g,eMap);

	assignProperty(dMap, v0, 3);
	assignProperty(dMap, v1, 1);
	assignProperty(eMap, e1, 'a');
	assignProperty(eMap, e2, 'b');
	SEQAN_TASSERT(getProperty(dMap, v0) == 3)
	SEQAN_TASSERT(getProperty(dMap, v1) == 1)
	SEQAN_TASSERT(getProperty(eMap, e2) == 'b')
	SEQAN_TASSERT(getProperty(eMap, e1) == 'a')
	property(dMap, v1) = 2;
	property(eMap, e2) = 'c';
	SEQAN_TASSERT(getProperty(dMap, v1) == 2)
	SEQAN_TASSERT(getProperty(eMap, e2) == 'c')

	String<int> const dMap2(dMap);
	SEQAN_TASSERT(getProperty(dMap2, v0) == 3)
	SEQAN_TASSERT(getProperty(dMap2, v1) == 2)
	SEQAN_TASSERT(property(dMap2, v1) == 2)

	clear(g);
	addVertex(g);addVertex(g);addVertex(g);

	char names[] = {'r', 's','t'};
	String<char> nameMap;
	resizeVertexMap(g,nameMap, names);
	SEQAN_TASSERT(getProperty(nameMap, v0) == 'r')
	SEQAN_TASSERT(getProperty(nameMap, v1) == 's')
}


//////////////////////////////////////////////////////////////////////////////

void Test_Property() {
//____________________________________________________________________________
// Graph properties
	typedef Pair<char, int> TPair;
	typedef Directed<TPair> TEdges;
	typedef Graph<TEdges> TGraph;
	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;

	// Create a simple graph
	TGraph g;
	TVertexDescriptor v0 = addVertex(g);
	TVertexDescriptor v1 = addVertex(g);
	addVertex(g);
	TEdgeDescriptor e1 =addEdge(g,0,2);
	TEdgeDescriptor e2 =addEdge(g,v0,v1);

	// First Variant: Explicit internal map with member Ids
	InternalMap<TPair, 1> eMap1; // This property map is used to access the first member
	InternalMap<TPair, 2> eMap2; // This property map is used to access the second member
	resizeEdgeMap(g,eMap1);
	resizeEdgeMap(g,eMap2);
	assignProperty(eMap1, e1, 'a');
	assignProperty(eMap2, e1, 20);
	assignProperty(eMap1, e2, 'b');
	assignProperty(eMap2, e2, 50);
	SEQAN_TASSERT(getProperty(eMap1, e1) == 'a')
	SEQAN_TASSERT(getProperty(eMap2, e1) == 20)
	SEQAN_TASSERT(getProperty(eMap1, e2) == 'b')
	SEQAN_TASSERT(getProperty(eMap2, e2) == 50)
	// Note: That these properties are stored inside the cargo of each edge
	SEQAN_TASSERT(getCargo(e1).i1 == 'a')
	SEQAN_TASSERT(getCargo(e1).i2 == 20)
	SEQAN_TASSERT(getCargo(e2).i1 == 'b')
	SEQAN_TASSERT(getCargo(e2).i2 == 50)
	assignProperty(eMap1, e1, 'c');
	assignProperty(eMap2, e1, 10);
	SEQAN_TASSERT(getProperty(eMap1, e1) == 'c')
	SEQAN_TASSERT(getProperty(eMap2, e1) == 10)
	SEQAN_TASSERT(property(eMap1, e1) == 'c')
	SEQAN_TASSERT(property(eMap2, e1) == 10)
	InternalMap<TPair, 1> const eMap3(eMap1);
	InternalMap<TPair, 2> const eMap31(eMap2);
	SEQAN_TASSERT(getProperty(eMap3, e1) == 'c')
	SEQAN_TASSERT(getProperty(eMap3, e2) == 'b')
	SEQAN_TASSERT(getProperty(eMap31, e1) == 10)
	SEQAN_TASSERT(property(eMap31, e1) == 10)
	SEQAN_TASSERT(property(eMap3, e2) == 'b')
	// Create a simple graph with unsigned int cargo
	typedef EdgeDescriptor<Graph<Directed<unsigned int> > >::Type TEdgeDescriptor2;
	Graph<Directed<unsigned int> > g2;
	addVertex(g2);
	addVertex(g2);
	TEdgeDescriptor2 edge1 =addEdge(g2,v0,v0);
	addEdge(g2,0,1);
	InternalMap<unsigned int> edgeMap;
	resizeEdgeMap(g2,edgeMap);
	assignProperty(edgeMap, edge1 ,3);
	SEQAN_TASSERT(getProperty(edgeMap, edge1) == 3)
	SEQAN_TASSERT(property(edgeMap, edge1) == 3)
	InternalMap<unsigned int> const edgeMap2(edgeMap);
	SEQAN_TASSERT(getProperty(edgeMap2, edge1) == 3)
	SEQAN_TASSERT(property(edgeMap2, edge1) == 3)

	// Second Variant: Pointer to member using a class
	InternalPointerMap<char TPair:: *, &TPair::i1> eMap4;
	InternalPointerMap<int TPair:: *, &TPair::i2> eMap5;
	resizeEdgeMap(g,eMap4);
	resizeEdgeMap(g,eMap5);
	assignProperty(eMap4, e1, 'c');
	assignProperty(eMap5, e1, 10);
	assignProperty(eMap4, e2, 'd');
	assignProperty(eMap5, e2, 30);
	SEQAN_TASSERT(getProperty(eMap4, e1) == 'c')
	SEQAN_TASSERT(getProperty(eMap5, e1) == 10)
	SEQAN_TASSERT(getProperty(eMap4, e2) == 'd')
	SEQAN_TASSERT(getProperty(eMap5, e2) == 30)
	property(eMap4,e1)='z';
	property(eMap5,e1)=100;
	SEQAN_TASSERT(getProperty(eMap4, e1) == 'z')
	SEQAN_TASSERT(getProperty(eMap5, e1) == 100)
	InternalPointerMap<char TPair:: *, &TPair::i1> const eMap6(eMap4);
	SEQAN_TASSERT(getProperty(eMap6, e1) == 'z')
	SEQAN_TASSERT(getProperty(eMap6, e2) == 'd')
	SEQAN_TASSERT(property(eMap6, e2) == 'd')
	
	// Third Variant: Raw pointer to member
	char TPair:: * pseudo_map = &TPair::i1;
	assignProperty(pseudo_map, e1, 'z');
	assignProperty(pseudo_map, e2, 'w');
	SEQAN_TASSERT(getProperty(pseudo_map, e1) == 'z')
	SEQAN_TASSERT(getProperty(pseudo_map, e2) == 'w')
	property(pseudo_map,e1)='k';
	SEQAN_TASSERT(getProperty(pseudo_map, e1) == 'k')


	// Test shortcuts
	unsigned int weights[] = {4,8};
	Graph<Directed<void> > g10;
	addVertex(g10);addVertex(g10);addVertex(g10);
	addEdge(g10,0,1);addEdge(g10,0,2);
	String<int> weightMap;
	resizeEdgeMap(g10, weightMap, weights);
	SEQAN_TASSERT(getProperty(weightMap, findEdge(g10, 0, 1)) == 4)
	SEQAN_TASSERT(getProperty(weightMap, findEdge(g10, 0, 2)) == 8)
}

//////////////////////////////////////////////////////////////////////////////

void Test_GraphProperties() {
	Test_ExternalProperty<Directed<char> >();
	Test_ExternalProperty<Undirected<char> >();
	Test_ExternalProperty<Tree<char> >();
	Test_ExternalProperty<Automaton<char> >();	
	Test_Property();


	debug::verifyCheckpoints("projects/library/seqan/graph_types/graph_property.h");
}


}

#endif


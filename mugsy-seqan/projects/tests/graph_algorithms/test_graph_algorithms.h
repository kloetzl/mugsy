#ifndef SEQAN_HEADER_TEST_GRAPH_ALGORITHMS_H
#define SEQAN_HEADER_TEST_GRAPH_ALGORITHMS_H

namespace SEQAN_NAMESPACE_MAIN
{

template<typename TGraph>
inline void 
_createRandomGraph(TGraph& g) {
	clear(g);
	unsigned int nVertices = (mtRand() % 50) + 10;
	for(unsigned int i=0; i<nVertices; ++i) addVertex(g);
	unsigned int maxEdges = (mtRand() % 100) + 2 * nVertices;
	unsigned int nEdges = 0;
	for(unsigned int i=0; i<maxEdges; ++i) {
		unsigned int source = 0;
		unsigned int target = 0;
		do {
			source = (mtRand() % nVertices);
			target = (mtRand() % nVertices);
		} while (source == target);
		if (findEdge(g, source, target) == 0) {
			++nEdges;
			addEdge(g, source, target);
		}
	}
}


//////////////////////////////////////////////////////////////////////////////
void Test_HeapTree() {
	// Keyless heap
	String<int> test;
	appendValue(test, 4);appendValue(test, 1);appendValue(test, 3);
	appendValue(test, 2);appendValue(test, 16);appendValue(test, 9);
	appendValue(test, 10);appendValue(test, 14);
	appendValue(test, 8);appendValue(test, 7);

	HeapTree<int, std::greater<int> > firstHeap;
	SEQAN_TASSERT(empty(firstHeap) == true)
	SEQAN_TASSERT(length(firstHeap) == 0)
	buildHeap(firstHeap, begin(test), end(test));
	HeapTree<int, std::greater<int> > firstHeapTest(firstHeap);
	SEQAN_TASSERT(length(firstHeapTest) == length(firstHeap))
	clear(firstHeapTest);
	firstHeapTest = firstHeap;
	SEQAN_TASSERT(length(firstHeapTest) == length(firstHeap))
	SEQAN_TASSERT(empty(firstHeap) == false)
	SEQAN_TASSERT(length(firstHeap) == 10)
	clear(firstHeap);
	SEQAN_TASSERT(empty(firstHeap) == true)
	SEQAN_TASSERT(length(firstHeap) == 0)
	buildHeap(firstHeap, begin(test), end(test));
	SEQAN_TASSERT(empty(firstHeap) == false)
	SEQAN_TASSERT(length(firstHeap) == 10)
	SEQAN_TASSERT(heapRoot(firstHeap) == 16)
	SEQAN_TASSERT(heapExtractRoot(firstHeap) == 16)
	SEQAN_TASSERT(heapRoot(firstHeap) == 14)
	SEQAN_TASSERT(heapExtractRoot(firstHeap) == 14)
	SEQAN_TASSERT(heapRoot(firstHeap) == 10)
	SEQAN_TASSERT(heapExtractRoot(firstHeap) == 10)
	SEQAN_TASSERT(heapRoot(firstHeap) == 9)
	SEQAN_TASSERT(heapExtractRoot(firstHeap) == 9)
	SEQAN_TASSERT(heapRoot(firstHeap) == 8)
	SEQAN_TASSERT(heapExtractRoot(firstHeap) == 8)
	SEQAN_TASSERT(heapRoot(firstHeap) == 7)
	SEQAN_TASSERT(heapExtractRoot(firstHeap) == 7)
	SEQAN_TASSERT(heapRoot(firstHeap) == 4)
	SEQAN_TASSERT(heapExtractRoot(firstHeap) == 4)
	SEQAN_TASSERT(heapRoot(firstHeap) == 3)
	SEQAN_TASSERT(heapExtractRoot(firstHeap) == 3)
	SEQAN_TASSERT(heapRoot(firstHeap) == 2)
	SEQAN_TASSERT(heapExtractRoot(firstHeap) == 2)
	SEQAN_TASSERT(heapRoot(firstHeap) == 1)
	SEQAN_TASSERT(heapExtractRoot(firstHeap) == 1)
	SEQAN_TASSERT(empty(firstHeap) == true)
	SEQAN_TASSERT(length(firstHeap) == 0)


	// Keyed Heap
	typedef Pair<char, int> TKeyValue;
	typedef HeapTree<TKeyValue, std::greater<int>, KeyedHeap<> > TKeyedHeap;
	TKeyedHeap keyedHeap;
	String<TKeyValue> keyedTest;
	appendValue(keyedTest, TKeyValue('i', 4));appendValue(keyedTest, TKeyValue('a', 1));
	appendValue(keyedTest, TKeyValue('b', 3));appendValue(keyedTest, TKeyValue('c', 2));
	appendValue(keyedTest, TKeyValue('j', 16));appendValue(keyedTest, TKeyValue('m', 9));
	appendValue(keyedTest, TKeyValue('o', 10));appendValue(keyedTest, TKeyValue('p', 14));
	appendValue(keyedTest, TKeyValue('h', 8));appendValue(keyedTest, TKeyValue('z', 7));

	TKeyedHeap firstKeyedHeap;
	SEQAN_TASSERT(empty(firstKeyedHeap) == true)
	SEQAN_TASSERT(length(firstKeyedHeap) == 0)
	buildHeap(firstKeyedHeap, begin(keyedTest), end(keyedTest));
	TKeyedHeap firstKeyedHeapTest(firstKeyedHeap);
	SEQAN_TASSERT(length(firstKeyedHeapTest) == length(firstKeyedHeap))
	clear(firstKeyedHeapTest);
	firstKeyedHeapTest = firstKeyedHeap;
	SEQAN_TASSERT(length(firstKeyedHeapTest) == length(firstKeyedHeap))

	SEQAN_TASSERT(empty(firstKeyedHeap) == false)
	SEQAN_TASSERT(length(firstKeyedHeap) == 10)
	clear(firstKeyedHeap);
	SEQAN_TASSERT(empty(firstKeyedHeap) == true)
	SEQAN_TASSERT(length(firstKeyedHeap) == 0)
	buildHeap(firstKeyedHeap, begin(keyedTest), end(keyedTest));
	SEQAN_TASSERT(empty(firstKeyedHeap) == false)
	SEQAN_TASSERT(length(firstKeyedHeap) == 10)
	SEQAN_TASSERT(heapRoot(firstKeyedHeap).i1 == 'j')
	SEQAN_TASSERT(heapRoot(firstKeyedHeap).i2 == 16)
	SEQAN_TASSERT(heapExtractRoot(firstKeyedHeap).i2 == 16)
	SEQAN_TASSERT(heapRoot(firstKeyedHeap).i1 == 'p')
	SEQAN_TASSERT(heapRoot(firstKeyedHeap).i2 == 14)
	SEQAN_TASSERT(heapExtractRoot(firstKeyedHeap).i2 == 14)
	SEQAN_TASSERT(heapRoot(firstKeyedHeap).i1 == 'o')
	SEQAN_TASSERT(heapRoot(firstKeyedHeap).i2 == 10)
	SEQAN_TASSERT(heapExtractRoot(firstKeyedHeap).i2 == 10)
	SEQAN_TASSERT(heapRoot(firstKeyedHeap).i2 == 9)
	SEQAN_TASSERT(heapExtractRoot(firstKeyedHeap).i2 == 9)
	SEQAN_TASSERT(heapRoot(firstKeyedHeap).i2 == 8)
	SEQAN_TASSERT(heapExtractRoot(firstKeyedHeap).i2 == 8)
	SEQAN_TASSERT(heapRoot(firstKeyedHeap).i2 == 7)
	SEQAN_TASSERT(heapExtractRoot(firstKeyedHeap).i2 == 7)
	SEQAN_TASSERT(heapRoot(firstKeyedHeap).i2 == 4)
	SEQAN_TASSERT(heapExtractRoot(firstKeyedHeap).i2 == 4)
	SEQAN_TASSERT(heapRoot(firstKeyedHeap).i2 == 3)
	SEQAN_TASSERT(heapExtractRoot(firstKeyedHeap).i2 == 3)
	SEQAN_TASSERT(heapRoot(firstKeyedHeap).i2 == 2)
	SEQAN_TASSERT(heapExtractRoot(firstKeyedHeap).i2 == 2)
	SEQAN_TASSERT(heapRoot(firstKeyedHeap).i2 == 1)
	SEQAN_TASSERT(heapExtractRoot(firstKeyedHeap).i2 == 1)
	SEQAN_TASSERT(empty(firstKeyedHeap) == true)
	SEQAN_TASSERT(length(firstKeyedHeap) == 0)
	
	// Change value only for keyed heaps
	buildHeap(firstKeyedHeap, begin(keyedTest), end(keyedTest));
	heapChangeValue(firstKeyedHeap, 'c', 20);
	SEQAN_TASSERT(heapRoot(firstKeyedHeap).i1 == 'c')
	SEQAN_TASSERT(heapRoot(firstKeyedHeap).i2 == 20)
	heapChangeValue(firstKeyedHeap, 'c', 2);
	SEQAN_TASSERT(heapRoot(firstKeyedHeap).i1 == 'j')
	SEQAN_TASSERT(heapRoot(firstKeyedHeap).i2 == 16)
	heapInsert(firstKeyedHeap, TKeyValue('x', 100));
	SEQAN_TASSERT(heapRoot(firstKeyedHeap).i1 == 'x')
	SEQAN_TASSERT(heapRoot(firstKeyedHeap).i2 == 100)
	heapChangeValue(firstKeyedHeap, 'c', 200);
	SEQAN_TASSERT(heapRoot(firstKeyedHeap).i1 == 'c')
	SEQAN_TASSERT(heapExtractRoot(firstKeyedHeap).i2 == 200)

	// Test Heap Sort
	typedef int TValue;
	String<TValue> result1;
	String<TValue> result2;
	
	for(unsigned int i=0; i<1000; ++i) {
		TValue val = (mtRand() % 10000) - 5000;
		appendValue(result1, val);
		appendValue(result2, val);
	}
	std::sort(begin(result1, Standard()), end(result1, Standard()));
	heapSort(begin(result2, Standard()), end(result2, Standard()));
	if (result1 != result2) {
		for(unsigned int i=0; i<length(result1);++i) std::cout << value(result1, i) << ',';
		std::cout << std::endl;
		for(unsigned int i=0; i<length(result2);++i) std::cout << value(result2, i) << ',';
		std::cout << std::endl;
		std::cout << "Error" << std::endl;
		exit(0);
	}

	// Test Heap Sort
	typedef int TValue;
	clear(result1);
	clear(result2);
	
	for(unsigned int i=0; i<1000; ++i) {
		TValue val = (mtRand() % 10000) - 5000;
		appendValue(result1, val);
		appendValue(result2, val);
	}
	std::sort(begin(result1, Standard()), end(result1, Standard()), std::greater<TValue>());
	heapSort(begin(result2, Standard()), end(result2, Standard()), std::greater<TValue>());
	if (result1 != result2) {
		for(unsigned int i=0; i<length(result1);++i) std::cout << value(result1, i) << ',';
		std::cout << std::endl;
		for(unsigned int i=0; i<length(result2);++i) std::cout << value(result2, i) << ',';
		std::cout << std::endl;
		std::cout << "Error" << std::endl;
		exit(0);
	}
}

//////////////////////////////////////////////////////////////////////////////
void Test_BreadthFirstSearch() {
//____________________________________________________________________________
// Breadth-First Search
	typedef Graph<Undirected<> > TGraph;
	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef Size<TGraph>::Type TSize;

	//Number of edges
	TSize numEdges = 10;
	//Source, Target, Source, Target, Source, ...
	TVertexDescriptor edges[] = {0,1, 0,4, 1,5, 2,5, 2,6, 2,3, 3,6, 3,7, 5,6, 6,7};
	
	//Create the graph
	TGraph g;
	addEdges(g,edges, numEdges);

	// Predecessor and distance map
	String<unsigned int> predMap;
	String<unsigned int> distMap;

	// Bfs
	breadth_first_search(g, 1, predMap, distMap);
	
	SEQAN_TASSERT(getProperty(distMap, 0) == 1)
	SEQAN_TASSERT(getProperty(distMap, 1) == 0)
	SEQAN_TASSERT(getProperty(distMap, 2) == 2)
	SEQAN_TASSERT(getProperty(distMap, 3) == 3)
	SEQAN_TASSERT(getProperty(distMap, 4) == 2)
	SEQAN_TASSERT(getProperty(distMap, 5) == 1)
	SEQAN_TASSERT(getProperty(distMap, 6) == 2)
	SEQAN_TASSERT(getProperty(distMap, 7) == 3)
	SEQAN_TASSERT(getProperty(predMap, 0) == 1)
	SEQAN_TASSERT(getProperty(predMap, 1) == getNil<TVertexDescriptor>())
	SEQAN_TASSERT(getProperty(predMap, 2) == 5)
	SEQAN_TASSERT(getProperty(predMap, 3) == 6)
	SEQAN_TASSERT(getProperty(predMap, 4) == 0)
	SEQAN_TASSERT(getProperty(predMap, 5) == 1)
	SEQAN_TASSERT(getProperty(predMap, 6) == 5)
	SEQAN_TASSERT(getProperty(predMap, 7) == 6)
}

//////////////////////////////////////////////////////////////////////////////

void Test_DepthFirstSearch() {
//____________________________________________________________________________
// Depth-First Search
	typedef VertexDescriptor<Graph<> >::Type TVertexDescriptor;
	typedef EdgeDescriptor<Graph<> >::Type TEdgeDescriptor;
	typedef Size<Graph<> >::Type TSize;

	//Number of edges
	TSize numEdges = 8;
	//Source, Target, Source, Target, Source, ...
	TVertexDescriptor edges[] = {0,3, 0,1, 1,4, 2,4, 2,5, 3,1, 4,3, 5,5};

	//Create the graph
	Graph<> g;
	addEdges(g,edges, numEdges);

	// Predecessor and distance map
	String<unsigned int> predMap;
	String<unsigned int> discoveryTimeMap;
	String<unsigned int> finishingTimeMap;

	// Dfs
	depth_first_search(g, predMap, discoveryTimeMap, finishingTimeMap);

	SEQAN_TASSERT(getProperty(discoveryTimeMap, 0) == 1)
	SEQAN_TASSERT(getProperty(discoveryTimeMap, 1) == 2)
	SEQAN_TASSERT(getProperty(discoveryTimeMap, 2) == 9)
	SEQAN_TASSERT(getProperty(discoveryTimeMap, 3) == 4)
	SEQAN_TASSERT(getProperty(discoveryTimeMap, 4) == 3)
	SEQAN_TASSERT(getProperty(discoveryTimeMap, 5) == 10)
	SEQAN_TASSERT(getProperty(finishingTimeMap, 0) == 8)
	SEQAN_TASSERT(getProperty(finishingTimeMap, 1) == 7)
	SEQAN_TASSERT(getProperty(finishingTimeMap, 2) == 12)
	SEQAN_TASSERT(getProperty(finishingTimeMap, 3) == 5)
	SEQAN_TASSERT(getProperty(finishingTimeMap, 4) == 6)
	SEQAN_TASSERT(getProperty(finishingTimeMap, 5) == 11)
	SEQAN_TASSERT(getProperty(predMap, 0) == getNil<TVertexDescriptor>())
	SEQAN_TASSERT(getProperty(predMap, 1) == 0)
	SEQAN_TASSERT(getProperty(predMap, 2) == getNil<TVertexDescriptor>())
	SEQAN_TASSERT(getProperty(predMap, 3) == 4)
	SEQAN_TASSERT(getProperty(predMap, 4) == 1)
	SEQAN_TASSERT(getProperty(predMap, 5) == 2)
}

//////////////////////////////////////////////////////////////////////////////

void Test_TopologicalSort() {
//____________________________________________________________________________
// Topological Sort
	typedef VertexDescriptor<Graph<> >::Type TVertexDescriptor;
	typedef EdgeDescriptor<Graph<> >::Type TEdgeDescriptor;
	typedef Size<Graph<> >::Type TSize;

	//Number of edges
	TSize numEdges = 9;
	//Source, Target, Source, Target, Source, ...
	TVertexDescriptor edges[] = {0,3, 0,1, 1,2, 3,2, 5,7, 5,6, 6,7, 6,3, 8,7};
	
	//Create the graph
	Graph<> g;
	addEdges(g,edges, numEdges);

	// Predecessor and distance map
	String<TVertexDescriptor> order;
	
	// Topological sort
	topological_sort(g, order);

	SEQAN_TASSERT(getValue(order, 0) == 8)
	SEQAN_TASSERT(getValue(order, 1) == 5)
	SEQAN_TASSERT(getValue(order, 2) == 6)
	SEQAN_TASSERT(getValue(order, 3) == 7)
	SEQAN_TASSERT(getValue(order, 4) == 4)
	SEQAN_TASSERT(getValue(order, 5) == 0)
	SEQAN_TASSERT(getValue(order, 6) == 3)
	SEQAN_TASSERT(getValue(order, 7) == 1)
	SEQAN_TASSERT(getValue(order, 8) == 2)
}

//////////////////////////////////////////////////////////////////////////////

void Test_StronglyConnectedComponents() {
//____________________________________________________________________________
// Strongly-Connected-Components

	typedef VertexDescriptor<Graph<> >::Type TVertexDescriptor;
	typedef EdgeDescriptor<Graph<> >::Type TEdgeDescriptor;
	typedef Size<Graph<> >::Type TSize;

	//Number of edges
	TSize numEdges = 14;
	//Source, Target, Source, Target, Source, ...
	TVertexDescriptor edges[] = {1,0, 0,4, 2,1, 4,1, 5,1, 6,2, 3,2, 2,3, 7,3, 5,4, 6,5, 5,6, 7,6, 7,7};

	//Create the graph
	Graph<> g;
	addEdges(g,edges, numEdges);

	// Predecessor and distance map
	String<unsigned int> component;
	
	// Strongly Connected Components
	strongly_connected_components(g, component);

	SEQAN_TASSERT(getValue(component, 0) == 3)
	SEQAN_TASSERT(getValue(component, 1) == 3)
	SEQAN_TASSERT(getValue(component, 2) == 2)
	SEQAN_TASSERT(getValue(component, 3) == 2)
	SEQAN_TASSERT(getValue(component, 4) == 3)
	SEQAN_TASSERT(getValue(component, 5) == 1)
	SEQAN_TASSERT(getValue(component, 6) == 1)
	SEQAN_TASSERT(getValue(component, 7) == 0)
}

//////////////////////////////////////////////////////////////////////////////

void Test_ConnectedComponents() {
//____________________________________________________________________________
// Connected-Components

	typedef Graph<Undirected<> > TGraph;
	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef Size<TGraph>::Type TSize;

	//Number of edges
	TSize numEdges = 3;
	//Source, Target, Source, Target, Source, ...
	TVertexDescriptor edges[] = {0,1, 2,3, 3,4};

	//Create the graph
	TGraph g;
	addEdges(g,edges, numEdges);

	//Components
	String<unsigned int> component;
	
	//Connected Components
	connected_components(g, component);

	SEQAN_TASSERT(getValue(component, 0) == 0)
	SEQAN_TASSERT(getValue(component, 1) == 0)
	SEQAN_TASSERT(getValue(component, 2) == 1)
	SEQAN_TASSERT(getValue(component, 3) == 1)
	SEQAN_TASSERT(getValue(component, 4) == 1)
}


//////////////////////////////////////////////////////////////////////////////

void Test_PrimsAlgorithm() {
//____________________________________________________________________________
// Prim's algorithm
	typedef Graph<Undirected<> > TGraph;
	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef Size<TGraph>::Type TSize;

	//Number of edges
	TSize numEdges = 14;
	//Source, Target, Source, Target, Source, ...
	TVertexDescriptor edges[] = {0,1, 0,6, 1,2, 1,6, 2,3, 2,4, 2,8, 3,5, 3,8, 4,6, 4,7, 5,8, 6,7, 7,8};
	unsigned int weights[] =    {4,   8,   8,   11,  7,   2,   4,   9,   14,  7,   6,   10,  1,   2  };

	//Create the graph
	TGraph g;
	addEdges(g,edges, numEdges);
	String<unsigned int> weightMap;
	resizeEdgeMap(g, weightMap, weights);

	// Tree and predecessor map 
	String<TVertexDescriptor> predMap;

	prims_algorithm(g, 0, weightMap, predMap);

	SEQAN_TASSERT(getProperty(predMap, 0) == getNil<TVertexDescriptor>())
	SEQAN_TASSERT(getProperty(predMap, 1) == 0)
	SEQAN_TASSERT(getProperty(predMap, 2) == 1)
	SEQAN_TASSERT(getProperty(predMap, 3) == 2)
	SEQAN_TASSERT(getProperty(predMap, 4) == 2)
	SEQAN_TASSERT(getProperty(predMap, 5) == 3)
	SEQAN_TASSERT(getProperty(predMap, 6) == 7)
	SEQAN_TASSERT(getProperty(predMap, 7) == 8)
	SEQAN_TASSERT(getProperty(predMap, 8) == 2)
}


//////////////////////////////////////////////////////////////////////////////

void Test_KruskalsAlgorithm() {
//____________________________________________________________________________
// Kruskal's algorithm
	typedef Graph<Undirected<> > TGraph;
	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef Size<TGraph>::Type TSize;

	//Number of edges
	TSize numEdges = 14;
	//Source, Target, Source, Target, Source, ...
	TVertexDescriptor edges[] = {0,1, 0,6, 1,2, 1,6, 2,3, 2,4, 2,8, 3,5, 3,8, 4,6, 4,7, 5,8, 6,7, 7,8};
	unsigned int weights[] =    {4,   8,   8,   11,  7,   2,   4,   9,   14,  7,   6,   10,  1,   2  };

	//Create the graph
	TGraph g;
	addEdges(g,edges, numEdges);
	String<unsigned int> weightMap;
	resizeEdgeMap(g, weightMap, weights);

	// Tree edges
	String<TVertexDescriptor> treeEdges;
	kruskals_algorithm(g, 0, weightMap, treeEdges);

	SEQAN_TASSERT(getValue(treeEdges, 0) == 6)
	SEQAN_TASSERT(getValue(treeEdges, 1) == 7)
	SEQAN_TASSERT(getValue(treeEdges, 2) == 2)
	SEQAN_TASSERT(getValue(treeEdges, 3) == 4)
	SEQAN_TASSERT(getValue(treeEdges, 4) == 7)
	SEQAN_TASSERT(getValue(treeEdges, 5) == 8)
	SEQAN_TASSERT(getValue(treeEdges, 6) == 0)
	SEQAN_TASSERT(getValue(treeEdges, 7) == 1)
	SEQAN_TASSERT(getValue(treeEdges, 8) == 2)
	SEQAN_TASSERT(getValue(treeEdges, 9) == 8)
	SEQAN_TASSERT(getValue(treeEdges, 10) == 2)
	SEQAN_TASSERT(getValue(treeEdges, 11) == 3)
	SEQAN_TASSERT(getValue(treeEdges, 12) == 0)
	SEQAN_TASSERT(getValue(treeEdges, 13) == 6)
	SEQAN_TASSERT(getValue(treeEdges, 14) == 3)
	SEQAN_TASSERT(getValue(treeEdges, 15) == 5)
}

//////////////////////////////////////////////////////////////////////////////

void Test_MST_All() {
	//while (true) {
		Graph<Undirected<> > myGraph;
		_createRandomGraph(myGraph);
		String<unsigned int> initialWeights;
		for(unsigned int i = 0; i<numEdges(myGraph); ++i) appendValue(initialWeights, (mtRand() % 1000));
		String<unsigned int> weightMapInput;
		resizeEdgeMap(myGraph, weightMapInput, initialWeights);
		clear(initialWeights);
	
		// Prim1
		String<unsigned int> predMapOut;
		prims_algorithm_spaceEfficient(myGraph, 0, weightMapInput, predMapOut);
		typedef std::set<EdgeDescriptor<Graph<Undirected<> > >::Type> TEdgeSet;
		TEdgeSet edgeSet1;
		_collect_Edges(myGraph,predMapOut,0,edgeSet1);
		unsigned int sum1 = 0;
		for(TEdgeSet::const_iterator pos = edgeSet1.begin(); pos != edgeSet1.end(); ++pos) {
			sum1 += getProperty(weightMapInput, *pos);
		}

		// Prim2
		String<unsigned int> predMapOut2;
		prims_algorithm(myGraph, 0, weightMapInput, predMapOut2);
		TEdgeSet edgeSet2;
		_collect_Edges(myGraph,predMapOut2,0,edgeSet2);
		unsigned int sum2 = 0;
		for(TEdgeSet::const_iterator pos = edgeSet2.begin(); pos != edgeSet2.end(); ++pos) {
			sum2 += getProperty(weightMapInput, *pos);
		}

		unsigned int sum3 = 0;
		if (sum1 != 0) {
			// Kruskal
			String<unsigned int> treeEdges;
			kruskals_algorithm(myGraph, 0, weightMapInput, treeEdges);
			for(unsigned int i = 0; i < length(treeEdges); i = i + 2) {
				sum3 += getProperty(weightMapInput, findEdge(myGraph, value(treeEdges, i), value(treeEdges, i + 1)));
			}
		}

		//std::cout << sum1 << ',';
		SEQAN_TASSERT(sum1 == sum2)
		SEQAN_TASSERT(sum2 == sum3)
	//}

}

//////////////////////////////////////////////////////////////////////////////

void Test_DagShortestPath() {
//____________________________________________________________________________
// DAG-Shortest Paths
	typedef VertexDescriptor<Graph<> >::Type TVertexDescriptor;
	typedef EdgeDescriptor<Graph<> >::Type TEdgeDescriptor;
	typedef Size<Graph<> >::Type TSize;

	//Number of edges
	TSize numEdges = 10;
	//Source, Target, Source, Target, Source, ...
	TVertexDescriptor edges[] = {0,2, 0,1, 1,3, 1,2, 2,5, 2,4, 2,3, 3,5, 3,4, 4,5};
	int weights[] =             {3,   5,   6,   2,   2,   4,   7,   1,   -1,  -2};
	
	//Create the graph
	Graph<> g;
	addEdges(g,edges, numEdges);

	String<int> weightMap;
	resizeEdgeMap(g, weightMap, weights);

	// Predecessor map and distance map
	String<unsigned int> predMap;
	String<unsigned int> distMap;

	// DAG-Shortest path(Graph, sourceVertex_vertex, weightMap, predMap, distMap)
	dag_shortest_path(g,1,weightMap,predMap,distMap);
	
	SEQAN_TASSERT(getProperty(distMap, 1) == 0)
	SEQAN_TASSERT(getProperty(distMap, 2) == 2)
	SEQAN_TASSERT(getProperty(distMap, 3) == 6)
	SEQAN_TASSERT(getProperty(distMap, 4) == 5)
	SEQAN_TASSERT(getProperty(distMap, 5) == 3)
	SEQAN_TASSERT(getProperty(predMap, 0) == getNil<TVertexDescriptor>())
	SEQAN_TASSERT(getProperty(predMap, 1) == getNil<TVertexDescriptor>())
	SEQAN_TASSERT(getProperty(predMap, 2) == 1)
	SEQAN_TASSERT(getProperty(predMap, 3) == 1)
	SEQAN_TASSERT(getProperty(predMap, 4) == 3)
	SEQAN_TASSERT(getProperty(predMap, 5) == 4)
}

//////////////////////////////////////////////////////////////////////////////

void Test_BellmanFord() {
//____________________________________________________________________________
// Bellman-Ford

	typedef VertexDescriptor<Graph<> >::Type TVertexDescriptor;
	typedef EdgeDescriptor<Graph<> >::Type TEdgeDescriptor;
	typedef Size<Graph<> >::Type TSize;

	//Number of edges
	TSize numEdges = 10;
	//Source, Target, Source, Target, Source, ...
	TVertexDescriptor edges[] = {0,1, 0,3, 1,2, 1,3, 2,4, 3,1, 3,2, 3,4, 4,0, 4,2};
	unsigned int weights[] =    {10,  5,   1,   2,   4,   3,   9,   2,   7,   6};

	//Create the graph
	Graph<> g;
	addEdges(g,edges, numEdges);
	
	String<unsigned int> weightMap;
	resizeEdgeMap(g,weightMap, weights);

	// Out parameters of Bellman-Ford: Predecessor map and distance map
	String<unsigned int> predMap;
	String<unsigned int> distMap;

	// Bellman-Ford
	bool noNegativeCycle = bellman_ford_algorithm(g,0,weightMap,predMap,distMap);

	SEQAN_TASSERT(getProperty(distMap, 0) == 0)
	SEQAN_TASSERT(getProperty(distMap, 1) == 8)
	SEQAN_TASSERT(getProperty(distMap, 2) == 9)
	SEQAN_TASSERT(getProperty(distMap, 3) == 5)
	SEQAN_TASSERT(getProperty(distMap, 4) == 7)
	SEQAN_TASSERT(getProperty(predMap, 0) == getNil<TVertexDescriptor>())
	SEQAN_TASSERT(getProperty(predMap, 1) == 3)
	SEQAN_TASSERT(getProperty(predMap, 2) == 1)
	SEQAN_TASSERT(getProperty(predMap, 3) == 0)
	SEQAN_TASSERT(getProperty(predMap, 4) == 3)
	SEQAN_TASSERT(noNegativeCycle == true)
}

//////////////////////////////////////////////////////////////////////////////

void Test_Dijkstra() {
//____________________________________________________________________________
// Dijkstra

	typedef VertexDescriptor<Graph<> >::Type TVertexDescriptor;
	typedef EdgeDescriptor<Graph<> >::Type TEdgeDescriptor;
	typedef Size<Graph<> >::Type TSize;

	//Number of edges
	TSize numEdges = 10;
	//Source, Target, Source, Target, Source, ...
	TVertexDescriptor edges[] = {0,1, 0,3, 1,2, 1,3, 2,4, 3,1, 3,2, 3,4, 4,0, 4,2};
	unsigned int weights[] =    {10,  5,   1,   2,   4,   3,   9,   2,   7,   6};

	//Create the graph
	Graph<> g;
	addEdges(g,edges, numEdges);

	String<unsigned int> weightMap;
	resizeEdgeMap(g , weightMap, weights);

	// Out parameters of Dijkstra: Predecessor map and distance map
	String<unsigned int> predMap;
	String<unsigned int> distMap;

	// Dijkstra
	dijkstra(g,0,weightMap,predMap,distMap);

	SEQAN_TASSERT(getProperty(distMap, 0) == 0)
	SEQAN_TASSERT(getProperty(distMap, 1) == 8)
	SEQAN_TASSERT(getProperty(distMap, 2) == 9)
	SEQAN_TASSERT(getProperty(distMap, 3) == 5)
	SEQAN_TASSERT(getProperty(distMap, 4) == 7)
	SEQAN_TASSERT(getProperty(predMap, 0) == getNil<TVertexDescriptor>())
	SEQAN_TASSERT(getProperty(predMap, 1) == 3)
	SEQAN_TASSERT(getProperty(predMap, 2) == 1)
	SEQAN_TASSERT(getProperty(predMap, 3) == 0)
	SEQAN_TASSERT(getProperty(predMap, 4) == 3)
}

//////////////////////////////////////////////////////////////////////////////

void Test_AllPairsShortestPath() {
//____________________________________________________________________________
// All-Pairs Shortest Path

	typedef VertexDescriptor<Graph<> >::Type TVertexDescriptor;
	typedef EdgeDescriptor<Graph<> >::Type TEdgeDescriptor;
	typedef Size<Graph<> >::Type TSize;

	//Number of edges
	TSize numEdges = 9;
	//Source, Target, Source, Target, Source, ...
	TVertexDescriptor edges[] = {0,1, 0,2, 0,4, 1,3, 1,4, 2,1, 3,0, 3,2, 4,3};
	int weights[] =    {3,   8,   -4,  1,   7,   4,   2,   -5,  6};

	//Create the graph
	Graph<> g;
	addEdges(g,edges, numEdges);
	
	String<int> weightMap;	
	resizeEdgeMap(g,weightMap, weights);

	// Out parameter
	String<int> distMat;
	String<TVertexDescriptor> predMat;

	// All-Pairs shortest path
	all_pairs_shortest_path(g,weightMap, distMat, predMat);

	unsigned int len = (unsigned int) sqrt((double) length(distMat));
	SEQAN_TASSERT(getValue(distMat, 0*len + 0) == 0)
	SEQAN_TASSERT(getValue(distMat, 0*len + 1) == 1)
	SEQAN_TASSERT(getValue(distMat, 0*len + 2) == -3)
	SEQAN_TASSERT(getValue(distMat, 0*len + 3) == 2)
	SEQAN_TASSERT(getValue(distMat, 0*len + 4) == -4)
	SEQAN_TASSERT(getValue(distMat, 1*len + 0) == 3)
	SEQAN_TASSERT(getValue(distMat, 1*len + 1) == 0)
	SEQAN_TASSERT(getValue(distMat, 1*len + 2) == -4)
	SEQAN_TASSERT(getValue(distMat, 1*len + 3) == 1)
	SEQAN_TASSERT(getValue(distMat, 1*len + 4) == -1)
	SEQAN_TASSERT(getValue(distMat, 2*len + 0) == 7)
	SEQAN_TASSERT(getValue(distMat, 2*len + 1) == 4)
	SEQAN_TASSERT(getValue(distMat, 2*len + 2) == 0)
	SEQAN_TASSERT(getValue(distMat, 2*len + 3) == 5)
	SEQAN_TASSERT(getValue(distMat, 2*len + 4) == 3)
	SEQAN_TASSERT(getValue(distMat, 3*len + 0) == 2)
	SEQAN_TASSERT(getValue(distMat, 3*len + 1) == -1)
	SEQAN_TASSERT(getValue(distMat, 3*len + 2) == -5)
	SEQAN_TASSERT(getValue(distMat, 3*len + 3) == 0)
	SEQAN_TASSERT(getValue(distMat, 3*len + 4) == -2)
	SEQAN_TASSERT(getValue(distMat, 4*len + 0) == 8)
	SEQAN_TASSERT(getValue(distMat, 4*len + 1) == 5)
	SEQAN_TASSERT(getValue(distMat, 4*len + 2) == 1)
	SEQAN_TASSERT(getValue(distMat, 4*len + 3) == 6)
	SEQAN_TASSERT(getValue(distMat, 4*len + 4) == 0)
	SEQAN_TASSERT(getValue(predMat, 0*len + 0) == getNil<TVertexDescriptor>())
	SEQAN_TASSERT(getValue(predMat, 0*len + 1) == 2)
	SEQAN_TASSERT(getValue(predMat, 0*len + 2) == 3)
	SEQAN_TASSERT(getValue(predMat, 0*len + 3) == 4)
	SEQAN_TASSERT(getValue(predMat, 0*len + 4) == 0)
	SEQAN_TASSERT(getValue(predMat, 1*len + 0) == 3)
	SEQAN_TASSERT(getValue(predMat, 1*len + 1) == getNil<TVertexDescriptor>())
	SEQAN_TASSERT(getValue(predMat, 1*len + 2) == 3)
	SEQAN_TASSERT(getValue(predMat, 1*len + 3) == 1)
	SEQAN_TASSERT(getValue(predMat, 1*len + 4) == 0)
	SEQAN_TASSERT(getValue(predMat, 2*len + 0) == 3)
	SEQAN_TASSERT(getValue(predMat, 2*len + 1) == 2)
	SEQAN_TASSERT(getValue(predMat, 2*len + 2) == getNil<TVertexDescriptor>())
	SEQAN_TASSERT(getValue(predMat, 2*len + 3) == 1)
	SEQAN_TASSERT(getValue(predMat, 2*len + 4) == 0)
	SEQAN_TASSERT(getValue(predMat, 3*len + 0) == 3)
	SEQAN_TASSERT(getValue(predMat, 3*len + 1) == 2)
	SEQAN_TASSERT(getValue(predMat, 3*len + 2) == 3)
	SEQAN_TASSERT(getValue(predMat, 3*len + 3) == getNil<TVertexDescriptor>())
	SEQAN_TASSERT(getValue(predMat, 3*len + 4) == 0)
	SEQAN_TASSERT(getValue(predMat, 4*len + 0) == 3)
	SEQAN_TASSERT(getValue(predMat, 4*len + 1) == 2)
	SEQAN_TASSERT(getValue(predMat, 4*len + 2) == 3)
	SEQAN_TASSERT(getValue(predMat, 4*len + 3) == 4)
	SEQAN_TASSERT(getValue(predMat, 4*len + 4) == getNil<TVertexDescriptor>())
}

//////////////////////////////////////////////////////////////////////////////

void Test_FloydWarshall() {
//____________________________________________________________________________
// Floyd Warshall

	typedef VertexDescriptor<Graph<> >::Type TVertexDescriptor;
	typedef EdgeDescriptor<Graph<> >::Type TEdgeDescriptor;
	typedef Size<Graph<> >::Type TSize;

	//Number of edges
	TSize numEdges = 9;
	//Source, Target, Source, Target, Source, ...
	TVertexDescriptor edges[] = {0,1, 0,2, 0,4, 1,3, 1,4, 2,1, 3,0, 3,2, 4,3};
	int weights[] =    {3,   8,   -4,  1,   7,   4,   2,   -5,  6};

	//Create the graph
	Graph<> g;
	addEdges(g,edges, numEdges);
	
	String<int> weightMap;	
	resizeEdgeMap(g,weightMap, weights);

	// Out parameter
	String<int> distMat;
	String<TVertexDescriptor> predMat;

	// Floyd-Warshall
	floyd_warshall(g,weightMap, distMat, predMat);

	unsigned int len = (unsigned int) sqrt((double) length(distMat));
	SEQAN_TASSERT(getValue(distMat, 0*len + 0) == 0)
	SEQAN_TASSERT(getValue(distMat, 0*len + 1) == 1)
	SEQAN_TASSERT(getValue(distMat, 0*len + 2) == -3)
	SEQAN_TASSERT(getValue(distMat, 0*len + 3) == 2)
	SEQAN_TASSERT(getValue(distMat, 0*len + 4) == -4)
	SEQAN_TASSERT(getValue(distMat, 1*len + 0) == 3)
	SEQAN_TASSERT(getValue(distMat, 1*len + 1) == 0)
	SEQAN_TASSERT(getValue(distMat, 1*len + 2) == -4)
	SEQAN_TASSERT(getValue(distMat, 1*len + 3) == 1)
	SEQAN_TASSERT(getValue(distMat, 1*len + 4) == -1)
	SEQAN_TASSERT(getValue(distMat, 2*len + 0) == 7)
	SEQAN_TASSERT(getValue(distMat, 2*len + 1) == 4)
	SEQAN_TASSERT(getValue(distMat, 2*len + 2) == 0)
	SEQAN_TASSERT(getValue(distMat, 2*len + 3) == 5)
	SEQAN_TASSERT(getValue(distMat, 2*len + 4) == 3)
	SEQAN_TASSERT(getValue(distMat, 3*len + 0) == 2)
	SEQAN_TASSERT(getValue(distMat, 3*len + 1) == -1)
	SEQAN_TASSERT(getValue(distMat, 3*len + 2) == -5)
	SEQAN_TASSERT(getValue(distMat, 3*len + 3) == 0)
	SEQAN_TASSERT(getValue(distMat, 3*len + 4) == -2)
	SEQAN_TASSERT(getValue(distMat, 4*len + 0) == 8)
	SEQAN_TASSERT(getValue(distMat, 4*len + 1) == 5)
	SEQAN_TASSERT(getValue(distMat, 4*len + 2) == 1)
	SEQAN_TASSERT(getValue(distMat, 4*len + 3) == 6)
	SEQAN_TASSERT(getValue(distMat, 4*len + 4) == 0)
	SEQAN_TASSERT(getValue(predMat, 0*len + 0) == getNil<TVertexDescriptor>())
	SEQAN_TASSERT(getValue(predMat, 0*len + 1) == 2)
	SEQAN_TASSERT(getValue(predMat, 0*len + 2) == 3)
	SEQAN_TASSERT(getValue(predMat, 0*len + 3) == 4)
	SEQAN_TASSERT(getValue(predMat, 0*len + 4) == 0)
	SEQAN_TASSERT(getValue(predMat, 1*len + 0) == 3)
	SEQAN_TASSERT(getValue(predMat, 1*len + 1) == getNil<TVertexDescriptor>())
	SEQAN_TASSERT(getValue(predMat, 1*len + 2) == 3)
	SEQAN_TASSERT(getValue(predMat, 1*len + 3) == 1)
	SEQAN_TASSERT(getValue(predMat, 1*len + 4) == 0)
	SEQAN_TASSERT(getValue(predMat, 2*len + 0) == 3)
	SEQAN_TASSERT(getValue(predMat, 2*len + 1) == 2)
	SEQAN_TASSERT(getValue(predMat, 2*len + 2) == getNil<TVertexDescriptor>())
	SEQAN_TASSERT(getValue(predMat, 2*len + 3) == 1)
	SEQAN_TASSERT(getValue(predMat, 2*len + 4) == 0)
	SEQAN_TASSERT(getValue(predMat, 3*len + 0) == 3)
	SEQAN_TASSERT(getValue(predMat, 3*len + 1) == 2)
	SEQAN_TASSERT(getValue(predMat, 3*len + 2) == 3)
	SEQAN_TASSERT(getValue(predMat, 3*len + 3) == getNil<TVertexDescriptor>())
	SEQAN_TASSERT(getValue(predMat, 3*len + 4) == 0)
	SEQAN_TASSERT(getValue(predMat, 4*len + 0) == 3)
	SEQAN_TASSERT(getValue(predMat, 4*len + 1) == 2)
	SEQAN_TASSERT(getValue(predMat, 4*len + 2) == 3)
	SEQAN_TASSERT(getValue(predMat, 4*len + 3) == 4)
	SEQAN_TASSERT(getValue(predMat, 4*len + 4) == getNil<TVertexDescriptor>())
}


//////////////////////////////////////////////////////////////////////////////

void Test_TransitiveClosure() {
//____________________________________________________________________________
// Transitive Closure

	typedef VertexDescriptor<Graph<> >::Type TVertexDescriptor;
	typedef EdgeDescriptor<Graph<> >::Type TEdgeDescriptor;
	typedef Size<Graph<> >::Type TSize;

	//Number of edges
	TSize numEdges = 5;
	//Source, Target, Source, Target, Source, ...
	TVertexDescriptor edges[] = {3,0, 1,2, 2,1, 1,3, 3,2};
		
	//Create the graph
	Graph<> g;
	addEdges(g,edges, numEdges);

	// Transitive-Closure
	String<bool> closure;
	transitive_closure(g,closure);
	
	unsigned int len = (unsigned int) sqrt((double) length(closure));
	SEQAN_TASSERT(getValue(closure, 0*len + 0) == 1)
	SEQAN_TASSERT(getValue(closure, 0*len + 1) == 0)
	SEQAN_TASSERT(getValue(closure, 0*len + 2) == 0)
	SEQAN_TASSERT(getValue(closure, 0*len + 3) == 0)
	SEQAN_TASSERT(getValue(closure, 1*len + 0) == 1)
	SEQAN_TASSERT(getValue(closure, 1*len + 1) == 1)
	SEQAN_TASSERT(getValue(closure, 1*len + 2) == 1)
	SEQAN_TASSERT(getValue(closure, 1*len + 3) == 1)
	SEQAN_TASSERT(getValue(closure, 2*len + 0) == 1)
	SEQAN_TASSERT(getValue(closure, 2*len + 1) == 1)
	SEQAN_TASSERT(getValue(closure, 2*len + 2) == 1)
	SEQAN_TASSERT(getValue(closure, 2*len + 3) == 1)
	SEQAN_TASSERT(getValue(closure, 3*len + 0) == 1)
	SEQAN_TASSERT(getValue(closure, 3*len + 1) == 1)
	SEQAN_TASSERT(getValue(closure, 3*len + 2) == 1)
	SEQAN_TASSERT(getValue(closure, 3*len + 3) == 1)
}

//////////////////////////////////////////////////////////////////////////////

void Test_FordFulkerson() {
//____________________________________________________________________________
// Ford-Fulkerson
	typedef Graph<Directed<> > TGraph;
	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
	typedef Size<TGraph>::Type TSize;

	//Number of edges
	TSize numEdges = 10;
	//Source, Target, Source, Target, Source, ...
	TVertexDescriptor edges[] = {0,1, 0,4, 1,2, 1,4, 2,3, 2,4, 4,1, 4,5, 5,2, 5,3};
	unsigned int capacity[] =    {16,  13,  12,  10,  20,  9,   4,   14,  7,   4};

	//Create the graph
	Graph<> g;
	addEdges(g,edges, numEdges);
	String<unsigned int> capMap;	
	resizeEdgeMap(g,capMap, capacity);

	// Out-parameter
	String<unsigned int> flow;	
	unsigned int valF = ford_fulkerson(g, 0, 3, capMap, flow);
	
	SEQAN_TASSERT(valF == 23)
	TEdgeIterator itEdge(g);
	for(;!atEnd(itEdge);goNext(itEdge)) {
		SEQAN_TASSERT(getProperty(flow, getValue(itEdge)) <= getProperty(capMap, getValue(itEdge)))
	}
}

//////////////////////////////////////////////////////////////////////////////

void Test_PathGrowingAlgorithm() {
//____________________________________________________________________________
// Path growing algorithm
	typedef Graph<Undirected<> > TGraph;
	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
	typedef Size<TGraph>::Type TSize;

	//Number of edges
	TSize numEdges = 14;
	//Source, Target, Source, Target, Source, ...
	TVertexDescriptor edges[] = {0,7, 0,5, 0,8, 1,5, 1,6, 2,6, 2,5, 2,8, 3,8, 3,4, 4,5, 4,6, 4,7, 4,8};
	unsigned int weights[] =    {20,  5,   19,  6,   7,   10,  3,   9,   12,  11,  13,  12,  9,   12};

	//Create the graph
	TGraph g;
	addEdges(g,edges, numEdges);
	String<unsigned int> weightMap;
	resizeEdgeMap(g, weightMap, weights);

	// Path growing algorithm
	String<bool> edgeMap;	

	// EdgeMap indicates whether an edge is selected or not
	unsigned int weight = path_growing_algorithm(g, weightMap, edgeMap);

	SEQAN_TASSERT(weight == 49)
	SEQAN_TASSERT(getProperty(edgeMap, findEdge(g, 0, 7)) == true)
	SEQAN_TASSERT(getProperty(edgeMap, findEdge(g, 1, 6)) == true)
	SEQAN_TASSERT(getProperty(edgeMap, findEdge(g, 2, 8)) == true)
	SEQAN_TASSERT(getProperty(edgeMap, findEdge(g, 4, 5)) == true)
}


//////////////////////////////////////////////////////////////////////////////

void Test_LongestIncreasingSubsequence() {
	String<char> seq1("zeitgeist");
	String<unsigned int, Block<> > pos1;
	longestIncreasingSubsequence(seq1,pos1);
	// Trace is backwards
	SEQAN_TASSERT(seq1[pos1[4]] == 'e')
	SEQAN_TASSERT(seq1[pos1[3]] == 'g')
	SEQAN_TASSERT(seq1[pos1[2]] == 'i')
	SEQAN_TASSERT(seq1[pos1[1]] == 's')
	SEQAN_TASSERT(seq1[pos1[0]] == 't')

	String<unsigned int> seq;
	appendValue(seq, 5); appendValue(seq, 3); appendValue(seq, 4);
	appendValue(seq, 9); appendValue(seq, 6); appendValue(seq, 2);
	appendValue(seq, 1); appendValue(seq, 8); appendValue(seq, 7);
	appendValue(seq, 10);
	String<unsigned int, Block<> > pos;
	longestIncreasingSubsequence(seq,pos);
	SEQAN_TASSERT(seq[pos[4]] == 3)
	SEQAN_TASSERT(seq[pos[3]] == 4)
	SEQAN_TASSERT(seq[pos[2]] == 6)
	SEQAN_TASSERT(seq[pos[1]] == 7)
	SEQAN_TASSERT(seq[pos[0]] == 10)
}

//////////////////////////////////////////////////////////////////////////////

void Test_LongestCommonSubsequence() {
	String<char> seq1("abacx");
	String<char> seq2("baabca");
	String<std::pair<unsigned int, unsigned int> > pos;
	longestCommonSubsequence(seq1, seq2, 100, pos);
	SEQAN_TASSERT(seq1[pos[2].first] == 'b')
	SEQAN_TASSERT(seq2[pos[2].second] == 'b')
	SEQAN_TASSERT(seq1[pos[1].first] == 'a')
	SEQAN_TASSERT(seq2[pos[1].second] == 'a')
	SEQAN_TASSERT(seq1[pos[0].first] == 'c')
	SEQAN_TASSERT(seq2[pos[0].second] == 'c')
}


//////////////////////////////////////////////////////////////////////////////

void Test_HeaviestIncreasingSubsequence() {
	String<char> seq1("zeitgeist");
	String<unsigned int> weights1;
	String<unsigned int> pos1;
	fill(weights1, length(seq1), 1);
	unsigned int w = heaviestIncreasingSubsequence(seq1, weights1, pos1);
	// Trace is backwards
	SEQAN_TASSERT(w == 5)
	SEQAN_TASSERT(seq1[pos1[4]] == 'e')
	SEQAN_TASSERT(seq1[pos1[3]] == 'g')
	SEQAN_TASSERT(seq1[pos1[2]] == 'i')
	SEQAN_TASSERT(seq1[pos1[1]] == 's')
	SEQAN_TASSERT(seq1[pos1[0]] == 't')
	//// Output
	//for(int i = length(pos1)-1; i>=0; --i) {
	//	std::cout << seq1[pos1[i]] <<  ',';
	//}
	//std::cout << std::endl;

	// Alter weights
	clear(pos1);
	assignProperty(weights1, 2, 10);
	w = heaviestIncreasingSubsequence(seq1, weights1, pos1);
	SEQAN_TASSERT(w == 13)
	SEQAN_TASSERT(seq1[pos1[3]] == 'e')
	SEQAN_TASSERT(seq1[pos1[2]] == 'i')
	SEQAN_TASSERT(seq1[pos1[1]] == 's')
	SEQAN_TASSERT(seq1[pos1[0]] == 't')
	//// Output
	//for(int i = length(pos1)-1; i>=0; --i) {
	//	std::cout << seq1[pos1[i]] <<  ',';
	//}
	//std::cout << std::endl;

	String<unsigned int> seq;
	appendValue(seq, 1); appendValue(seq, 0);
	appendValue(seq, 1); appendValue(seq, 0);
	appendValue(seq, 1); appendValue(seq, 0);
	String<unsigned int> weights;
	appendValue(weights, 15); appendValue(weights, 10);
	appendValue(weights, 10); appendValue(weights, 10);
	appendValue(weights, 10); appendValue(weights, 15);
	String<unsigned int> pos;
	w = heaviestIncreasingSubsequence(seq, weights, pos);
	SEQAN_TASSERT(w == 20)
	SEQAN_TASSERT(seq[pos[1]] == 0)
	SEQAN_TASSERT(seq[pos[0]] == 1)
	//// Output
	//for(int i = length(pos)-1; i>=0; --i) {
	//	std::cout << seq[pos[i]] <<  ',';
	//}
	//std::cout << std::endl;
}


//////////////////////////////////////////////////////////////////////////////

void Test_HmmAlgorithms() {
	typedef double TProbability;
	typedef Dna TAlphabet;
	typedef Size<TAlphabet>::Type TSize;
	typedef Graph<Hmm<TAlphabet, TProbability> > THmm;
	typedef VertexDescriptor<THmm>::Type TVertexDescriptor;
	typedef EdgeDescriptor<THmm>::Type TEdgeDescriptor;
	TSize alph_size = ValueSize<TAlphabet>::VALUE;
	
	Dna dnaA = Dna('A');
	Dna dnaC = Dna('C');
	Dna dnaG = Dna('G');
	Dna dnaT = Dna('T');

	THmm hmm;
	TVertexDescriptor state1 = addVertex(hmm);
	emissionProbability(hmm, state1, dnaA) = 0.2;
	emissionProbability(hmm, state1, dnaC) = 0.2;
	emissionProbability(hmm, state1, dnaG) = 0.3;
	emissionProbability(hmm, state1, dnaT) = 0.3;
	String<TProbability> emis;
	resize(emis, alph_size);
	value(emis, ordValue(dnaA)) = 0.5;
	value(emis, ordValue(dnaC)) = 0.5;
	value(emis, ordValue(dnaG)) = 0.0;
	value(emis, ordValue(dnaT)) = 0.0;
	TVertexDescriptor state2 = addVertex(hmm, emis);
	TVertexDescriptor state3 = addVertex(hmm, emis);
	assignEmissionProbability(hmm, state3, dnaA, 0.3);
	assignEmissionProbability(hmm, state3, dnaC, 0.3);
	assignEmissionProbability(hmm, state3, dnaG, 0.2);
	assignEmissionProbability(hmm, state3, dnaT, 0.2);
	addEdge(hmm, state1, state1, 0.95);
	addEdge(hmm, state1, state3, 0.05);
	TVertexDescriptor begState = addVertex(hmm);
	TVertexDescriptor eState = addVertex(hmm);
	addEdge(hmm, begState, state1, 1.0);
	addEdge(hmm, state3, eState, 0.5);
	addEdge(hmm, eState, eState, 1.0);
	assignBeginState(hmm, begState);
	assignEndState(hmm, eState);
	removeVertex(hmm, state2);

	// Algorithms
	String<Dna> sequence = "AC";
	String<TVertexDescriptor> path;
	viterbiAlgorithm(hmm, sequence, path);
	forwardAlgorithm(hmm, sequence);
	backwardAlgorithm(hmm, sequence);
}


//////////////////////////////////////////////////////////////////////////////

void Test_GraphAlgorithms() {
	mtRandInit();
//____________________________________________________________________________
// Graph Algorithms
	// Test HeapTree
	Test_HeapTree();

	// Elementary graph algorithms
	Test_BreadthFirstSearch();
	Test_DepthFirstSearch();
	Test_TopologicalSort();
	Test_StronglyConnectedComponents();
	Test_ConnectedComponents();

	// Minimum Spanning Trees
	Test_PrimsAlgorithm();
	Test_KruskalsAlgorithm();
	Test_MST_All();

	// Single-Source shortest paths
	Test_DagShortestPath();
	Test_BellmanFord();
	Test_Dijkstra();

	// All-Pairs Shortest paths
	Test_AllPairsShortestPath();
	Test_FloydWarshall();
	Test_TransitiveClosure();

	//Maximum Flow
	Test_FordFulkerson();

	//Matching
	Test_PathGrowingAlgorithm();

	// Lis, lcs, his, hcs
	Test_LongestIncreasingSubsequence();
	Test_LongestCommonSubsequence();
	Test_HeaviestIncreasingSubsequence();
	// ToDo: Generic heaviest common subsequence

	// Hmm algorithms
	Test_HmmAlgorithms();

	debug::verifyCheckpoints("projects/library/seqan/graph_algorithms/graph_algorithm.h");
	debug::verifyCheckpoints("projects/library/seqan/graph_algorithms/graph_algorithm_hmm.h");
	debug::verifyCheckpoints("projects/library/seqan/graph_algorithms/graph_algorithm_lis_his.h");
}


}

#endif


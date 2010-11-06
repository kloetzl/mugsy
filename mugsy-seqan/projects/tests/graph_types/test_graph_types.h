#ifndef SEQAN_HEADER_TEST_GRAPH_TYPES_H
#define SEQAN_HEADER_TEST_GRAPH_TYPES_H

namespace SEQAN_NAMESPACE_MAIN
{

void Test_Directed() {
//____________________________________________________________________________
// Graph without edge cargo but with edge ids

	typedef Graph<Directed<> > StandardGraph;
	typedef VertexDescriptor<StandardGraph>::Type TVertexDescriptor;
	typedef EdgeDescriptor<StandardGraph>::Type TEdgeDescriptor;
	
	StandardGraph g;
	SEQAN_TASSERT(numVertices(g) == 0)
	SEQAN_TASSERT(numEdges(g) == 0)
	SEQAN_TASSERT(empty(g) == true)

	// Add vertex
	TVertexDescriptor v0 = addVertex(g);
	SEQAN_TASSERT(v0 == 0)
	SEQAN_TASSERT(outDegree(g, v0) == 0)	
	SEQAN_TASSERT(inDegree(g, 0) == 0)
	SEQAN_TASSERT(degree(g, 0) == 0)
	SEQAN_TASSERT(numVertices(g) == 1)
	SEQAN_TASSERT(empty(g) == false)
	
	// Add edge
	TEdgeDescriptor e1 =addEdge(g,v0,v0);
	SEQAN_TASSERT(findEdge(g, v0, v0) == e1)
	SEQAN_TASSERT(_getVertexString(g)[0] == e1)
	SEQAN_TASSERT(getIdUpperBound(_getVertexIdManager(g)) == 1)
	SEQAN_TASSERT(getIdUpperBound(_getEdgeIdManager(g)) == 1)
	SEQAN_TASSERT(targetVertex(g, e1) == 0)
	SEQAN_TASSERT(sourceVertex(g, e1) == 0)  //Expensive in standard graph!
	SEQAN_TASSERT(numEdges(g) == 1)
	SEQAN_TASSERT(outDegree(g, v0) == 1)	
	SEQAN_TASSERT(inDegree(g, v0) == 1)
	SEQAN_TASSERT(degree(g, v0) == 2)
	
	// Add further edges and vertices
	TVertexDescriptor v1 = addVertex(g);
	TEdgeDescriptor e2 =addEdge(g,0,1);
	SEQAN_TASSERT(v1 == 1)
	SEQAN_TASSERT(numVertices(g) == 2)
	SEQAN_TASSERT(targetVertex(g, e2) == 1)
	SEQAN_TASSERT(sourceVertex(g, e2) == 0)
	SEQAN_TASSERT(numEdges(g) == 2)
	SEQAN_TASSERT(outDegree(g, v0) == 2)	
	SEQAN_TASSERT(inDegree(g, 1) == 1)
	SEQAN_TASSERT(inDegree(g, 0) == 1)	
	SEQAN_TASSERT(degree(g, 0) == 3)
		
	// Add more vertices and edges
	addVertex(g);  //2
	TVertexDescriptor v3 = addVertex(g);  //3
	addVertex(g);  //4
	addEdge(g,3,4);
	TEdgeDescriptor my_edge = addEdge(g,3,1);
	addEdge(g,3,0);
	SEQAN_TASSERT(v3 == 3)
	SEQAN_TASSERT(numVertices(g) == 5)
	SEQAN_TASSERT(targetVertex(g, e2) == 1)
	SEQAN_TASSERT(sourceVertex(g, e2) == 0)
	SEQAN_TASSERT(targetVertex(g, my_edge) == 1)
	SEQAN_TASSERT(sourceVertex(g, my_edge) == 3)
	SEQAN_TASSERT(numEdges(g) == 5)
	SEQAN_TASSERT(outDegree(g, v3) == 3)	
	
	// Output
	std::cout << g << std::endl;

	// Remove edges
	removeEdge(g,my_edge);
	removeEdge(g,0,1);
	SEQAN_TASSERT(numEdges(g) == 3)

	// Remove vertices 
	TEdgeDescriptor e3 = addEdge(g,3,3);
	addEdge(g,1,3);
	addEdge(g,0,3);
	addEdge(g,0,4);
	SEQAN_TASSERT(outDegree(g, 0) == 3)
	SEQAN_TASSERT(outDegree(g, 1) == 1)
	SEQAN_TASSERT(targetVertex(g, e3) == 3)
	SEQAN_TASSERT(sourceVertex(g, e3) == 3)
	removeVertex(g, v3);
	SEQAN_TASSERT(outDegree(g, 0) == 2)
	SEQAN_TASSERT(outDegree(g, 1) == 0)
	SEQAN_TASSERT(numVertices(g) == 4)
	SEQAN_TASSERT(numEdges(g) == 2)

	// Clear graph
	clearEdges(g);
	SEQAN_TASSERT(numVertices(g) == 4)
	SEQAN_TASSERT(numEdges(g) == 0)
	addEdge(g,2,0);
	addEdge(g,4,1);
	clearVertices(g);
	SEQAN_TASSERT(numVertices(g) == 0)
	SEQAN_TASSERT(numEdges(g) == 0)
	addVertex(g);addVertex(g);addVertex(g);
	addVertex(g);addVertex(g);
	addEdge(g,2,0);
	addEdge(g,4,1);
	clear(g);
	SEQAN_TASSERT(numVertices(g) == 0)
	SEQAN_TASSERT(numEdges(g) == 0)
	addVertex(g);addVertex(g);addVertex(g);
	addVertex(g);addVertex(g);
	addEdge(g,2,0);
	addEdge(g,4,1);
	addEdge(g,4,2);
	removeVertex(g,3);
	SEQAN_TASSERT(numVertices(g) == 4)
	SEQAN_TASSERT(numEdges(g) == 3)
	SEQAN_TASSERT(outDegree(g, 4) == 2)
	SEQAN_TASSERT(inDegree(g, 4) == 0)

	// Transpose
	transpose(g); 
	SEQAN_TASSERT(numVertices(g) == 4)
	SEQAN_TASSERT(numEdges(g) == 3)
	SEQAN_TASSERT(outDegree(g, 4) == 0)
	SEQAN_TASSERT(inDegree(g, 4) == 2)
	StandardGraph g_copy(g);
	SEQAN_TASSERT(numVertices(g_copy) == 4)
	SEQAN_TASSERT(numEdges(g_copy) == 3)
	SEQAN_TASSERT(outDegree(g_copy, 4) == 0)
	SEQAN_TASSERT(inDegree(g_copy, 4) == 2)
	addVertex(g_copy);
	addEdge(g_copy,3,0);
	g_copy = g;
	SEQAN_TASSERT(numVertices(g_copy) == 4)
	SEQAN_TASSERT(numEdges(g_copy) == 3)
	SEQAN_TASSERT(outDegree(g_copy, 4) == 0)
	SEQAN_TASSERT(inDegree(g_copy, 4) == 2)
	//Copies the graph and transposes just the copy
	transpose(g,g_copy);  // g does not change!
	SEQAN_TASSERT(numVertices(g_copy) == 4)
	SEQAN_TASSERT(numEdges(g_copy) == 3)
	SEQAN_TASSERT(outDegree(g_copy, 4) == 2)
	SEQAN_TASSERT(inDegree(g_copy, 4) == 0)
	removeVertex(g,0);



	// Adjacency matrix
	String<unsigned int> mat;
	getAdjacencyMatrix(g, mat);
	unsigned int len = (unsigned int) std::sqrt((double) length(mat));
	SEQAN_TASSERT(getValue(mat, 1*len+4) == 1)
	SEQAN_TASSERT(getValue(mat, 2*len+4) == 1)
	SEQAN_TASSERT(getValue(mat, 2*len+2) == 0)

//____________________________________________________________________________
//Graph with edge cargo and edge ids
	typedef Pair<char, int> TPair;
	typedef Directed<TPair> TEdges;
	typedef VertexDescriptor<Graph<TEdges> >::Type TVertexDescriptor2;
	typedef EdgeDescriptor<Graph<TEdges> >::Type TEdgeDescriptor2;

	Graph<TEdges> g2;
	SEQAN_TASSERT(numVertices(g2) == 0)
	SEQAN_TASSERT(numEdges(g2) == 0)
	TVertexDescriptor2 ver0 = addVertex(g2);
	SEQAN_TASSERT(ver0 == 0)
	SEQAN_TASSERT(numVertices(g2) == 1)
	TVertexDescriptor2 ver1 = addVertex(g2);
	SEQAN_TASSERT(ver1 == 1)
	SEQAN_TASSERT(numVertices(g2) == 2)
	TEdgeDescriptor2 ed1 =addEdge(g2,ver0,ver0, TPair('a',3));
	TEdgeDescriptor2 ed2 =addEdge(g2,0,1);
	SEQAN_TASSERT((getCargo(ed1)).i1 == 'a')
	SEQAN_TASSERT((getCargo(ed1)).i2 == 3)
	SEQAN_TASSERT(targetVertex(g2, ed1) == v0)
	SEQAN_TASSERT(targetVertex(g2, ed1) == 0)
	SEQAN_TASSERT(sourceVertex(g2, ed1) == 0)
	SEQAN_TASSERT(targetVertex(g2, ed2) == 1)
	SEQAN_TASSERT(numEdges(g2) == 2)
	assignCargo(ed2, TPair('b',4));
	SEQAN_TASSERT((getCargo(ed1)).i1 == 'a')
	SEQAN_TASSERT((getCargo(ed1)).i2 == 3)
	SEQAN_TASSERT((getCargo(ed2)).i1 == 'b')
	SEQAN_TASSERT((getCargo(ed2)).i2 == 4)
	cargo(ed1) = TPair('c',1);
	SEQAN_TASSERT((getCargo(ed1)).i1 == 'c')
	SEQAN_TASSERT((getCargo(ed1)).i2 == 1)
	addVertex(g2);
	addVertex(g2);
	addVertex(g2);
	TEdgeDescriptor2 ed4 =addEdge(g2,1,4);
	cargo(ed4) = TPair('z',100);
	removeVertex(g2, 2);
	Graph<TEdges> g2_copy(g2);
	SEQAN_TASSERT(numVertices(g2_copy) == 4)
	SEQAN_TASSERT(numEdges(g2_copy) == 3)
	clearEdges(g2_copy);
	SEQAN_TASSERT(numVertices(g2_copy) == 4)
	SEQAN_TASSERT(numEdges(g2_copy) == 0)
	clearVertices(g2_copy);
	SEQAN_TASSERT(numVertices(g2_copy) == 0)
	addVertex(g2_copy);addVertex(g2_copy);
	addEdge(g2_copy,0,1);
	clear(g2_copy);
	SEQAN_TASSERT(numVertices(g2_copy) == 0)
	addVertex(g2_copy);addVertex(g2_copy);
	addEdge(g2_copy,0,1);
	SEQAN_TASSERT(numEdges(g2) == 3)
	SEQAN_TASSERT(outDegree(g2, 0) == 2)
	SEQAN_TASSERT(inDegree(g2, 0) == 1)
	transpose(g2, g2_copy);
	SEQAN_TASSERT(outDegree(g2_copy, 0) == 1)
	SEQAN_TASSERT(inDegree(g2_copy, 0) == 2)
	SEQAN_TASSERT(numEdges(g2_copy) == 3)
	TEdgeDescriptor2 edgCargo = addEdge(g2, 0, 0, TPair('m',3));
	SEQAN_TASSERT((getCargo(edgCargo)).i1 == 'm')
	SEQAN_TASSERT((getCargo(edgCargo)).i2 == 3)

//____________________________________________________________________________
//Graph without edge cargo and without edge ids
	typedef Directed<void, WithoutEdgeId> TEdges3;
	typedef VertexDescriptor<Graph<TEdges3> >::Type TVertexDescriptor3;
	typedef EdgeDescriptor<Graph<TEdges3> >::Type TEdgeDescriptor3;

	Graph<TEdges3> g3;
	addVertex(g3);addVertex(g3);addVertex(g3);
	addVertex(g3);addVertex(g3);
	addEdge(g3,1,4);
	SEQAN_TASSERT(numVertices(g3) == 5)
	SEQAN_TASSERT(numEdges(g3) == 1)
	TEdgeDescriptor3 edge3 = addEdge(g3,0,4);
	//SEQAN_TASSERT(_getId(edge3) == 0);
	SEQAN_TASSERT(getCargo(edge3) == (void*) 0);
	addEdge(g3,0,2);
	addEdge(g3,0,0);
	removeEdge(g3,0,4);
	removeEdge(g3,0,2);
	SEQAN_TASSERT(numEdges(g3) == 2)

	// Multigraph
	StandardGraph multiG;
	addVertex(multiG);addVertex(multiG);addVertex(multiG);
	TEdgeDescriptor edgeD1 = addEdge(multiG, 1, 2);
	TEdgeDescriptor edgeD2 = addEdge(multiG, 1, 2);
	TEdgeDescriptor edgeD3 = addEdge(multiG, 1, 2);
	removeEdge(multiG, edgeD2);
	SEQAN_TASSERT(sourceVertex(multiG,edgeD1) == 1)
	SEQAN_TASSERT(sourceVertex(multiG,edgeD2) == 0) // EdgeDescriptor invalid
	SEQAN_TASSERT(sourceVertex(multiG,edgeD3) == 1)
}

//////////////////////////////////////////////////////////////////////////////

void Test_Undirected() {
//____________________________________________________________________________
// Graph without edge cargo but with edge ids

	typedef Graph<Undirected<void> > StandardGraph;
	typedef VertexDescriptor<StandardGraph>::Type TVertexDescriptor;
	typedef EdgeDescriptor<StandardGraph>::Type TEdgeDescriptor;
	
	StandardGraph g;
	SEQAN_TASSERT(numVertices(g) == 0)
	SEQAN_TASSERT(numEdges(g) == 0)
	SEQAN_TASSERT(empty(g) == true)

	// Add vertex
	TVertexDescriptor v0 = addVertex(g);
	SEQAN_TASSERT(v0 == 0)
	SEQAN_TASSERT(outDegree(g, v0) == 0)	
	SEQAN_TASSERT(inDegree(g, 0) == 0)
	SEQAN_TASSERT(degree(g, 0) == 0)
	SEQAN_TASSERT(numVertices(g) == 1)
	SEQAN_TASSERT(empty(g) == false)
	
	// Add edge
	// TEdgeDescriptor e1 =addEdge(g,v0,v0);  // Self edges are not allowed in undirected graphs
	TVertexDescriptor v1 = addVertex(g);
	TEdgeDescriptor e =addEdge(g,0,1);
	SEQAN_TASSERT(findEdge(g, 0, 1) == e)
	SEQAN_TASSERT(_getVertexString(g)[0] == e)
	SEQAN_TASSERT(getIdUpperBound(_getVertexIdManager(g)) == 2)
	SEQAN_TASSERT(getIdUpperBound(_getEdgeIdManager(g)) == 1)
	SEQAN_TASSERT(v1 == 1)
	SEQAN_TASSERT(numVertices(g) == 2)
	SEQAN_TASSERT(targetVertex(g, e) == 1)
	SEQAN_TASSERT(sourceVertex(g, e) == 0)
	SEQAN_TASSERT(numEdges(g) == 1)
	SEQAN_TASSERT(outDegree(g, v0) == 1)	
	SEQAN_TASSERT(inDegree(g, 1) == 1)
	SEQAN_TASSERT(inDegree(g, 0) == 1)	
	SEQAN_TASSERT(degree(g, 0) == 1)

	// Add more vertices and edges
	addVertex(g);  //2
	TVertexDescriptor v3 = addVertex(g);  //3
	addVertex(g);  //4
	addEdge(g,3,4);
	TEdgeDescriptor my_edge = addEdge(g,3,1);
	addEdge(g,3,0);
	SEQAN_TASSERT(v3 == 3)
	SEQAN_TASSERT(numVertices(g) == 5)
	SEQAN_TASSERT(targetVertex(g, my_edge) == 3)
	SEQAN_TASSERT(sourceVertex(g, my_edge) == 1)
	SEQAN_TASSERT(numEdges(g) == 4)
	SEQAN_TASSERT(outDegree(g, v3) == 3)
	SEQAN_TASSERT(inDegree(g, v3) == 3)
	SEQAN_TASSERT(degree(g, v3) == 3)

	// Output
	std::cout << g << std::endl;

	// Remove edges
	removeEdge(g,my_edge);
	removeEdge(g,0,1);
	SEQAN_TASSERT(numEdges(g) == 2)

	
	// Remove vertices 
	addVertex(g);  //5
	addEdge(g,5,2);
	addEdge(g,2,3);
	addEdge(g,1,3);
	addEdge(g,1,4);
	SEQAN_TASSERT(outDegree(g, 3) == 4)
	SEQAN_TASSERT(outDegree(g, 4) == 2)
	removeVertex(g, v3);
	SEQAN_TASSERT(outDegree(g, 4) == 1)
	SEQAN_TASSERT(outDegree(g, 0) == 0)
	SEQAN_TASSERT(numVertices(g) == 5)
	SEQAN_TASSERT(numEdges(g) == 2)

	// Clear graph
	clearEdges(g);
	SEQAN_TASSERT(numVertices(g) == 5)
	SEQAN_TASSERT(numEdges(g) == 0)
	addEdge(g,2,0);
	addEdge(g,4,1);
	clearVertices(g);
	SEQAN_TASSERT(numVertices(g) == 0)
	SEQAN_TASSERT(numEdges(g) == 0)
	addVertex(g);addVertex(g);addVertex(g);
	addVertex(g);addVertex(g);
	addEdge(g,2,0);
	addEdge(g,4,1);
	clear(g);
	SEQAN_TASSERT(numVertices(g) == 0)
	SEQAN_TASSERT(numEdges(g) == 0)
	addVertex(g);addVertex(g);addVertex(g);
	addVertex(g);addVertex(g);
	addEdge(g,2,0);
	addEdge(g,4,1);
	addEdge(g,4,2);
	removeVertex(g,3);
	SEQAN_TASSERT(numVertices(g) == 4)
	SEQAN_TASSERT(numEdges(g) == 3)
	SEQAN_TASSERT(outDegree(g, 4) == 2)
	SEQAN_TASSERT(inDegree(g, 4) == 2)

	// Transpose
	transpose(g); 
	SEQAN_TASSERT(numVertices(g) == 4)
	SEQAN_TASSERT(numEdges(g) == 3)
	SEQAN_TASSERT(outDegree(g, 4) == 2)
	SEQAN_TASSERT(inDegree(g, 4) == 2)
	StandardGraph g_copy(g);
	SEQAN_TASSERT(numVertices(g_copy) == 4)
	SEQAN_TASSERT(numEdges(g_copy) == 3)
	SEQAN_TASSERT(outDegree(g_copy, 4) == 2)
	SEQAN_TASSERT(inDegree(g_copy, 4) == 2)
	addVertex(g_copy);
	addEdge(g_copy,3,0);
	g_copy = g;
	SEQAN_TASSERT(numVertices(g_copy) == 4)
	SEQAN_TASSERT(numEdges(g_copy) == 3)
	SEQAN_TASSERT(outDegree(g_copy, 4) == 2)
	SEQAN_TASSERT(inDegree(g_copy, 4) == 2)
	//Copies the graph and transposes just the copy
	transpose(g,g_copy);  // g does not change!
	SEQAN_TASSERT(numVertices(g_copy) == 4)
	SEQAN_TASSERT(numEdges(g_copy) == 3)
	SEQAN_TASSERT(outDegree(g_copy, 4) == 2)
	SEQAN_TASSERT(inDegree(g_copy, 4) == 2)

	// Adjacency matrix
	String<unsigned int> mat;
	getAdjacencyMatrix(g, mat);
	unsigned int len = (unsigned int) std::sqrt((double) length(mat));
	SEQAN_TASSERT(getValue(mat,0*len+2) == 1)
	SEQAN_TASSERT(getValue(mat,3*len+2) == 0)
	SEQAN_TASSERT(getValue(mat,0*len+2) == getValue(mat,2*len+0))
	SEQAN_TASSERT(getValue(mat,1*len+4) == getValue(mat,4*len+1))
	SEQAN_TASSERT(getValue(mat,2*len+4) == getValue(mat,4*len+2))

//____________________________________________________________________________
//Graph with edge cargo and edge ids
	typedef Pair<char, int> TPair;
	typedef Undirected<TPair> TEdges;
	typedef VertexDescriptor<Graph<TEdges> >::Type TVertexDescriptor2;
	typedef EdgeDescriptor<Graph<TEdges> >::Type TEdgeDescriptor2;

	Graph<TEdges> g2;
	SEQAN_TASSERT(numVertices(g2) == 0)
	SEQAN_TASSERT(numEdges(g2) == 0)
	TVertexDescriptor2 ver0 = addVertex(g2);
	SEQAN_TASSERT(ver0 == 0)
	SEQAN_TASSERT(numVertices(g2) == 1)
	TVertexDescriptor2 ver1 = addVertex(g2);
	SEQAN_TASSERT(ver1 == 1)
	SEQAN_TASSERT(numVertices(g2) == 2)
	TEdgeDescriptor2 ed1 =addEdge(g2,0,1);
	SEQAN_TASSERT(targetVertex(g2, ed1) == 1)
	SEQAN_TASSERT(sourceVertex(g2, ed1) == 0)
	SEQAN_TASSERT(numEdges(g2) == 1)
	assignCargo(ed1, TPair('a',3));
	SEQAN_TASSERT((getCargo(ed1)).i1 == 'a')
	SEQAN_TASSERT((getCargo(ed1)).i2 == 3)
	cargo(ed1) = TPair('c',1);
	SEQAN_TASSERT((getCargo(ed1)).i1 == 'c')
	SEQAN_TASSERT((getCargo(ed1)).i2 == 1)
	addVertex(g2);
	addVertex(g2);
	addVertex(g2);
	TEdgeDescriptor2 ed4 =addEdge(g2,1,4);
	cargo(ed4) = TPair('z',100);
	removeVertex(g2, 2);
	Graph<TEdges> g2_copy(g2);
	SEQAN_TASSERT(numVertices(g2_copy) == 4)
	SEQAN_TASSERT(numEdges(g2_copy) == 2)
	clearEdges(g2_copy);
	SEQAN_TASSERT(numVertices(g2_copy) == 4)
	SEQAN_TASSERT(numEdges(g2_copy) == 0)
	clearVertices(g2_copy);
	SEQAN_TASSERT(numVertices(g2_copy) == 0)
	addVertex(g2_copy);addVertex(g2_copy);
	addEdge(g2_copy,0,1);
	clear(g2_copy);
	SEQAN_TASSERT(numVertices(g2_copy) == 0)
	addVertex(g2_copy);addVertex(g2_copy);
	addEdge(g2_copy,0,1);
	transpose(g2, g2_copy);
	SEQAN_TASSERT(outDegree(g2_copy, 0) == 1)
	SEQAN_TASSERT(inDegree(g2_copy, 0) == 1)
	SEQAN_TASSERT(numEdges(g2_copy) == 2)
	TEdgeDescriptor2 edgCargo = addEdge(g2, 0, 3, TPair('m',3));
	SEQAN_TASSERT((getCargo(edgCargo)).i1 == 'm')
	SEQAN_TASSERT((getCargo(edgCargo)).i2 == 3)

//____________________________________________________________________________
//Graph without edge cargo and without edge ids
	typedef Undirected<void, WithoutEdgeId> TEdges3;
	typedef VertexDescriptor<Graph<TEdges3> >::Type TVertexDescriptor3;
	typedef EdgeDescriptor<Graph<TEdges3> >::Type TEdgeDescriptor3;

	Graph<TEdges3> g3;
	addVertex(g3);addVertex(g3);addVertex(g3);
	addVertex(g3);addVertex(g3);
	addEdge(g3,1,4);
	SEQAN_TASSERT(numVertices(g3) == 5)
	SEQAN_TASSERT(numEdges(g3) == 1)
	TEdgeDescriptor3 edge3 = addEdge(g3,0,4);
	SEQAN_TASSERT(_getId(edge3) == 0);
	SEQAN_TASSERT(getCargo(edge3) == (void*) 0);
	addEdge(g3,0,2);
	addEdge(g3,0,1);
	removeEdge(g3,0,4);
	removeEdge(g3,0,2);
	SEQAN_TASSERT(numEdges(g3) == 2)
	removeInEdges(g3,1);
	SEQAN_TASSERT(numEdges(g3) == 0)
	

//____________________________________________________________________________
// Undirected graph iterators
	typedef Graph<Undirected<> > TGraphIter;
	typedef VertexDescriptor<TGraphIter>::Type TVertexDescriptorIter;
	typedef EdgeDescriptor<TGraphIter>::Type TEdgeDescriptorIter;
	
	TGraphIter gIter;
	addVertex(gIter);addVertex(gIter);addVertex(gIter);addVertex(gIter);
	addVertex(gIter);addVertex(gIter);addVertex(gIter);addVertex(gIter);
	removeVertex(gIter,0);
	removeVertex(gIter,5);
	addEdge(gIter,2,7);
	addEdge(gIter,2,3);
	addEdge(gIter,2,4);
	addEdge(gIter,4,3);
	addEdge(gIter,3,6);
	addEdge(gIter,4,6);

	typedef Iterator<TGraphIter, OutEdgeIterator>::Type TOutEdgeIterator;
	TOutEdgeIterator itOutEdge(gIter,3);
	// Both ways are fast for undirected graphs
	SEQAN_TASSERT(sourceVertex(gIter, getValue(itOutEdge))==3)
	SEQAN_TASSERT(targetVertex(gIter, getValue(itOutEdge))==6)
	SEQAN_TASSERT(sourceVertex(itOutEdge)==3)
	SEQAN_TASSERT(targetVertex(itOutEdge)==6)
	SEQAN_TASSERT(sourceVertex(gIter, value(itOutEdge))==3)
	SEQAN_TASSERT(targetVertex(gIter, *itOutEdge)==6)
	SEQAN_TASSERT(atEnd(itOutEdge)==false)
	SEQAN_TASSERT(atBegin(itOutEdge)==true)
	goNext(itOutEdge);
	SEQAN_TASSERT(atEnd(itOutEdge)==false)
	SEQAN_TASSERT(atBegin(itOutEdge)==false)
	SEQAN_TASSERT(sourceVertex(itOutEdge)==3)
	SEQAN_TASSERT(targetVertex(itOutEdge)==4)
	++itOutEdge;
	itOutEdge++;
	SEQAN_TASSERT(atEnd(itOutEdge)==true)
	SEQAN_TASSERT(atBegin(itOutEdge)==false)
	goPrevious(itOutEdge);
	SEQAN_TASSERT(sourceVertex(gIter, getValue(itOutEdge))==2)
	SEQAN_TASSERT(targetVertex(gIter, getValue(itOutEdge))==3)
	--itOutEdge;
	itOutEdge--; 
	SEQAN_TASSERT(sourceVertex(gIter, getValue(itOutEdge))==3)
	SEQAN_TASSERT(targetVertex(gIter, getValue(itOutEdge))==6)
	itOutEdge--;
	itOutEdge--;
	SEQAN_TASSERT(atBegin(itOutEdge)==true)
	TOutEdgeIterator itEdge2(itOutEdge);
	TOutEdgeIterator itEdge3;
	itEdge3 = itOutEdge;
	SEQAN_TASSERT(itOutEdge == itEdge2)
	SEQAN_TASSERT(itEdge2 == itEdge3)
	goEnd(itOutEdge);
	SEQAN_TASSERT(itEdge2 != itOutEdge)
	goEnd(itEdge2);
	SEQAN_TASSERT(itEdge2 == itOutEdge)
	goBegin(itEdge2);
	SEQAN_TASSERT(atBegin(itEdge2)==true)
	SEQAN_TASSERT(&gIter == &hostGraph(itOutEdge))

	
	typedef Iterator<TGraphIter, EdgeIterator>::Type TEdgeIterator;
	TEdgeIterator itEdge(gIter);
	SEQAN_TASSERT(sourceVertex(gIter, getValue(itEdge))==2)
	SEQAN_TASSERT(targetVertex(gIter, getValue(itEdge))==4)
	SEQAN_TASSERT(atBegin(itEdge)==true)
	SEQAN_TASSERT(atEnd(itEdge)==false)
	goNext(itEdge);
	SEQAN_TASSERT(sourceVertex(gIter, getValue(itEdge))==2)
	SEQAN_TASSERT(targetVertex(gIter, getValue(itEdge))==3)
	SEQAN_TASSERT(atBegin(itEdge)==false)
	SEQAN_TASSERT(atEnd(itEdge)==false)
	goNext(itEdge);
	SEQAN_TASSERT(sourceVertex(gIter, getValue(itEdge))==2)
	SEQAN_TASSERT(targetVertex(gIter, getValue(itEdge))==7)
	++itEdge;
	SEQAN_TASSERT(sourceVertex(gIter, getValue(itEdge))==3)
	SEQAN_TASSERT(targetVertex(gIter, getValue(itEdge))==6)
	itEdge++;
	SEQAN_TASSERT(sourceVertex(gIter, getValue(itEdge))==3)
	SEQAN_TASSERT(targetVertex(gIter, getValue(itEdge))==4)
	goNext(itEdge);
	SEQAN_TASSERT(sourceVertex(gIter, getValue(itEdge))==4)
	SEQAN_TASSERT(targetVertex(gIter, getValue(itEdge))==6)
	goNext(itEdge);
	SEQAN_TASSERT(atBegin(itEdge)==false)
	SEQAN_TASSERT(atEnd(itEdge)==true)
	goPrevious(itEdge);	
	SEQAN_TASSERT(sourceVertex(gIter, getValue(itEdge))==4)
	SEQAN_TASSERT(targetVertex(gIter, getValue(itEdge))==6)
	--itEdge;	
	SEQAN_TASSERT(sourceVertex(gIter, getValue(itEdge))==3)
	SEQAN_TASSERT(targetVertex(gIter, getValue(itEdge))==4)
	itEdge--;	
	SEQAN_TASSERT(sourceVertex(gIter, getValue(itEdge))==3)
	SEQAN_TASSERT(targetVertex(gIter, getValue(itEdge))==6)
	goPrevious(itEdge);	
	SEQAN_TASSERT(sourceVertex(gIter, getValue(itEdge))==2)
	SEQAN_TASSERT(targetVertex(gIter, getValue(itEdge))==7)
	goPrevious(itEdge);	
	SEQAN_TASSERT(sourceVertex(gIter, getValue(itEdge))==2)
	SEQAN_TASSERT(targetVertex(gIter, getValue(itEdge))==3)
	goPrevious(itEdge);	
	SEQAN_TASSERT(sourceVertex(gIter, getValue(itEdge))==2)
	SEQAN_TASSERT(targetVertex(gIter, getValue(itEdge))==4)
	SEQAN_TASSERT(atBegin(itEdge)==true)
	SEQAN_TASSERT(atEnd(itEdge)==false)

	// Multigraph
	StandardGraph multiG;
	addVertex(multiG);addVertex(multiG);addVertex(multiG);
	addEdge(multiG, 1, 2);
	TEdgeDescriptor edgeD2 = addEdge(multiG, 1, 2);
	addEdge(multiG, 1, 2);
	removeEdge(multiG, edgeD2);
	SEQAN_TASSERT(numEdges(multiG) == 2)
}


//////////////////////////////////////////////////////////////////////////////

void Test_Automaton() {
//____________________________________________________________________________
// Standard automaton: No edge cargo

	typedef Graph<Automaton<Dna> > StandardAutomaton;
	typedef VertexDescriptor<StandardAutomaton>::Type TVertexDescriptor;
	typedef EdgeDescriptor<StandardAutomaton>::Type TEdgeDescriptor;
	
	StandardAutomaton g;
	SEQAN_TASSERT(numVertices(g) == 0)
	SEQAN_TASSERT(numEdges(g) == 0)
	SEQAN_TASSERT(empty(g) == true)

	// Add vertex
	createRoot(g);
	TVertexDescriptor v0 = getRoot(g);
	SEQAN_TASSERT(v0 == 0)
	SEQAN_TASSERT(outDegree(g, v0) == 0)	
	SEQAN_TASSERT(inDegree(g, 0) == 0)
	SEQAN_TASSERT(degree(g, 0) == 0)
	SEQAN_TASSERT(numVertices(g) == 1)
	SEQAN_TASSERT(empty(g) == false)

	// Add edge
	TEdgeDescriptor e1 =addEdge(g,v0,v0,'a');
	SEQAN_TASSERT(findEdge(g, 0, 'a') == e1)
	SEQAN_TASSERT(&_getVertexString(g)[0].data_edge[0] == e1)
	SEQAN_TASSERT(getIdUpperBound(_getVertexIdManager(g)) == 1)
	SEQAN_TASSERT(getIdUpperBound(_getEdgeIdManager(g)) == 1)
	SEQAN_TASSERT(_getId(e1) == 0)
	SEQAN_TASSERT(_getId(e1) == 0)
	SEQAN_TASSERT(targetVertex(g, e1) == 0)
	SEQAN_TASSERT(sourceVertex(g, e1) == 0) 
	SEQAN_TASSERT(numEdges(g) == 1)
	SEQAN_TASSERT(outDegree(g, v0) == 1)	
	SEQAN_TASSERT(inDegree(g, v0) == 1)
	SEQAN_TASSERT(degree(g, v0) == 2)
	
	// Add further edges and vertices
	TVertexDescriptor v1 = addVertex(g);
	TEdgeDescriptor e2 =addEdge(g,0,1,'g');
	SEQAN_TASSERT(_getId(e2) == 1)
	SEQAN_TASSERT(v1 == 1)
	SEQAN_TASSERT(numVertices(g) == 2)
	SEQAN_TASSERT(targetVertex(g, e2) == 1)
	SEQAN_TASSERT(sourceVertex(g, e2) == 0)
	SEQAN_TASSERT(numEdges(g) == 2)
	SEQAN_TASSERT(outDegree(g, v0) == 2)	
	SEQAN_TASSERT(inDegree(g, 1) == 1)
	SEQAN_TASSERT(inDegree(g, 0) == 1)	
	SEQAN_TASSERT(degree(g, 0) == 3)

	// Add more vertices and edges
	addVertex(g);  //2
	TVertexDescriptor v3 = addVertex(g);  //3
	addVertex(g);  //4
	addEdge(g,3,4,'g');
	TEdgeDescriptor my_edge = addEdge(g,3,1,'c');
	SEQAN_TASSERT(_getId(my_edge) == 3)
	addEdge(g,3,0,'t');
	SEQAN_TASSERT(v3 == 3)
	SEQAN_TASSERT(numVertices(g) == 5)
	SEQAN_TASSERT(targetVertex(g, e2) == 1)
	SEQAN_TASSERT(sourceVertex(g, e2) == 0)
	SEQAN_TASSERT(targetVertex(g, my_edge) == 1)
	SEQAN_TASSERT(sourceVertex(g, my_edge) == 3)
	SEQAN_TASSERT(numEdges(g) == 5)
	SEQAN_TASSERT(outDegree(g, v3) == 3)	

	// Output
	std::cout << g << std::endl;

	// Remove edges
	removeEdge(g,3,1,'c');
	removeEdge(g,0,1,'g');
	SEQAN_TASSERT(numEdges(g) == 3)

	// Remove vertices 
	TEdgeDescriptor e3 = addEdge(g,3,3,'a');
	addEdge(g,1,3,'a');
	addEdge(g,0,3,'c');
	addEdge(g,0,4,'t');
	SEQAN_TASSERT(outDegree(g, 0) == 3)
	SEQAN_TASSERT(outDegree(g, 1) == 1)
	SEQAN_TASSERT(targetVertex(g, e3) == 3)
	SEQAN_TASSERT(sourceVertex(g, e3) == 3)
	SEQAN_TASSERT(numEdges(g) == 7)
	removeVertex(g, v3);
	SEQAN_TASSERT(outDegree(g, 0) == 2)
	SEQAN_TASSERT(outDegree(g, 1) == 0)
	SEQAN_TASSERT(numVertices(g) == 4)
	SEQAN_TASSERT(numEdges(g) == 2)

	// Clear graph
	clearEdges(g);
	SEQAN_TASSERT(numVertices(g) == 4)
	SEQAN_TASSERT(numEdges(g) == 0)
	addEdge(g,2,0,'a');
	addEdge(g,4,1,'c');
	clearVertices(g);
	SEQAN_TASSERT(numVertices(g) == 0)
	SEQAN_TASSERT(numEdges(g) == 0)
	addVertex(g);addVertex(g);addVertex(g);
	addVertex(g);addVertex(g);
	addEdge(g,2,0,'t');
	addEdge(g,4,1,'g');
	clear(g);
	SEQAN_TASSERT(numVertices(g) == 0)
	SEQAN_TASSERT(numEdges(g) == 0)
	addVertex(g);addVertex(g);addVertex(g);
	addVertex(g);addVertex(g);
	addEdge(g,2,0,'c');
	addEdge(g,4,1,'g');
	addEdge(g,4,2,'t');
	removeVertex(g,3);
	SEQAN_TASSERT(numVertices(g) == 4)
	SEQAN_TASSERT(numEdges(g) == 3)
	SEQAN_TASSERT(outDegree(g, 4) == 2)
	SEQAN_TASSERT(inDegree(g, 4) == 0)

	//Transposes the graph in-place
	transpose(g); 
	SEQAN_TASSERT(numVertices(g) == 4)
	SEQAN_TASSERT(numEdges(g) == 3)
	SEQAN_TASSERT(outDegree(g, 4) == 0)
	SEQAN_TASSERT(inDegree(g, 4) == 2)
	StandardAutomaton g_copy(g);
	SEQAN_TASSERT(numVertices(g_copy) == 4)
	SEQAN_TASSERT(numEdges(g_copy) == 3)
	SEQAN_TASSERT(outDegree(g_copy, 4) == 0)
	SEQAN_TASSERT(inDegree(g_copy, 4) == 2)
	addVertex(g_copy);
	addEdge(g_copy,3,0,'a');
	g_copy = g;
	SEQAN_TASSERT(numVertices(g_copy) == 4)
	SEQAN_TASSERT(numEdges(g_copy) == 3)
	SEQAN_TASSERT(outDegree(g_copy, 4) == 0)
	SEQAN_TASSERT(inDegree(g_copy, 4) == 2)
	//Copies the graph and transposes just the copy
	transpose(g,g_copy);  // g does not change!
	SEQAN_TASSERT(numVertices(g_copy) == 4)
	SEQAN_TASSERT(numEdges(g_copy) == 3)
	SEQAN_TASSERT(outDegree(g_copy, 4) == 2)
	SEQAN_TASSERT(inDegree(g_copy, 4) == 0)
	removeVertex(g,0);

	// Adjacency matrix
	String<unsigned int> mat;
	getAdjacencyMatrix(g, mat);
	unsigned int len = (unsigned int) std::sqrt((double) length(mat));
	SEQAN_TASSERT(getValue(mat,1*len+4) == 1)
	SEQAN_TASSERT(getValue(mat,2*len+4) == 1)
	SEQAN_TASSERT(getValue(mat,0*len+2) == 0)

	// Test iterators
	typedef Iterator<StandardAutomaton, VertexIterator>::Type TVertexIterator;
	TVertexIterator itVert(g);
	SEQAN_TASSERT(getValue(itVert) == 1)
	++itVert;
	SEQAN_TASSERT(getValue(itVert) == 2)
	itVert++;
	SEQAN_TASSERT(getValue(itVert) == 4)
	goNext(itVert);
	SEQAN_TASSERT(atEnd(itVert) == true)

	addEdge(g,1,2,'T');
	typedef Iterator<StandardAutomaton, OutEdgeIterator>::Type TOutEdgeIterator;
	TOutEdgeIterator itEdge(g,1);
	// Slow
	SEQAN_TASSERT(sourceVertex(g, getValue(itEdge))==1)
	SEQAN_TASSERT(targetVertex(g, getValue(itEdge))==4)
	// Fast
	SEQAN_TASSERT(sourceVertex(itEdge)==1)
	SEQAN_TASSERT(targetVertex(itEdge)==4)
	
	SEQAN_TASSERT(sourceVertex(g, value(itEdge))==1)
	SEQAN_TASSERT(targetVertex(g, *itEdge)==4)
	SEQAN_TASSERT(atEnd(itEdge)==false)
	SEQAN_TASSERT(atBegin(itEdge)==true)
	goNext(itEdge);
	SEQAN_TASSERT(sourceVertex(itEdge)==1)
	SEQAN_TASSERT(targetVertex(itEdge)==2)
	++itEdge;
	itEdge++;
	SEQAN_TASSERT(atEnd(itEdge)==true)
	SEQAN_TASSERT(atBegin(itEdge)==false)
	goPrevious(itEdge);
	SEQAN_TASSERT(sourceVertex(g, getValue(itEdge))==1)
	SEQAN_TASSERT(targetVertex(g, getValue(itEdge))==2)
	--itEdge;
	itEdge++; 
	SEQAN_TASSERT(sourceVertex(g, getValue(itEdge))==1)
	SEQAN_TASSERT(targetVertex(g, getValue(itEdge))==2)
	itEdge--;
	itEdge--;
	SEQAN_TASSERT(atBegin(itEdge)==true)
	TOutEdgeIterator itEdge2(itEdge);
	TOutEdgeIterator itEdge3;
	itEdge3 = itEdge;
	SEQAN_TASSERT(itEdge == itEdge2)
	SEQAN_TASSERT(itEdge2 == itEdge3)
	goEnd(itEdge);
	SEQAN_TASSERT(itEdge2 != itEdge)
	goEnd(itEdge2);
	SEQAN_TASSERT(itEdge2 == itEdge)
	goBegin(itEdge2);
	SEQAN_TASSERT(atBegin(itEdge2)==true)
	SEQAN_TASSERT(&g == &hostGraph(itEdge))

	typedef Iterator<StandardAutomaton, EdgeIterator>::Type TEdgeIterator;
	TEdgeIterator itEd(g);
	SEQAN_TASSERT(sourceVertex(g, getValue(itEd))==1)
	SEQAN_TASSERT(targetVertex(g, getValue(itEd))==4)
	SEQAN_TASSERT(sourceVertex(g, value(itEd))==1)
	SEQAN_TASSERT(targetVertex(g, *itEd)==4)
	SEQAN_TASSERT(atEnd(itEd)==false)
	SEQAN_TASSERT(atBegin(itEd)==true)
	goNext(itEd);
	SEQAN_TASSERT(sourceVertex(g, getValue(itEd))==1)
	SEQAN_TASSERT(targetVertex(g, getValue(itEd))==2)
	SEQAN_TASSERT(atEnd(itEd)==false)
	SEQAN_TASSERT(atBegin(itEd)==false)
	++itEd;
	SEQAN_TASSERT(atEnd(itEd)==false)
	SEQAN_TASSERT(atBegin(itEd)==false)
	// Slow
	SEQAN_TASSERT(sourceVertex(g, getValue(itEd))==2)
	SEQAN_TASSERT(targetVertex(g, getValue(itEd))==4)
	// Fast
	SEQAN_TASSERT(sourceVertex(itEd)==2)
	SEQAN_TASSERT(targetVertex(itEd)==4)
	itEd++;
	itEd++;
	SEQAN_TASSERT(atEnd(itEd)==true)
	SEQAN_TASSERT(atBegin(itEd)==false)
	goPrevious(itEd);
	SEQAN_TASSERT(sourceVertex(g, getValue(itEd))==2)
	SEQAN_TASSERT(targetVertex(g, getValue(itEd))==4)
	--itEd;
	SEQAN_TASSERT(sourceVertex(g, getValue(itEd))==1)
	SEQAN_TASSERT(targetVertex(g, getValue(itEd))==2)
	TEdgeIterator itEd2(g);
	TEdgeIterator itEd3;
	goBegin(itEd);
	itEd3 = itEd;
	SEQAN_TASSERT(itEd == itEd2)
	SEQAN_TASSERT(itEd2 == itEd3)
	goEnd(itEd);
	SEQAN_TASSERT(itEd2 != itEd)
	goEnd(itEd2);
	SEQAN_TASSERT(itEd2 == itEd)
	goBegin(itEd2);
	SEQAN_TASSERT(itEd2 != itEd)
	SEQAN_TASSERT(&hostGraph(itEd) == &g)

	typedef Iterator<StandardAutomaton, AdjacencyIterator>::Type TAdjacencyIterator;
	TAdjacencyIterator itAd(g,1);
	SEQAN_TASSERT(getValue(itAd) == 4)
	SEQAN_TASSERT(&hostGraph(itAd) == &g)
	SEQAN_TASSERT(value(itAd) == 4)
	SEQAN_TASSERT(*itAd == 4)
	SEQAN_TASSERT(atEnd(itAd)==false)
	SEQAN_TASSERT(atBegin(itAd)==true)
	goNext(itAd);
	SEQAN_TASSERT(getValue(itAd)==2)
	SEQAN_TASSERT(atEnd(itAd)==false)
	SEQAN_TASSERT(atBegin(itAd)==false)
	++itAd;
	SEQAN_TASSERT(atEnd(itAd)==true)
	SEQAN_TASSERT(atBegin(itAd)==false)
	goBegin(itAd);
	itAd++;
	itAd++;
	itAd++;
	SEQAN_TASSERT(atEnd(itAd)==true)
	SEQAN_TASSERT(atBegin(itAd)==false)
	goPrevious(itAd);
	SEQAN_TASSERT(getValue(itAd)==2)
	--itAd;
	SEQAN_TASSERT(getValue(itAd)==4)
	SEQAN_TASSERT(atEnd(itAd)==false)
	SEQAN_TASSERT(atBegin(itAd)==true)
	goEnd(itAd);
	itAd--;
	SEQAN_TASSERT(getValue(itAd)==2)
	goBegin(itAd);
	TAdjacencyIterator itAd2(itAd);
	TAdjacencyIterator itAd3;
	itAd3 = itAd;
	SEQAN_TASSERT(itAd == itAd2)
	SEQAN_TASSERT(itAd2 == itAd3)
	goEnd(itAd);
	SEQAN_TASSERT(itAd2 != itAd)
	goEnd(itAd2);
	SEQAN_TASSERT(itAd2 == itAd)
	goBegin(itAd2);
	SEQAN_TASSERT(itAd2 != itAd)



//____________________________________________________________________________
// Automaton - Different alphabet
	typedef VertexDescriptor<Graph<Automaton<char> > >::Type VertexDescriptorType;
	typedef EdgeDescriptor<Graph<Automaton<char> > >::Type EdgeDescriptorType;
	Graph<Automaton<char> > automaton;
	VertexDescriptorType rootVertex = addVertex(automaton); // A = 0
	addVertex(automaton); // B = 1
	addVertex(automaton); // C = 2
	addVertex(automaton); // D = 3
	addVertex(automaton); // E = 4
	addVertex(automaton); // F = 5
	addEdge(automaton,0,1,'2');
	addEdge(automaton,1,0,'1');
	addEdge(automaton,4,0,'6');
	addEdge(automaton,0,3,'7');
	addEdge(automaton,1,1,'3');
	addEdge(automaton,1,2,'4');
	addEdge(automaton,5,1,'8');
	addEdge(automaton,2,5,'5');
	addEdge(automaton,3,4,'2');
	addEdge(automaton,5,3,'7');

	VertexDescriptorType succ;
	succ = getSuccessor(automaton,rootVertex,'7');
	SEQAN_TASSERT(succ == 3)
	// Throws an error in debug mode because edge does not exist
	//succ = getSuccessor(automaton,rootVertex,'6');
	succ = getSuccessor(automaton,succ,'2');
	SEQAN_TASSERT(succ == 4)
	succ = getSuccessor(automaton,succ,'6');
	SEQAN_TASSERT(succ == 0)
	// If no map is specified it is assumed that an edge cargo exists!!!
	succ = getSuccessor(automaton,succ,'2');
	SEQAN_TASSERT(succ == 1)

	// Now using shortcuts
	SEQAN_TASSERT(canParseString(automaton, rootVertex, "7262"))
	SEQAN_TASSERT(!canParseString(automaton, rootVertex, "726C"))
	SEQAN_TASSERT(canParseString(automaton, "7262"))
	SEQAN_TASSERT(!canParseString(automaton, "726C"))
	succ = parseString(automaton,rootVertex,"7262");
	SEQAN_TASSERT(succ == 1)
	std::string str = "7262";
	succ = parseString(automaton,rootVertex, str.begin(), str.end());
	SEQAN_TASSERT(succ == 1)
	String<char> str2("7262");
	succ = parseString(automaton,rootVertex, begin(str2), end(str2));
	SEQAN_TASSERT(succ == 1)
	String<char> input("7262");
	SEQAN_TASSERT(canParseString(automaton, rootVertex, input))
	SEQAN_TASSERT(canParseString(automaton, input))
	succ = parseString(automaton,rootVertex, input);
	SEQAN_TASSERT(succ == 1)

	// Additional cargo
	typedef Graph<Automaton<Dna, short> > TGraph9;
	typedef VertexDescriptor<TGraph9>::Type TVertexDescriptor9;
	typedef EdgeDescriptor<TGraph9>::Type TEdgeDescriptor9;
	typedef Size<TGraph9>::Type TSize9;

	TGraph9 g9;
	addVertex(g9);
	addVertex(g9);
	Dna aDna('a');
	Dna gDna('g');
	TEdgeDescriptor9 edg1 = addEdge(g9,0,1,aDna,12);
	TEdgeDescriptor9 edg2 = addEdge(g9,1,0,gDna,21);
	TGraph9 g10;
	transpose(g9, g10);
	TEdgeDescriptor9 edg1_10 = findEdge(g10, 1, aDna);
	TEdgeDescriptor9 edg1_11 = findEdge(g10, 0, gDna);
	SEQAN_TASSERT(getCargo(edg1)==12)
	SEQAN_TASSERT(getCargo(edg2)==21)
	SEQAN_TASSERT(getCargo(edg1_10)==12)
	SEQAN_TASSERT(getCargo(edg1_11)==21)

	// Multigraph
	StandardAutomaton multiG;
	addVertex(multiG);addVertex(multiG);addVertex(multiG);
	addEdge(multiG, 1, 2, 'a');
	TEdgeDescriptor edgeD2 = addEdge(multiG, 1, 2, 'c');
	addEdge(multiG, 1, 2, 'g');
	removeEdge(multiG, edgeD2);
	SEQAN_TASSERT(numEdges(multiG) == 2)
}

//////////////////////////////////////////////////////////////////////////////

void Test_WordGraph() {
//____________________________________________________________________________
// Standard automaton: No edge cargo

	typedef Graph<Automaton<Dna, String<Dna>, WordGraph<> > > TWordGraph;
	typedef VertexDescriptor<TWordGraph>::Type TVertexDescriptor;
	typedef EdgeDescriptor<TWordGraph>::Type TEdgeDescriptor;
	

	TWordGraph g;
	SEQAN_TASSERT(numVertices(g) == 0)
	SEQAN_TASSERT(numEdges(g) == 0)
	SEQAN_TASSERT(empty(g) == true)

	// Add vertex
	TVertexDescriptor v0 = addVertex(g);
	SEQAN_TASSERT(v0 == 0)
	SEQAN_TASSERT(outDegree(g, v0) == 0)	
	SEQAN_TASSERT(inDegree(g, 0) == 0)
	SEQAN_TASSERT(degree(g, 0) == 0)
	SEQAN_TASSERT(numVertices(g) == 1)
	SEQAN_TASSERT(empty(g) == false)
	addVertex(g);
	addVertex(g);
	TVertexDescriptor v3 = addVertex(g);
	SEQAN_TASSERT(isRoot(g, 0) == true)
	SEQAN_TASSERT(getRoot(g) == 0)
	assignRoot(g,3);
	SEQAN_TASSERT(getRoot(g) == 3)
	SEQAN_TASSERT(isRoot(g, 0) == false)
	SEQAN_TASSERT(isRoot(g, 3) == true)
	root(g) = 2;
	SEQAN_TASSERT(getRoot(g) == 2)
	SEQAN_TASSERT(isRoot(g, 3) == false)
	SEQAN_TASSERT(isRoot(g, 2) == true)

	// Add edge
	TEdgeDescriptor e1 =addEdge(g,v0,v3,"ag");
	SEQAN_TASSERT(findEdge(g,v0,'a') == e1)
	SEQAN_TASSERT(_getId(e1) == 0)
	// First letter -> edge label, all other letters into the cargo
	SEQAN_TASSERT(getCargo(e1) == "g")
	SEQAN_TASSERT(targetVertex(g, e1) == 3)
	SEQAN_TASSERT(sourceVertex(g, e1) == 0) 
	SEQAN_TASSERT(numEdges(g) == 1)
	SEQAN_TASSERT(outDegree(g, v0) == 1)	
	SEQAN_TASSERT(inDegree(g, v0) == 0)
	SEQAN_TASSERT(degree(g, v0) == 1)

	// Add further edges and vertices
	addVertex(g);
	TVertexDescriptor v5 = addVertex(g);
	TEdgeDescriptor e2 =addEdge(g,0,5,"g");
	SEQAN_TASSERT(_getId(e2) == 1) 
	SEQAN_TASSERT(v5 == 5)
	SEQAN_TASSERT(numVertices(g) == 6)
	SEQAN_TASSERT(targetVertex(g, e2) == 5)
	SEQAN_TASSERT(sourceVertex(g, e2) == 0)
	SEQAN_TASSERT(numEdges(g) == 2)
	removeEdge(g,0,5,String<Dna>("g"));
	SEQAN_TASSERT(numEdges(g) == 1)
	e2 =addEdge(g,0,5,"g");
	SEQAN_TASSERT(outDegree(g, v0) == 2)	
	SEQAN_TASSERT(inDegree(g, 5) == 1)
	SEQAN_TASSERT(degree(g, 0) == 2)
	SEQAN_TASSERT(getSuccessor(g, 0, "g") == 5)
	SEQAN_TASSERT(getSuccessor(g, 0, String<Dna>("ag")) == 3)  // The whole edge label or just the first letter
	SEQAN_TASSERT(getSuccessor(g, 0, "a") == getNil<TVertexDescriptor>())
	addVertex(g);
	addVertex(g);
	addEdge(g,3,1,"aggg");
	addEdge(g,3,4,"gg");
	addEdge(g,5,2,"aggg");
	addEdge(g,5,7,"g");
	addEdge(g,7,6,"g");
	SEQAN_TASSERT(parseString(g, 0, "agaggg") == 1)
	SEQAN_TASSERT(parseString(g, 0, "aga") == 3)  // Does not reach 1
	SEQAN_TASSERT(parseString(g, 0, 'g') == 5)
	SEQAN_TASSERT(parseString(g, 0, "ggg") == 6)
	SEQAN_TASSERT(parseString(g, 0, "gaggg") == 2)
	SEQAN_TASSERT(parseString(g, 0, "gagggg") == 2)
	assignRoot(g,0);

	// Output
	std::cout << g << std::endl;

	assignRoot(g,2);
	TWordGraph g_tmp(g);
	SEQAN_TASSERT(numVertices(g_tmp) == 8)
	SEQAN_TASSERT(parseString(g_tmp, 0, "agaggg") == 1)
	SEQAN_TASSERT(inDegree(g_tmp, 5) == 1)
	SEQAN_TASSERT(degree(g_tmp, 0) == 2)
	SEQAN_TASSERT(isRoot(g_tmp, 2) == true)
	TWordGraph g_assign;
	g_assign = g;
	SEQAN_TASSERT(numVertices(g_assign) == 8)
	SEQAN_TASSERT(parseString(g_assign, 0, "agaggg") == 1)
	SEQAN_TASSERT(inDegree(g_assign, 5) == 1)
	SEQAN_TASSERT(degree(g_assign, 0) == 2)

	// Transpose
	transpose(g, g_tmp);
	SEQAN_TASSERT(numVertices(g_tmp) == 8)
	SEQAN_TASSERT(parseString(g_tmp, 2, "aggg") == 5)
	SEQAN_TASSERT(inDegree(g_tmp, 5) == 2)
	SEQAN_TASSERT(outDegree(g_tmp, 0) == 0)
	SEQAN_TASSERT(isRoot(g_tmp, 2) == true)
}


//////////////////////////////////////////////////////////////////////////////

void Test_Tree() {
//____________________________________________________________________________
// Tree without edge cargo

	typedef Graph<Tree<void> > TTree;
	typedef VertexDescriptor<TTree>::Type TVertexDescriptor;
	typedef EdgeDescriptor<TTree>::Type TEdgeDescriptor;
	
	TTree g;
	SEQAN_TASSERT(empty(g) == true)
	createRoot(g);
	TVertexDescriptor rootV = getRoot(g);
	SEQAN_TASSERT(rootV == 0)
	SEQAN_TASSERT(isRoot(g, rootV) == true)
	SEQAN_TASSERT(root(g) == rootV)
	SEQAN_TASSERT(empty(g) == false)
	TVertexDescriptor childC1 = addChild(g,rootV);
	String<TVertexDescriptor> leaves;
	collectLeaves(g, rootV, leaves);
	TEdgeDescriptor childC1e = findEdge(g, rootV, childC1);
	SEQAN_TASSERT(_getVertexString(g)[0] == childC1e)
	SEQAN_TASSERT(getIdUpperBound(_getVertexIdManager(g)) == 2)
	SEQAN_TASSERT(getIdUpperBound(_getEdgeIdManager(g)) == 2)
	SEQAN_TASSERT(targetVertex(g, childC1e) == childC1) // Target in a tree = child
	SEQAN_TASSERT(sourceVertex(g, childC1e) == rootV)  // Source in a tree = parent
	SEQAN_TASSERT(childVertex(g, childC1e) == childC1)  // Shortcuts
	SEQAN_TASSERT(parentVertex(g, childC1e) == rootV)
	SEQAN_TASSERT(parentVertex(g, childC1) == rootV)
	SEQAN_TASSERT(empty(g) == false)
	SEQAN_TASSERT(outDegree(g, rootV) == 1)
	TVertexDescriptor childC2 = addChild(g,rootV);
	TVertexDescriptor childC3 = addChild(g,rootV);
	clear(leaves);
	collectLeaves(g, rootV, leaves);
	SEQAN_TASSERT(length(leaves) == 3)
	SEQAN_TASSERT(outDegree(g, rootV) == 3)
	SEQAN_TASSERT(childC1 == 1)
	SEQAN_TASSERT(childC2 == 2)
	SEQAN_TASSERT(childC3 == 3)
	TVertexDescriptor childC2C1 = addChild(g,childC2);
	TVertexDescriptor childC2C1C1 = addChild(g,childC2C1);
	TVertexDescriptor childC2C1C1C1 = addChild(g,childC2C1C1);
	TVertexDescriptor childC2C1C1C2 = addChild(g,childC2C1C1);
	TVertexDescriptor childC4 = addChild(g,rootV);
	SEQAN_TASSERT(inDegree(g, childC2C1) == 1) 
	SEQAN_TASSERT(outDegree(g, childC2C1) == 1)
	SEQAN_TASSERT(degree(g, childC2C1) == 2)
	SEQAN_TASSERT(numEdges(g) == 8)
	SEQAN_TASSERT(numVertices(g) == 9)
	SEQAN_TASSERT(numTreeEdges(g) == numVertices(g) - 1)
	TEdgeDescriptor childC2C1C1e = findEdge(g, childC2C1C1, childC2C1);
	
	// Output
	std::cout << g << std::endl;

	SEQAN_TASSERT(g.data_parent[0]==getNil<TVertexDescriptor>())
	SEQAN_TASSERT(g.data_parent[1]==0)
	SEQAN_TASSERT(g.data_parent[2]==0)
	SEQAN_TASSERT(g.data_parent[3]==0)
	SEQAN_TASSERT(g.data_parent[4]==2)
	SEQAN_TASSERT(g.data_parent[5]==4)
	SEQAN_TASSERT(g.data_parent[6]==5)
	SEQAN_TASSERT(g.data_parent[7]==5)
	SEQAN_TASSERT(g.data_parent[8]==0)
	_rebuildParentMap(g);
	SEQAN_TASSERT(g.data_parent[0]==getNil<TVertexDescriptor>())
	SEQAN_TASSERT(g.data_parent[1]==0)
	SEQAN_TASSERT(g.data_parent[2]==0)
	SEQAN_TASSERT(g.data_parent[3]==0)
	SEQAN_TASSERT(g.data_parent[4]==2)
	SEQAN_TASSERT(g.data_parent[5]==4)
	SEQAN_TASSERT(g.data_parent[6]==5)
	SEQAN_TASSERT(g.data_parent[7]==5)
	SEQAN_TASSERT(g.data_parent[8]==0)
	SEQAN_TASSERT(childVertex(g, childC2C1C1e) == childC2C1C1)  
	SEQAN_TASSERT(parentVertex(g, childC2C1C1e) == childC2C1)
	removeChild(g, rootV, childC2);
	SEQAN_TASSERT(g.data_parent[0]==getNil<TVertexDescriptor>())
	SEQAN_TASSERT(g.data_parent[1]==0)
	SEQAN_TASSERT(g.data_parent[2]==getNil<TVertexDescriptor>())
	SEQAN_TASSERT(g.data_parent[3]==0)
	SEQAN_TASSERT(g.data_parent[4]==getNil<TVertexDescriptor>())
	SEQAN_TASSERT(g.data_parent[5]==getNil<TVertexDescriptor>())
	SEQAN_TASSERT(g.data_parent[6]==getNil<TVertexDescriptor>())
	SEQAN_TASSERT(g.data_parent[7]==getNil<TVertexDescriptor>())
	SEQAN_TASSERT(g.data_parent[8]==0)
	SEQAN_TASSERT(numEdges(g) == 3)
	SEQAN_TASSERT(numVertices(g) == 4)
	SEQAN_TASSERT(empty(g) == false)
	SEQAN_TASSERT(inDegree(g, rootV) == 0) 
	SEQAN_TASSERT(outDegree(g, rootV) == 3)
	SEQAN_TASSERT(degree(g, rootV) == 3)
	childC2 = addChild(g,rootV);
	childC2C1 = addChild(g,childC2);
	childC2C1C1 = addChild(g,childC2C1);
	childC2C1C1C1 = addChild(g,childC2C1C1);
	childC2C1C1C2 = addChild(g,childC2C1C1);
	removeAllChildren(g, rootV);
	SEQAN_TASSERT(empty(g) == false)
	SEQAN_TASSERT(numEdges(g) == 0)
	SEQAN_TASSERT(numTreeEdges(g) == 0)
	SEQAN_TASSERT(numVertices(g) == 1) // Just the root
	SEQAN_TASSERT(inDegree(g, rootV) == 0) 
	SEQAN_TASSERT(outDegree(g, rootV) == 0)
	SEQAN_TASSERT(degree(g, rootV) == 0)
	addChild(g,rootV);addChild(g,rootV);
	SEQAN_TASSERT(empty(g) == false)
	SEQAN_TASSERT(numEdges(g) == 2)
	clearEdges(g);
	SEQAN_TASSERT(numEdges(g) == 0)
	SEQAN_TASSERT(numVertices(g) == 3)
	SEQAN_TASSERT(empty(g) == false)
	addChild(g,rootV);addChild(g,rootV);
	clearVertices(g);
	SEQAN_TASSERT(empty(g) == true)
	SEQAN_TASSERT(numVertices(g) == 0)
	SEQAN_TASSERT(numEdges(g) == 0)
	createRoot(g);
	childC1 = addChild(g,rootV);
	SEQAN_TASSERT(empty(g) == false)
	SEQAN_TASSERT(numEdges(g) == 1)
	childC3 = addChild(g,rootV);
	childC2 = addChild(g,rootV);
	childC2C1 = addChild(g,childC2);
	childC2C1C1 = addChild(g,childC2C1);
	childC2C1C1C1 = addChild(g,childC2C1C1);
	childC2C1C1C2 = addChild(g,childC2C1C1);
	childC4 = addChild(g,rootV);
	String<unsigned int> mat; 	// Adjacency matrix
	getAdjacencyMatrix(g, mat);
	unsigned int len = (unsigned int) std::sqrt((double) length(mat));
	SEQAN_TASSERT(getValue(mat, 0*len+8) == 1)
	SEQAN_TASSERT(getValue(mat, 8*len+0) == 0)
	SEQAN_TASSERT(getValue(mat, 3*len+0) == 0)
	SEQAN_TASSERT(getValue(mat, 0*len+3) == 1)
	SEQAN_TASSERT(getValue(mat, 0*len+4) == 0)
	SEQAN_TASSERT(numEdges(g) == 8)
	SEQAN_TASSERT(numVertices(g) == 9)
	transpose(g); 
	SEQAN_TASSERT(numEdges(g) == 8)
	SEQAN_TASSERT(numVertices(g) == 9)
	TTree g_copy(g);
	SEQAN_TASSERT(numEdges(g) == 8)
	SEQAN_TASSERT(numVertices(g) == 9)
	clear(g_copy);
	g_copy = g;
	SEQAN_TASSERT(numEdges(g) == 8)
	SEQAN_TASSERT(numVertices(g) == 9)
	transpose(g,g_copy);  
	g = g_copy;
	_rebuildParentMap(g);
	SEQAN_TASSERT(numEdges(g) == 8)
	SEQAN_TASSERT(numVertices(g) == 9)
	SEQAN_TASSERT(g.data_parent[0]==getNil<TVertexDescriptor>())
	SEQAN_TASSERT(g.data_parent[1]==0)
	SEQAN_TASSERT(g.data_parent[2]==0)
	SEQAN_TASSERT(g.data_parent[3]==0)
	SEQAN_TASSERT(g.data_parent[4]==3)
	SEQAN_TASSERT(g.data_parent[5]==4)
	SEQAN_TASSERT(g.data_parent[6]==5)
	SEQAN_TASSERT(g.data_parent[7]==5)
	SEQAN_TASSERT(g.data_parent[8]==0)
	removeOutEdges(g,childC2C1C1);
	SEQAN_TASSERT(g.data_parent[0]==getNil<TVertexDescriptor>())
	SEQAN_TASSERT(g.data_parent[1]==0)
	SEQAN_TASSERT(g.data_parent[2]==0)
	SEQAN_TASSERT(g.data_parent[3]==0)
	SEQAN_TASSERT(g.data_parent[4]==3)
	SEQAN_TASSERT(g.data_parent[5]==4)
	SEQAN_TASSERT(g.data_parent[6]==getNil<TVertexDescriptor>())
	SEQAN_TASSERT(g.data_parent[7]==getNil<TVertexDescriptor>())
	SEQAN_TASSERT(g.data_parent[8]==0)
	SEQAN_TASSERT(numVertices(g) == 9)
	SEQAN_TASSERT(numEdges(g) == 6)
	removeVertex(g,childC2C1);
	SEQAN_TASSERT(g.data_parent[0]==getNil<TVertexDescriptor>())
	SEQAN_TASSERT(g.data_parent[1]==0)
	SEQAN_TASSERT(g.data_parent[2]==0)
	SEQAN_TASSERT(g.data_parent[3]==0)
	SEQAN_TASSERT(g.data_parent[5]==getNil<TVertexDescriptor>())
	SEQAN_TASSERT(g.data_parent[6]==getNil<TVertexDescriptor>())
	SEQAN_TASSERT(g.data_parent[7]==getNil<TVertexDescriptor>())
	SEQAN_TASSERT(g.data_parent[8]==0)
	SEQAN_TASSERT(numEdges(g) == 4)

	SEQAN_TASSERT(numVertices(g) == 8)
	removeInEdges(g,childC2);
	SEQAN_TASSERT(numEdges(g) == 3)
	SEQAN_TASSERT(numVertices(g) == 8)
	removeOutEdges(g,rootV);
	SEQAN_TASSERT(g.data_parent[0]==getNil<TVertexDescriptor>())
	SEQAN_TASSERT(g.data_parent[1]==getNil<TVertexDescriptor>())
	SEQAN_TASSERT(g.data_parent[2]==getNil<TVertexDescriptor>())
	SEQAN_TASSERT(g.data_parent[3]==getNil<TVertexDescriptor>())
	SEQAN_TASSERT(g.data_parent[5]==getNil<TVertexDescriptor>())
	SEQAN_TASSERT(g.data_parent[6]==getNil<TVertexDescriptor>())
	SEQAN_TASSERT(g.data_parent[7]==getNil<TVertexDescriptor>())
	SEQAN_TASSERT(g.data_parent[8]==getNil<TVertexDescriptor>())
	SEQAN_TASSERT(numVertices(g) == 8)
	SEQAN_TASSERT(numEdges(g) == 0)
	SEQAN_TASSERT(empty(g) == false) 
	addVertex(g);
	TEdgeDescriptor my_edge = addEdge(g,0,1);
	removeEdge(g,my_edge);

//____________________________________________________________________________
// Tree with cargo

	typedef Pair<char, int> TPair;
	typedef Tree<TPair> TEdges;
	typedef Graph<TEdges> TCargoGraph;
	typedef VertexDescriptor<TCargoGraph>::Type TVertexDescriptor2;
	typedef EdgeDescriptor<TCargoGraph>::Type TEdgeDescriptor2;

	TCargoGraph g2;
	createRoot(g2);
	SEQAN_TASSERT(numVertices(g2) == 1)
	SEQAN_TASSERT(numEdges(g2) == 0)
	TVertexDescriptor2 ver1 = addChild(g2, getRoot(g2), TPair('a',3));
	SEQAN_TASSERT(numChildren(g2, getRoot(g2)) == 1);
	SEQAN_TASSERT(ver1 == 1)
	SEQAN_TASSERT(numVertices(g2) == 2)
	TVertexDescriptor2 ver2 = addChild(g2, getRoot(g2));
	SEQAN_TASSERT(ver2 == 2)
	SEQAN_TASSERT(numVertices(g2) == 3)
	TEdgeDescriptor2 ed1 =findEdge(g2,getRoot(g2),ver1);
	TEdgeDescriptor2 ed2 =findEdge(g2,getRoot(g2),ver2);
	SEQAN_TASSERT((getCargo(ed1)).i1 == 'a')
	SEQAN_TASSERT((getCargo(ed1)).i2 == 3)
	SEQAN_TASSERT(targetVertex(g2, ed1) == ver1)
	SEQAN_TASSERT(sourceVertex(g2, ed1) == getRoot(g2))
	SEQAN_TASSERT(numEdges(g2) == 2)
	assignCargo(ed2, TPair('b',4));
	SEQAN_TASSERT((getCargo(ed1)).i1 == 'a')
	SEQAN_TASSERT((getCargo(ed1)).i2 == 3)
	SEQAN_TASSERT((getCargo(ed2)).i1 == 'b')
	SEQAN_TASSERT((getCargo(ed2)).i2 == 4)
	cargo(ed1) = TPair('c',1);
	SEQAN_TASSERT((getCargo(ed1)).i1 == 'c')
	SEQAN_TASSERT((getCargo(ed1)).i2 == 1)
	assignRoot(g2,1);
	SEQAN_TASSERT(getRoot(g2) == 1)

	//// File read
	//fstream strmKnut;
	//TTree gKnut;
	//strmKnut.open(TEST_PATH "my_tree2.dot", ios_base::in);
	//String<String<char> > nodeMap;
	//String<String<char> > edgeMap;
	//read(strmKnut, gKnut, nodeMap, edgeMap, DotDrawing());
	//strmKnut.close();

	//assignRoot(gKnut, 26);
	//typedef Iterator<TTree, DfsPreorder>::Type TDfsPreorder;
	//TDfsPreorder dfsIt(gKnut, 26);
	//for(;!atEnd(dfsIt);++dfsIt) {
	//	std::cout << *dfsIt << ": ";
	//	std::cout << getProperty(nodeMap, *dfsIt) << std::endl;
	//}
	//std::cout << gKnut << std::endl;
}



//////////////////////////////////////////////////////////////////////////////

void Test_Fragment() {
	// Test Fragment
	typedef String<char> TString;
	typedef StringSet<TString, Dependent<> > TStringSet;
	typedef	Id<TStringSet>::Type TId;
	typedef	Size<TStringSet>::Type TSize;
	
	TStringSet str;
	TString str0("annual");	assignValueById(str, str0);
	TString str1("anneal"); assignValueById(str, str1);

	// Fragment: SeqId1, Begin1, SeqId2, Begin2, Length of Fragment
	Fragment<> f(0,4,1,4,2);
	SEQAN_TASSERT(f.seqId1 == 0)
	SEQAN_TASSERT(f.begin1 == 4)
	SEQAN_TASSERT(f.seqId2 == 1)
	SEQAN_TASSERT(f.begin2 == 4)
	SEQAN_TASSERT(f.len == 2)
	SEQAN_TASSERT(fragmentBegin(f, 0) == 4)
	SEQAN_TASSERT(fragmentBegin(f, 1) == 4)
	SEQAN_TASSERT(fragmentLength(f, 0) == 2)
	SEQAN_TASSERT(fragmentLength(f, 1) == 2)
	SEQAN_TASSERT(sequenceId(f, 0) == 0)
	SEQAN_TASSERT(sequenceId(f, 1) == 1)
	SEQAN_TASSERT(label(f, str, 0) == "al")
	SEQAN_TASSERT(label(f, str, 1) == "al")
	TId id2;
	TSize pos2;
	getProjectedPosition(f, 0, 5, id2, pos2);
	SEQAN_TASSERT(pos2 == 5)
	SEQAN_TASSERT(id2 == 1)
	getProjectedPosition(f, 1, 5, id2, pos2);
	SEQAN_TASSERT(pos2 == 5)
	SEQAN_TASSERT(id2 == 0)

	// Reversable Fragment
	typedef Fragment<Size<Fragment<> >::Type, ExactReversableFragment<> > TRevFrag;
	TRevFrag fRev(0,4,1,4,2);
	SEQAN_TASSERT(fRev.seqId1 == 0)
	SEQAN_TASSERT(fRev.begin1 == 4)
	SEQAN_TASSERT(fRev.seqId2 == 1)
	SEQAN_TASSERT(fRev.begin2 == 4)
	SEQAN_TASSERT(fRev.len == 2)
	SEQAN_TASSERT(fragmentBegin(fRev, 0) == 4)
	SEQAN_TASSERT(fragmentBegin(fRev, 1) == 4)
	SEQAN_TASSERT(fragmentLength(fRev, 0) == 2)
	SEQAN_TASSERT(fragmentLength(fRev, 1) == 2)
	SEQAN_TASSERT(sequenceId(fRev, 0) == 0)
	SEQAN_TASSERT(sequenceId(fRev, 1) == 1)
	SEQAN_TASSERT(label(fRev, str, 0) == "al")
	SEQAN_TASSERT(label(fRev, str, 1) == "al")
	SEQAN_TASSERT(!isReversed(fRev))
	getProjectedPosition(fRev, 0, 5, id2, pos2);
	SEQAN_TASSERT(pos2 == 5)
	SEQAN_TASSERT(id2 == 1)
	getProjectedPosition(fRev, 1, 5, id2, pos2);
	SEQAN_TASSERT(pos2 == 5)
	SEQAN_TASSERT(id2 == 0)
	TRevFrag fRev2(0,4,1,4,2, true);
	SEQAN_TASSERT(fRev2.seqId1 == 0)
	SEQAN_TASSERT(fRev2.begin1 == 4)
	SEQAN_TASSERT(fRev2.seqId2 == 1)
	SEQAN_TASSERT(fRev2.begin2 == 4)
	SEQAN_TASSERT(fRev2.len == 2)
	SEQAN_TASSERT(fragmentBegin(fRev2, 0) == 4)
	SEQAN_TASSERT(fragmentBegin(fRev2, 1) == 4)
	SEQAN_TASSERT(fragmentLength(fRev2, 0) == 2)
	SEQAN_TASSERT(fragmentLength(fRev2, 1) == 2)
	SEQAN_TASSERT(sequenceId(fRev2, 0) == 0)
	SEQAN_TASSERT(sequenceId(fRev2, 1) == 1)
	SEQAN_TASSERT(label(fRev2, str, 0) == "al")
	SEQAN_TASSERT(label(fRev2, str, 1) == "al")
	SEQAN_TASSERT(isReversed(fRev2))
	getProjectedPosition(fRev2, 0, 5, id2, pos2);
	SEQAN_TASSERT(pos2 == 4)
	SEQAN_TASSERT(id2 == 1)
	getProjectedPosition(fRev2, 1, 5, id2, pos2);
	SEQAN_TASSERT(pos2 == 4)
	SEQAN_TASSERT(id2 == 0)

	// Fragment: SeqId1, Begin1, SeqId2, Begin2, Length of Fragment
	Fragment<> f2(0,0,1,4,1);
	SEQAN_TASSERT(f2.seqId1 == 0)
	SEQAN_TASSERT(f2.begin1 == 0)
	SEQAN_TASSERT(f2.seqId2 == 1)
	SEQAN_TASSERT(f2.begin2 == 4)
	SEQAN_TASSERT(f2.len == 1)
	SEQAN_TASSERT(fragmentBegin(f2, 0) == 0)
	SEQAN_TASSERT(fragmentBegin(f2, 1) == 4)
	SEQAN_TASSERT(fragmentLength(f2, 0) == 1)
	SEQAN_TASSERT(fragmentLength(f2, 1) == 1)
	SEQAN_TASSERT(fragmentLength(f2) == 1)
	SEQAN_TASSERT(label(f2, str, 0) == "a")
	SEQAN_TASSERT(label(f2, str, 1) == "a")
	getProjectedPosition(f2, 0, 0, id2, pos2);
	SEQAN_TASSERT(pos2 == 4)
	SEQAN_TASSERT(id2 == 1)
	getProjectedPosition(f2, 1, 4, id2, pos2);
	SEQAN_TASSERT(pos2 == 0)
	SEQAN_TASSERT(id2 == 0)
}

//////////////////////////////////////////////////////////////////////////////

void Test_Hmm() {
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

	// Create an empty HMM
	THmm hmm;
	SEQAN_TASSERT(numVertices(hmm) == 0)
	SEQAN_TASSERT(numEdges(hmm) == 0)
	SEQAN_TASSERT(empty(hmm) == true)
	clearEdges(hmm);
	clearVertices(hmm);
	SEQAN_TASSERT(numVertices(hmm) == 0)
	SEQAN_TASSERT(numEdges(hmm) == 0)
	SEQAN_TASSERT(empty(hmm) == true)

	// Add state1
	TVertexDescriptor state1 = addVertex(hmm);
	SEQAN_TASSERT(length(_getVertexString(hmm)) == 1)
	SEQAN_TASSERT(empty(hmm) == false)
	SEQAN_TASSERT(outDegree(hmm, state1) == 0)
	SEQAN_TASSERT(inDegree(hmm, state1) == 0)
	SEQAN_TASSERT(degree(hmm, state1) == 0)
	SEQAN_TASSERT(numVertices(hmm) == 1)
	SEQAN_TASSERT(numEdges(hmm) == 0)
	emissionProbability(hmm, state1, dnaA) = 0.2;
	SEQAN_TASSERT(getEmissionProbability(hmm, state1, dnaA) == 0.2)
	emissionProbability(hmm, state1, dnaC) = 0.2;
	emissionProbability(hmm, state1, dnaG) = 0.3;
	emissionProbability(hmm, state1, dnaT) = 0.3;
	SEQAN_TASSERT(getEmissionProbability(hmm, state1, dnaA) == 0.2)
	SEQAN_TASSERT(getEmissionProbability(hmm, state1, dnaG) == 0.3)

	// Add state2
	String<TProbability> emis;
	resize(emis, alph_size);
	value(emis, ordValue(dnaA)) = 0.5;
	value(emis, ordValue(dnaC)) = 0.5;
	value(emis, ordValue(dnaG)) = 0.0;
	value(emis, ordValue(dnaT)) = 0.0;
	TVertexDescriptor state2 = addVertex(hmm, emis);
	SEQAN_TASSERT(numVertices(hmm) == 2)
	SEQAN_TASSERT(numEdges(hmm) == 0)
	SEQAN_TASSERT(getEmissionProbability(hmm, state2, dnaC) == 0.5)

	// Add state3
	TVertexDescriptor state3 = addVertex(hmm, emis);
	assignEmissionProbability(hmm, state3, dnaA, 0.3);
	assignEmissionProbability(hmm, state3, dnaC, 0.3);
	assignEmissionProbability(hmm, state3, dnaG, 0.2);
	assignEmissionProbability(hmm, state3, dnaT, 0.2);
	SEQAN_TASSERT(numVertices(hmm) == 3)
	SEQAN_TASSERT(numEdges(hmm) == 0)
	SEQAN_TASSERT(getEmissionProbability(hmm, state3, dnaC) == 0.3)

	// Add edges (transitions)
	TEdgeDescriptor e = addEdge(hmm, state1, state1, 0.95);
	SEQAN_TASSERT(numEdges(hmm) == 1)
	removeEdge(hmm, e);
	SEQAN_TASSERT(numEdges(hmm) == 0)
	e = addEdge(hmm, state1, state1, 0.95);
	SEQAN_TASSERT(numEdges(hmm) == 1)
	removeOutEdges(hmm, state1);
	SEQAN_TASSERT(numEdges(hmm) == 0)
	e = addEdge(hmm, state1, state1, 0.95);
	removeEdge(hmm, sourceVertex(hmm, e), targetVertex(hmm, e));
	SEQAN_TASSERT(numEdges(hmm) == 0)
	e = addEdge(hmm, state1, state1, 0.95);
	removeInEdges(hmm, state1);
	SEQAN_TASSERT(numEdges(hmm) == 0)
	e = addEdge(hmm, state1, state1, 0.95);
	SEQAN_TASSERT(outDegree(hmm, state1) == 1)
	SEQAN_TASSERT(inDegree(hmm, state1) == 1)
	SEQAN_TASSERT(degree(hmm, state1) == 2)
	e = addEdge(hmm, state1, state3);
	assignTransitionProbability(hmm, e, 0.05);
	THmm hmm_tr = hmm;
	transpose(hmm_tr);
	SEQAN_TASSERT(getTransitionProbability(hmm_tr, state3, state1) == 0.05)
	clear(hmm_tr);
	transpose(hmm, hmm_tr);
	SEQAN_TASSERT(getTransitionProbability(hmm_tr, state3, state1) == 0.05)
	e = addEdge(hmm, state3, state3);
	transitionProbability(hmm, e) = 0.4;
	SEQAN_TASSERT(getTransitionProbability(hmm, e) == 0.4)
	e = addEdge(hmm, state3, state1);
	transitionProbability(hmm, state3, state1) = 0.1;
	e = addEdge(hmm, state2, state2);
	assignTransitionProbability(hmm, state2, state2, 1.0);
	SEQAN_TASSERT(numVertices(hmm) == 3)
	SEQAN_TASSERT(numEdges(hmm) == 5)
	SEQAN_TASSERT(getTransitionProbability(hmm, state3, state1) == 0.1)
	
	// Add begin and end state
	TVertexDescriptor begState = addVertex(hmm);
	TVertexDescriptor eState = addVertex(hmm);
	addEdge(hmm, begState, state1, 1.0);
	addEdge(hmm, state3, eState, 0.5);
	addEdge(hmm, eState, eState, 1.0);
	beginState(hmm) = state3;
	SEQAN_TASSERT(numVertices(hmm) == 5)
	SEQAN_TASSERT(getBeginState(hmm) == state3)
	assignBeginState(hmm, begState);
	SEQAN_TASSERT(getBeginState(hmm) == begState)
	endState(hmm) = state3;
	SEQAN_TASSERT(getEndState(hmm) == state3)
	assignEndState(hmm, eState);
	SEQAN_TASSERT(getEndState(hmm) == eState)
	
	// Output
	std::cout << hmm << std::endl;

	// Change model
	removeVertex(hmm, state2);
	THmm hmm_copy(hmm);
	SEQAN_TASSERT(numVertices(hmm_copy) == 4)
	SEQAN_TASSERT(getBeginState(hmm_copy) == begState)
	SEQAN_TASSERT(getEndState(hmm_copy) == eState)
	clear(hmm_copy);
	SEQAN_TASSERT(numVertices(hmm_copy) == 0)
	hmm_copy = hmm;
	SEQAN_TASSERT(numVertices(hmm_copy) == 4)
	SEQAN_TASSERT(getBeginState(hmm_copy) == begState)
	SEQAN_TASSERT(getEndState(hmm_copy) == eState)
	SEQAN_TASSERT(idCount(_getEdgeIdManager(hmm_copy)) == 7)

	// Test silent states
	TVertexDescriptor testState1 = addVertex(hmm, emis, true);
	TVertexDescriptor testState2 = addVertex(hmm, emis, false);
	SEQAN_TASSERT(isSilent(hmm, testState1) != isSilent(hmm, testState2))
	assignSilentStatus(hmm, testState1, false);
	SEQAN_TASSERT(isSilent(hmm, testState1) == isSilent(hmm, testState2))
	silentStatus(hmm, testState1) = true;
	SEQAN_TASSERT(isSilent(hmm, testState1) != isSilent(hmm, testState2))
}

//////////////////////////////////////////////////////////////////////////////

void Test_GraphTypes() {
	Test_Directed();	// Directed graphs
	Test_Undirected();	// Undirected graphs
	Test_Automaton();	// Automatons
	Test_WordGraph();	// Word Graph
	Test_Tree();		// Trees
	Test_Fragment();	// Fragment
	Test_Hmm();			// Hmm

	debug::verifyCheckpoints("projects/library/seqan/graph_types/graph_impl_directed.h");
	debug::verifyCheckpoints("projects/library/seqan/graph_types/graph_impl_undirected.h");
	debug::verifyCheckpoints("projects/library/seqan/graph_types/graph_impl_automaton.h");
	debug::verifyCheckpoints("projects/library/seqan/graph_types/graph_impl_wordgraph.h");
	debug::verifyCheckpoints("projects/library/seqan/graph_types/graph_impl_tree.h");
	debug::verifyCheckpoints("projects/library/seqan/graph_types/graph_impl_fragment.h");
	debug::verifyCheckpoints("projects/library/seqan/graph_types/graph_impl_hmm.h");
}


}

#endif


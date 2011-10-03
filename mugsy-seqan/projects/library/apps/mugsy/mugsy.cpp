
//This code has mixed conventions because it combines different
//original sources. I've done my best to use the Seqan conventions
//where I can, but there is quite a jumble between STL and Seqan data
//structures

#define SEQAN_PROFILE 
//#define SEQAN_PROFILE2 //more verbose. SEQAN_PROFILE must also be defined
//#define SEQAN_TEST
#define NDEBUG //define this to disable assert statements
//#define KEEPCHAINTMP

#define TIMING
#ifdef TIMING
#include <time.h>
time_t now;
time_t lasttime;
#endif

//#define DEBUGGING //SVA custom debugging
//#define DEBUGGING_GRAPH //SVA custom debugging
//#define DEBUGGING2 //SVA verbose custom debugging

//There is some overhead to capturing scoring info
//Undef to turn off reporting of SP scores
//#define SCORING 


#include <seqan/basic.h>
#include <seqan/graph_msa.h>
#include <seqan/graph_types.h>
#include <seqan/graph_align.h>
#include <seqan/modifier.h>
#include <seqan/refinement.h>

#include "rna_alphabet.h"
#include <seqan/modifier.h>
#include <seqan/misc/misc_cmdparser.h>

//#include "sangiuoli/mummer/trunk/MUMmer3.20/src/tigr/delta.hh"

#include <sstream>
#include <fstream>
#include <vector>
#include <set>
#include <queue>
#include <list>
#include <bitset>
#include <algorithm>

#include <cstdlib>
#include <errno.h>
//#include <stdio.h>
//#include <stdlib.h>
#include <libgen.h>

#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/strong_components.hpp>
#include <boost/graph/topological_sort.hpp>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>

#include <boost/graph/edmonds_karp_max_flow.hpp>
#include <boost/graph/kolmogorov_max_flow.hpp>
#include <boost/graph/push_relabel_max_flow.hpp>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/read_dimacs.hpp>
#include <boost/graph/graph_utility.hpp>

#ifdef SEQAN_PROFILE
SEQAN_PROTIMESTART(__myProfileTime); // Profiling
#endif

//Example transform in multiz-tba/trunk/transformcoords.cpp
#include "transformcoords.h"

using namespace seqan;
using namespace std;

struct s_offset{
  unsigned int offset;
  unsigned int spanlen;
  unsigned int seqlen;
  unsigned int orient;
};

struct iloc{
  int first;
  int second;
  int blocknum;
};


struct s_score{
  unsigned int numGapEx;
  unsigned int numGap;
  unsigned int numPairs;
  unsigned int numIdents;
  unsigned int alignLen;
  unsigned int totalLen;
  unsigned int alignScore;
  unsigned int seqCount;
  String<unsigned int> colCount;
  String<unsigned int> pairCount;
};

//////////////////////////////////////////////////////////////////////////////
namespace SEQAN_NAMESPACE_MAIN
{

  struct vectorsizecmp {
    bool operator()( const String<Fragment<> > & s1, const String<Fragment<> > & s2 ) const {
      return length(s1) > length(s2);
    }
  };

  template <typename TMap>
  class lcblencmp{
  public:
    lcblencmp(TMap & m)
      :myMap(&m)
    {}
    bool operator() ( const int i, const int j) const {
      assert(myMap != NULL);
      assert(myMap->find(i)!=myMap->end());
      assert(myMap->find(j)!=myMap->end());
      return (myMap->find(i)->second > myMap->find(j)->second);
    }
    TMap *myMap;
  };
  
  template<typename TGraph>
  class vertexdegreecmp
  {
  public:
    vertexdegreecmp(TGraph & g)
      :myGraph(&g)
    {}
    template<typename TVertexDescriptor>
    bool operator()( const TVertexDescriptor v1, const TVertexDescriptor v2 ) const {
      return degree(*myGraph,v1) > degree(*myGraph,v2);
    }
    TGraph * myGraph;
  };
  
  template<typename TPosScores>
  class edgeposscorecmp
  {
  public:
    edgeposscorecmp(TPosScores & p)
      :posscores(&p)
    {}
    template<typename TEdgeDescriptor>
    bool operator()( const TEdgeDescriptor &e1, const TEdgeDescriptor &e2 ) const {
      return posscores->find(e1)->second < posscores->find(e2)->second;
    }
    TPosScores * posscores;
  };
  
  template<typename TPos>
  class poscmp
  {
  public:
    poscmp()
    {}
    bool operator()( const TPos &e1, const TPos &e2 ) const {
      if(e1.first==e2.first){
	//return false if e2 is interval close
	/*
	if(e2.second == false){
	  return 0;
	}
	else{
	  return 1;
	}
	*/
	return e1.second < e2.second;
      }
      else{
	return e1.first < e2.first;
      }
    }
  };

  template<typename TGraph>
  class vertexposcmp
  {
  public:
    vertexposcmp(TGraph & g)
      :myGraph(&g)
    {}
    template<typename TVertexDescriptor>
    bool operator()( const TVertexDescriptor v1, const TVertexDescriptor v2 ) const {
      return fragmentBegin(*myGraph,v1) < fragmentBegin(*myGraph,v2);
    }
    TGraph * myGraph;
  };

  template<typename TScoreMap>
  class edgescorecmp
  {
  public:
    edgescorecmp(TScoreMap * s)
      :myScoreMap(s)
    {}
    edgescorecmp()
      :myScoreMap(NULL)
    {}
    template<typename TEdgeDescriptor>
    bool operator()( const TEdgeDescriptor e1, const TEdgeDescriptor e2 ) const {
      if(myScoreMap!=NULL && myScoreMap->size()>0){
	assert(myScoreMap->find(e1)!=myScoreMap->end());
	assert(myScoreMap->find(e2)!=myScoreMap->end());
	if(abs(cargo(e1)) == abs(cargo(e2))){
	  //Secondary sort on adjacency score
	  return myScoreMap->find(e1)->second < myScoreMap->find(e2)->second;
	}
	else{
	  //Primary sort on consistency score
	  return abs(cargo(e1))<abs(cargo(e2));
	}
      }
      else{
	//Primary sort on consistency score
	return abs(cargo(e1))<abs(cargo(e2));
      }
    }
    TScoreMap * myScoreMap;
  };

  template<typename TBlock, typename TSize=unsigned int>  
  class blockorder
  {
  public:
    blockorder(TSize s){
      currentSeq = s;
    }
    bool operator()( const TBlock * s1, const TBlock * s2 ) const {
      assert(s1->currentSeq == currentSeq);
      assert(s2->currentSeq == currentSeq);
      return s1->begCoord < s2->begCoord;
    }
    TSize currentSeq;
  };

  template<typename TComponent = unsigned int, typename TSize = unsigned int, 
	   typename TVertexDescriptor = unsigned int, typename TPos = unsigned int>
  class SVABlock
  {
  public:
    SVABlock()
    {}
    SVABlock(const SVABlock & s)
      :begCoord(s.begCoord),
       endCoord(s.endCoord),
       orient(s.orient),
       c(s.c),
       currentSeq(s.currentSeq),
       currV(s.currV)
    {}
    SVABlock(TComponent inc, TSize s, TPos b, TPos e, char o, TVertexDescriptor v)
      :begCoord(b),endCoord(e),orient(o),c(inc),currentSeq(s)
    {
      currV.push_back(v);
    }
    template<typename TGraph>
    void addVertex(TGraph & g, TVertexDescriptor v){
      assert(sequenceId(g,v)==currentSeq);
      if(fragmentBegin(g,v)<begCoord){
	begCoord=fragmentBegin(g,v);
      }
      if(fragmentBegin(g,v)+fragmentLength(g,v)>endCoord){
	endCoord=fragmentBegin(g,v)+fragmentLength(g,v);
      }
      currV.push_back(v);
    }
    
    TPos begCoord;
    TPos endCoord;
    char orient;//'+' or '-', consider changing to false, true to be consistent with fragment.reversed
    TComponent c;
    TSize currentSeq;
    std::vector<TVertexDescriptor> currV;
  };
}
void process_mem_usage(double& vm_usage, double& resident_set)
{
   using std::ios_base;
   using std::ifstream;
   using std::string;

   vm_usage     = 0.0;
   resident_set = 0.0;

   // 'file' stat seems to give the most reliable results
   //
   ifstream stat_stream("/proc/self/stat",ios_base::in);

   // dummy vars for leading entries in stat that we don't care about
   //
   string pid, comm, state, ppid, pgrp, session, tty_nr;
   string tpgid, flags, minflt, cminflt, majflt, cmajflt;
   string utime, stime, cutime, cstime, priority, nice;
   string O, itrealvalue, starttime;

   // the two fields we want
   //
   unsigned long vsize;
   long rss;

   stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
               >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
               >> utime >> stime >> cutime >> cstime >> priority >> nice
               >> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest

   long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
   vm_usage     = vsize / 1024.0;
   resident_set = rss * page_size_kb;
}


//////////////////////////////////////////////////////////////////////////////
// Connected components
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
//Return minimum distance seperating two intervals s1-e1 and s2-e2
//If intervals are overlapping distance is 0
template<typename TSize>
inline unsigned int getIntervalDist(TSize &s1, TSize &e1, TSize &s2, TSize &e2){
  //Overlapping or contained
  if(s1>s2 && s1<e2){
    return 0;
  }
  else{
    if(s2>s1 && s2<e1){
      return 0;
    }
    else{
      if(s1<s2){
	return s2-s1;
      }
      else{
	return s1-s2;
      }
    }
  }
}

template<typename TGraph, 
	 typename TVertexDescriptor,
	 typename TGenomeVertexMapIter>
inline unsigned int distance(TGraph const& g, 
			     TVertexDescriptor const & u,
			     std::pair<TGenomeVertexMapIter,TGenomeVertexMapIter> &vmapiter){
  typedef unsigned int TSize;
  TSize mindist = std::numeric_limits<TSize>::max();
  for (TGenomeVertexMapIter vit=vmapiter.first; vit!=vmapiter.second; ++vit){
    if(sequenceId(g,u)==sequenceId(g,vit->second)){
      TSize beg1 = fragmentBegin(g,u);
      TSize end1 = beg1+fragmentLength(g,u);
      TSize beg2 = fragmentBegin(g,vit->second);
      TSize end2 = beg2+fragmentLength(g,vit->second);
      TSize dist = getIntervalDist(beg1,end1,beg2,end2);
      mindist = dist<mindist ? dist : mindist;
    }
  }
  return mindist;
}


template<typename TSpec, 
	 typename TVertexDescriptor, 
	 typename TTokenMap, 
	 typename TComponents, 
	 typename TVal, 
	 typename TGenomeVertexMap, 
	 typename TNames,
	 typename TSize>
inline void
_cc_visit_g_ranked(Graph<TSpec> const& g,
		   TVertexDescriptor const u,
		   TTokenMap& tokenMap,
		   TComponents& components,
		   TVal label, 
		   TVal &maxlabel, //changed to reference so that i can reassign
		   std::vector<TGenomeVertexMap> & genomeMap,
		   TNames &genomeNames,
		   TSize &maxdist)
{
  //SEQAN_CHECKPOINT
	typedef typename Iterator<Graph<TSpec>, AdjacencyIterator>::Type TAdjacencyIterator; 
	typedef typename EdgeDescriptor<Graph<TSpec> >::Type TEdgeDescriptor; 
	typedef typename Iterator<Graph<TSpec>, OutEdgeIterator>::Type TOutEdgeIterator;
	//Add all edges from u to ccedges
	if(getProperty(tokenMap, u) == false){ 
	  //TODO support for genomeidx in addition to sequenceid
	  assert(sequenceId(g,u)<length(genomeNames));
#ifdef DEBUGGING
	  std::cout << "Connecting vertex " << u << " from sequence " << sequenceId(g,u) << std::endl;
#endif
	  assert(label<genomeMap.size());
	  //Multiple copies of this genome in the current component
	  TSize gname = genomeNames[sequenceId(g,u)];
	  if(genomeMap[label].find(genomeNames[sequenceId(g,u)]) != genomeMap[label].end()){
	    std::pair<typename TGenomeVertexMap::iterator,typename TGenomeVertexMap::iterator> vmapiter =genomeMap[label].equal_range(gname); 
	    unsigned int dist = distance(g,u,vmapiter);
#ifdef DEBUGGING
	    std::cout << "Multiple copies found for genome " << gname
		      << " seq " << sequenceId(g,u)  << " on vertex " << u << std::endl;
	    std::cout << "Minimum distance to a copy: " << dist << std::endl;
	    for(typename TGenomeVertexMap::iterator it=vmapiter.first;it!=vmapiter.second;it++){
	      std::cout << "Copy " << it->second  << " seq_id:" << sequenceId(g,it->second) << std::endl;
	      assert((TSize)genomeNames[sequenceId(g,it->second)]==gname);
	      assert(dist>=0);
	      assert(maxdist>0);
	      if(dist<=maxdist){
		assert(sequenceId(g,it->second)==sequenceId(g,u));
	      }
	    }
#endif
	    //Start a new component
	    if(dist>maxdist){
	      //Increment and set maxlabel
	      ++maxlabel;
	      label=maxlabel;
#ifdef DEBUGGING
	      std::cout << "Starting new component " << label << " for " << u << std::endl;
#endif
	      genomeMap.push_back(TGenomeVertexMap());
	      //Ensure genome is not present in the existing label
	      assert(label<genomeMap.size());
	      assert(label==genomeMap.size()-1);
	      assert(genomeMap[label].find(gname) == genomeMap[label].end());
	    }
	  }
#ifdef DEBUGGING
	  std::cout << "Component " << label << " V:" << u << " seq:" << sequenceId(g,u) << std::endl;
#endif
	  assert(label<genomeMap.size());
	  //Add vertex,u to component,label
	  genomeMap[label].insert(std::make_pair(gname,u));
	  assignProperty(tokenMap, u, true);
	  assignProperty(components, u, label);
	  
	  //Capture all edges for this vertex
	  std::vector<TEdgeDescriptor> ccedges;
	  for(TOutEdgeIterator itOut(g, u);!atEnd(itOut); ++itOut) {
	    //TODO if(!visited) shortcut
	    ccedges.push_back(*itOut);
	  }
	  assert(ccedges.size()==degree(g,u));
	  //Sort edges on consistency and visit most consistent edges first in a greedy fashion
	  sort(ccedges.begin(),ccedges.end(),edgescorecmp<std::map<TEdgeDescriptor,float> >());

	  for(typename std::vector<TEdgeDescriptor>::reverse_iterator cit = ccedges.rbegin();cit!=ccedges.rend();++cit){
	    TVertexDescriptor s = getSource(*cit);
	    TVertexDescriptor t = getTarget(*cit);
	    assert(s==u || t==u);
	    if(s!=u){
	      assert(getProperty(tokenMap,t)==true);
	      if (getProperty(tokenMap, s) == false) {
#ifdef DEBUGGING
		std::cout << " edge " << u << "-" << s;
		std::cout << std::endl;
#endif
		_cc_visit_g_ranked(g, s, tokenMap, components, label, maxlabel,genomeMap,genomeNames,maxdist);
	      }
	    }
	    else{
	      if(t!=u){
		assert(getProperty(tokenMap,s)==true);
		if (getProperty(tokenMap, t) == false) {
#ifdef DEBUGGING
		  std::cout << " edge " << u << "-" << t;
		  std::cout << std::endl;
#endif
		  _cc_visit_g_ranked(g, t, tokenMap, components, label, maxlabel,genomeMap,genomeNames,maxdist);
		}
	      }
	      else{
		assert(false);
	      }
	    }
	  }
	}
}
//////////////////////////////////////////////////////////////////////////////

/**
.Function.connected_components:
..cat:Graph
..summary:Decomposes an undirected graph into its connected components.
..signature:connected_components(g, components)
..param.g:In-parameter:An undirected graph.
...type:Spec.Undirected graph
..param.components:Out-parameter:A property map.
...remarks:Each vertex is mapped to a component id. If two vertices share the same id they are in the same component.
..returns: The number of components.
*/

template<typename TSpec, typename TComponents, typename TNames, typename TSize2>
typename Size<Graph<TSpec> >::Type
connected_components_by_genome_ranked_RECURSIVE(Graph<TSpec> const& g,
						TComponents& components,
						TNames &genomeNames,
						TSize2 maxdist)
{
  //SEQAN_CHECKPOINT

	typedef typename Size<Graph<TSpec> >::Type TSize;
	typedef typename Iterator<Graph<TSpec>, EdgeIterator>::Type TEdgeIterator;
	typedef typename Iterator<Graph<TSpec>, VertexIterator>::Type TVertexIterator;
	typedef typename VertexDescriptor<Graph<TSpec> >::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<Graph<TSpec> >::Type TEdgeDescriptor;

	typedef std::multimap<TSize,TVertexDescriptor> TGenomeVertexMap;

	typedef typename Iterator<Graph<TSpec>, AdjacencyIterator>::Type TAdjacencyIterator; 
	typedef typename Iterator<Graph<TSpec>, OutEdgeIterator>::Type TOutEdgeIterator;

	clear(components);
	resizeVertexMap(g,components);

#ifdef DEBUGGING	
	std::cout << "Calculating connected components on" << length(genomeNames) << "genomes" << std::endl;
#endif

	// Initialization
	String<bool> tokenMap;
	fill(tokenMap, getIdUpperBound(_getVertexIdManager(g)), false);
	
	// Genome tracker
	std::vector<TGenomeVertexMap> genomeMap(1);


	// Find connected components greedy on consistency score
	TSize label = 0;
	TSize maxlabel = label;
	TEdgeIterator itE(g);
	std::vector<TEdgeDescriptor> ccedges;
	std::set<TVertexDescriptor> visited;
	for(;!atEnd(itE);goNext(itE)){
	  ccedges.push_back(*itE);
	}
	//Sort on consistency score
	sort(ccedges.begin(),ccedges.end(),edgescorecmp<std::map<TEdgeDescriptor,float> >());
	TVertexDescriptor s,t;
	//From most consistent edges in G to least, determine CC
	for(typename std::vector<TEdgeDescriptor>::reverse_iterator cit = ccedges.rbegin();cit!=ccedges.rend();cit++){
	  s = getSource(*cit);
	  if (getProperty(tokenMap, s) == false) {
#ifdef DEBUGGING
	    std::cout << "Component" << maxlabel
		      << std::endl;
#endif
	    //Capture all vertices connected to s
	    _cc_visit_g_ranked(g, s, tokenMap, components, label, maxlabel, genomeMap,genomeNames,maxdist);
	    ++maxlabel;
	    label=maxlabel;
	    genomeMap.push_back(TGenomeVertexMap());
	  }
	  t = getTarget(*cit);
	  if (getProperty(tokenMap, t) == false) {
	    _cc_visit_g_ranked(g, t, tokenMap, components, label, maxlabel, genomeMap,genomeNames,maxdist);
	    ++maxlabel;
	    label=maxlabel;
	    genomeMap.push_back(TGenomeVertexMap());
	  }
	}
	//Capture all vertices with no edges degree==0
	TVertexIterator it(g);
	TVertexDescriptor u;
	for(;!atEnd(it);goNext(it)) {
	  u = getValue(it);
	  if (getProperty(tokenMap, u) == false) {
#ifdef DEBUGGING
	    std::cout << "Component" << maxlabel
		      << std::endl;
#endif
	    assert(degree(g,u)==0);
	    _cc_visit_g_ranked(g, u, tokenMap, components, label, maxlabel, genomeMap,genomeNames,maxdist);
	    ++maxlabel;
	    label=maxlabel;
	    genomeMap.push_back(TGenomeVertexMap());
	  }
	}
	return label;
}


//connected_components_by_genome_ranked()
//
//Connected components greedy on consistency score and ensuring one
//anchor per genome.
//
//Used to convert segment graph (V=genome segments on one genome) into
//anchor graph (V=genome segments on multiple genomes)
//
//Run DFS to determine connected components. Order traversal by edge
//score largest-smallest. Break and start a new component upon
//encountering a second anchor in a genome that has already been
//visited if the new anchor > maxdist from the other anchors already
//visited
template<typename TSpec, typename TComponents, typename TNames, typename TSize2>
typename Size<Graph<TSpec> >::Type
connected_components_by_genome_ranked(Graph<TSpec> const& g,
				      TComponents& components,
				      TNames &genomeNames,
				      TSize2 maxdist){
  //SEQAN_CHECKPOINT
  typedef typename Size<Graph<TSpec> >::Type TSize;
  typedef typename Iterator<Graph<TSpec>, EdgeIterator>::Type TEdgeIterator;
  typedef typename Iterator<Graph<TSpec>, VertexIterator>::Type TVertexIterator;
  typedef typename VertexDescriptor<Graph<TSpec> >::Type TVertexDescriptor;
  typedef typename EdgeDescriptor<Graph<TSpec> >::Type TEdgeDescriptor;
  
  typedef std::multimap<TSize,TVertexDescriptor> TGenomeVertexMap;
  
  typedef typename Iterator<Graph<TSpec>, AdjacencyIterator>::Type TAdjacencyIterator; 
  typedef typename Iterator<Graph<TSpec>, OutEdgeIterator>::Type TOutEdgeIterator;
  
  clear(components);
  resizeVertexMap(g,components);
  
#ifdef DEBUGGING	
  std::cout << "Calculating connected components" << std::endl;
#endif
  
  // Initialization
  String<bool> tokenMap;
  fill(tokenMap, getIdUpperBound(_getVertexIdManager(g)), false);
  fill(components,getIdUpperBound(_getVertexIdManager(g)), 0);
  for(unsigned int i=0;i<getIdUpperBound(_getVertexIdManager(g));++i){
    assignProperty(components,i,0);
    assert(getProperty(components,i)==0);
  }
  // Genome tracker
  std::vector<TGenomeVertexMap> genomeMap(1);
  
  //TODO
  //Initial CC with maxdist
  //Save all nodes,edges > maxdist 
  //Score edges by adjacency score
  //Break edges < cutoff || keep only highest scoring node per genome
  //Recompute CC
  
  // Connected components
  TEdgeIterator itE(g);
  std::vector<TEdgeDescriptor> ccedges;
  int maxlabel=-1;
  for(;!atEnd(itE);goNext(itE)){
    ccedges.push_back(*itE);
  }
  sort(ccedges.begin(),ccedges.end(),edgescorecmp<std::map<TEdgeDescriptor,float> >());
  std::vector<std::pair<TVertexDescriptor,TVertexDescriptor> > stack;
  //outer loop ensures we visit disconnected subgraphs, considering most consistent edges first
  for(typename std::vector<TEdgeDescriptor>::reverse_iterator cit = ccedges.rbegin();cit!=ccedges.rend();cit++){
    TVertexDescriptor s = getSource(*cit);
    TVertexDescriptor t = getTarget(*cit);
#ifdef DEBUGGING
    std::cout << "Edge score:" << cargo(*cit) << std::endl;
#endif
    if (getProperty(tokenMap, s) == false){
      if(getProperty(tokenMap,t) == false) {
	stack.push_back(std::make_pair(s, t));
	stack.push_back(std::make_pair(t, s));
      }
      else{
	stack.push_back(std::make_pair(s, t));
      }
    }
    else{
      if(getProperty(tokenMap,t) == false) {
	stack.push_back(std::make_pair(t, s));
      }
    }
    while(!stack.empty()){
      std::pair<TVertexDescriptor,TVertexDescriptor> & node = stack.back();
      TVertexDescriptor u = node.first;
      TVertexDescriptor prev = node.second;
      assert(sequenceId(g,u)<length(genomeNames));
      TSize gname = genomeNames[sequenceId(g,u)];
      stack.pop_back();
      if(getProperty(tokenMap,u)==false){
	assert(getProperty(components,u)==0);
#ifdef DEBUGGING
	std::cout << "New node " << u << " " << getProperty(components,u) 
		  << std::endl;
#endif
	//Encountered new node, assign label, track genome
	assignProperty(tokenMap, u, true);
	int label=-1;
	if(getProperty(tokenMap,prev)==false){
	  //Use new label
	  ++maxlabel;
#ifdef DEBUGGING
	  std::cout << "Starting new component " << maxlabel << " for " << u << std::endl;
#endif
	  label=maxlabel;
	  genomeMap.push_back(TGenomeVertexMap());
	}
	else{
	  //there is already a label
	  int prevlabel=getProperty(components,prev);
	  if(genomeMap[prevlabel].find(gname) != genomeMap[prevlabel].end()){
	    //there is already a genome, retrieve all the anchors in this genome to determine the anchors
	    std::pair<typename TGenomeVertexMap::iterator,typename TGenomeVertexMap::iterator> vmapiter =genomeMap[prevlabel].equal_range(gname); 
	    unsigned int dist = distance(g,u,vmapiter);
#ifdef DEBUGGING
	    std::cout << "Multiple copies found for genome " << gname
		      << " seq " << sequenceId(g,u)  << " on vertex " << u << std::endl;
	    std::cout << "Minimum distance to a copy: " << dist << std::endl;
	    for(typename TGenomeVertexMap::iterator it=vmapiter.first;it!=vmapiter.second;it++){
	      std::cout << "Copy " << it->second  << " seq_id:" << sequenceId(g,it->second) << std::endl;
	      assert((TSize)genomeNames[sequenceId(g,it->second)]==gname);
	      assert(dist>=0);
	      if(dist<=maxdist){
		assert(sequenceId(g,it->second)==sequenceId(g,u));
	      }
	    }
#endif
	    //Start a new component
	    if(dist>=maxdist){
	      //Increment and set maxlabel
	      ++maxlabel;
#ifdef DEBUGGING
	      std::cout << "Starting new component " << maxlabel << " for " << u << std::endl;
#endif
	      label=maxlabel;
	      genomeMap.push_back(TGenomeVertexMap());	      
	      //Ensure genome is not present in the existing label
	      assert((unsigned int)maxlabel<genomeMap.size());
	      assert(genomeMap[maxlabel].find(gname) == genomeMap[maxlabel].end());
	    }
	    else{
#ifdef DEBUGGING
	      std::cout << "Adding to component " << getProperty(components,prev) << " for " << u << std::endl;
#endif
	      label=prevlabel;
	    }
	  }
	  else{
#ifdef DEBUGGING
	      std::cout << "Adding to component " << getProperty(components,prev) << " for " << u << std::endl;
#endif
	      label=prevlabel;
	  }
	}

	assignProperty(components, u, label);
#ifdef DEBUGGING
	std::cout << "V:" << u << " component " << getProperty(components,u) << std::endl;
#endif
	//assert(getProperty(components,u)>=0);
	assert(getProperty(components,u) < genomeMap.size());

	genomeMap[label].insert(std::make_pair(gname,u));
      
	//Add all edges from u to ccedges
	std::vector<TEdgeDescriptor> ccedges;
	//Edge iterator
	TOutEdgeIterator itOut(g, u);
	for(;!atEnd(itOut); ++itOut) {
	  if(getProperty(tokenMap,getSource(*itOut)) || getProperty(tokenMap,getTarget(*itOut))){
	    ccedges.push_back(*itOut);
	  }
	}
	//assert(ccedges.size()==degree(g,u));
	sort(ccedges.begin(),ccedges.end(),edgescorecmp<std::map<TEdgeDescriptor,float> >());
	//traverse scores low->high so highest scores are last on the stack
	for(typename std::vector<TEdgeDescriptor>::iterator cit = ccedges.begin();cit!=ccedges.end();cit++){
	  TVertexDescriptor s = getSource(*cit);
	  TVertexDescriptor t = getTarget(*cit);
	  assert(s==u || t==u);
	  if(s!=u){
	    assert(getProperty(tokenMap,t)==true);
	    if (getProperty(tokenMap, s) == false) {
#ifdef DEBUGGING
	      std::cout << " edge " << u << "-" << s;
	      std::cout << std::endl;
	      std::cout << "Edge score:" << cargo(*cit) << std::endl;
#endif
	      stack.push_back(std::make_pair(s,u));
	    }
	  }
	  else{
	    if(t!=u){
	      assert(getProperty(tokenMap,s)==true);
	      if (getProperty(tokenMap, t) == false) {
#ifdef DEBUGGING
		std::cout << " edge " << u << "-" << t;
		std::cout << std::endl;
		std::cout << "Edge score:" << cargo(*cit) << std::endl;
#endif
		stack.push_back(std::make_pair(t,u));
	      }
	    }
	    else{
	      assert(false);
	    }
	  }
	}
      }
    }
  }
  //Capture all vertices with no edges degree==0
  TVertexIterator it(g);
  TVertexDescriptor u;
  for(;!atEnd(it);goNext(it)) {
    u = getValue(it);
    if (getProperty(tokenMap, u) == false) {
#ifdef DEBUGGING
      std::cout << "Component" << maxlabel
		<< std::endl;
#endif
      assert(degree(g,u)==0);
      ++maxlabel;
      assignProperty(components, u, maxlabel);
      genomeMap.push_back(TGenomeVertexMap());
      TSize gname = genomeNames[sequenceId(g,u)];
      assert(sequenceId(g,u)<length(genomeNames));
      genomeMap[getProperty(components,u)].insert(std::make_pair(gname,u));
    }
  }
  return maxlabel+1;
}

template<typename TSpec, typename TComponents>
typename Size<Graph<TSpec> >::Type
connected_components_ranked(Graph<TSpec> const& g,
				      TComponents& components){
  //SEQAN_CHECKPOINT
  typedef typename Size<Graph<TSpec> >::Type TSize;
  typedef typename Iterator<Graph<TSpec>, EdgeIterator>::Type TEdgeIterator;
  typedef typename Iterator<Graph<TSpec>, VertexIterator>::Type TVertexIterator;
  typedef typename VertexDescriptor<Graph<TSpec> >::Type TVertexDescriptor;
  typedef typename EdgeDescriptor<Graph<TSpec> >::Type TEdgeDescriptor;
  
  typedef typename Iterator<Graph<TSpec>, AdjacencyIterator>::Type TAdjacencyIterator; 
  typedef typename Iterator<Graph<TSpec>, OutEdgeIterator>::Type TOutEdgeIterator;
  
  clear(components);
  resizeVertexMap(g,components);
  
#ifdef DEBUGGING	
  std::cout << "Calculating connected components" << std::endl;
#endif
  
  // Initialization
  String<bool> tokenMap;
  fill(tokenMap, getIdUpperBound(_getVertexIdManager(g)), false);
  fill(components,getIdUpperBound(_getVertexIdManager(g)), 0);
  for(unsigned int i=0;i<getIdUpperBound(_getVertexIdManager(g));++i){
    assignProperty(components,i,0);
    assert(getProperty(components,i)==0);
  }
  TEdgeIterator itE(g);
  std::vector<TEdgeDescriptor> ccedges;
  int maxlabel=-1;
  for(;!atEnd(itE);goNext(itE)){
    ccedges.push_back(*itE);
  }
  sort(ccedges.begin(),ccedges.end(),edgescorecmp<std::map<TEdgeDescriptor,float> >());
  std::vector<std::pair<TVertexDescriptor,TVertexDescriptor> > stack;
  //outer loop ensures we visit disconnected subgraphs, considering most consistent edges first
  for(typename std::vector<TEdgeDescriptor>::reverse_iterator cit = ccedges.rbegin();cit!=ccedges.rend();cit++){
    TVertexDescriptor s = getSource(*cit);
    TVertexDescriptor t = getTarget(*cit);
    if (getProperty(tokenMap, s) == false){
      if(getProperty(tokenMap,t) == false) {
	stack.push_back(std::make_pair(s, t));
	stack.push_back(std::make_pair(t, s));
      }
      else{
	stack.push_back(std::make_pair(s, t));
      }
    }
    else{
      if(getProperty(tokenMap,t) == false) {
	stack.push_back(std::make_pair(t, s));
      }
    }
    while(!stack.empty()){
      std::pair<TVertexDescriptor,TVertexDescriptor> & node = stack.back();
      TVertexDescriptor u = node.first;
      TVertexDescriptor prev = node.second;
      stack.pop_back();
      if(getProperty(tokenMap,u)==false){
	assert(getProperty(components,u)==0);
#ifdef DEBUGGING
	std::cout << "New node " << u << " " << getProperty(components,u) 
		  << std::endl;
#endif
	assignProperty(tokenMap, u, true);
	int label=-1;
	if(getProperty(tokenMap,prev)==false){
	  //Use new label
	  ++maxlabel;
#ifdef DEBUGGING
	  std::cout << "Starting new component " << maxlabel << " for " << u << std::endl;
#endif
	  label=maxlabel;
	}
	else{
	  //there is already a label
	  int prevlabel=getProperty(components,prev);
#ifdef DEBUGGING
	  std::cout << "Adding to component " << getProperty(components,prev) << " for " << u << std::endl;
#endif
	  label=prevlabel;
	}

	assignProperty(components, u, label);
#ifdef DEBUGGING
	std::cout << "V:" << u << " component " << getProperty(components,u) << std::endl;
#endif
	//Add all edges from u to ccedges
	std::vector<TEdgeDescriptor> ccedges;
	//Edge iterator
	TOutEdgeIterator itOut(g, u);
	for(;!atEnd(itOut); ++itOut) {
	  if(getProperty(tokenMap,getSource(*itOut)) || getProperty(tokenMap,getTarget(*itOut))){
	    ccedges.push_back(*itOut);
	  }
	}
	//assert(ccedges.size()==degree(g,u));
	sort(ccedges.begin(),ccedges.end(),edgescorecmp<std::map<TEdgeDescriptor,float> >());
	for(typename std::vector<TEdgeDescriptor>::iterator cit = ccedges.begin();cit!=ccedges.end();cit++){
	  TVertexDescriptor s = getSource(*cit);
	  TVertexDescriptor t = getTarget(*cit);
	  assert(s==u || t==u);
	  if(s!=u){
	    assert(getProperty(tokenMap,t)==true);
	    if (getProperty(tokenMap, s) == false) {
#ifdef DEBUGGING
	      std::cout << " edge " << u << "-" << s;
	      std::cout << std::endl;
#endif
	      stack.push_back(std::make_pair(s,u));
	    }
	  }
	  else{
	    if(t!=u){
	      assert(getProperty(tokenMap,s)==true);
	      if (getProperty(tokenMap, t) == false) {
#ifdef DEBUGGING
		std::cout << " edge " << u << "-" << t;
		std::cout << std::endl;
#endif
		stack.push_back(std::make_pair(t,u));
	      }
	    }
	    else{
	      assert(false);
	    }
	  }
	}
      }
    }
  }
  //Capture all vertices with no edges degree==0
  TVertexIterator it(g);
  TVertexDescriptor u;
  for(;!atEnd(it);goNext(it)) {
    u = getValue(it);
    if (getProperty(tokenMap, u) == false) {
#ifdef DEBUGGING
      std::cout << "Component" << maxlabel
		<< std::endl;
#endif
      assert(degree(g,u)==0);
      ++maxlabel;
      assignProperty(components, u, maxlabel);
    }
  }
  return maxlabel+1;
}


//////////////////////////////////////////////////////////////////////////////

/**
.Function.alignmentEvaluation:
..summary:Given a multiple alignment, this function calculates all kinds of alignment statistics.
..cat:Graph
..signature:
alignmentEvaluation(graph, score_type, gapExCount, gapCount, pairCount, numPairs, len)
..param.graph:An alignment graph.
...type:Spec.Alignment Graph
..param.score_type:A score object.
...type:Class.Score
..param.gapExCount:Number of gap extensions.
..param.gapCount:Number of gaps.
..param.pairCount:Number of aligned pairs.
..param.numPairs:Counter for each pair.
..param.len:Alignment length.
..returns:Score of the alignment.
*/
template<typename TStringSet, typename TCargo, typename TSpec, typename TScore, typename TSize> 
//inline typename Value<TScore>::Type
s_score 
alignmentEvaluationCustom(Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
		    TScore const& score_type,
		    TSize& gapExCount,
		    TSize& gapCount,
		    TSize& pairCount,
		    TSize& pairIdent,
		    String<TSize>& numPairs,
		    String<TSize>& numIdentCols,
		    TSize& len,
		    TSize& totalLen)
{
  //SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename Value<TScore>::Type TScoreValue;
	typedef typename Value<typename Value<TStringSet>::Type>::Type TAlphabet;
	TSize alphSize = ValueSize<TAlphabet>::VALUE;
	s_score sscore;
	// Initialization;
	gapExCount = 0;
	gapCount = 0;
	pairCount = 0;
	clear(numPairs);

	// Convert the graph
	String<char> mat;
	convertAlignment(g, mat);
	char gapChar = gapValue<char>();

	TScoreValue gap = scoreGapExtend(score_type);
	TScoreValue gapOpen = scoreGapOpen(score_type);
	TSize nseq = length(stringSet(g));
	len = length(mat) / nseq;

	for(TSize i = 0; i<nseq; ++i) {
	  totalLen += length(stringSet(g)[i]);
	}

	fill(numIdentCols, nseq+1, 0);
	for(TSize j=0; j<=nseq; ++j) {
	  assert(numIdentCols[j]==0);
	}
	char c;
	for(TSize k=0;k<len; ++k) {
	  TSize numIdents=0;
	  for(TSize j=0; j<nseq; ++j) {
	    if (value(mat, j*len+k) != gapChar) {
	      if(numIdents==0){
		 c = TAlphabet(value(mat, j*len+k));
		 ++numIdents;
	      }
	      else{
		if(TAlphabet(value(mat, j*len+k))==c){
		  ++numIdents;
		}
		else{
		  numIdents=0;
		  break;
		}
	      }
	    }
	  }
	  assert(numIdents<=nseq);
	  numIdentCols[numIdents]++;
	}

	bool gapOpeni = false;
	bool gapOpenj = false;
	TScoreValue totalScore = 0;
	fill(numPairs, alphSize * alphSize, 0);
	for(TSize i = 0; i<nseq-1; ++i) {
	  for(TSize j=i+1; j<nseq; ++j) {
			for(TSize k=0;k<len; ++k) {
				if (value(mat, i*len+k) != gapChar) {
					if (value(mat, j*len + k) != gapChar) {
						gapOpeni = false;
						gapOpenj = false;
						++pairCount;
						if(TAlphabet(value(mat, i*len+k)) == TAlphabet(value(mat, j*len + k))){
						  ++pairIdent;
						}
						TSize index1 = ordValue(TAlphabet(value(mat, i*len+k)));
						TSize index2 = ordValue(TAlphabet(value(mat, j*len + k)));
						value(numPairs, index1 * alphSize + index2) += 1;
						totalScore += score(const_cast<TScore&>(score_type), TAlphabet(value(mat, i*len+k)), TAlphabet(value(mat, j*len + k)));
					} else {
						if (gapOpenj) {
							++gapExCount;
							totalScore += gap;
						} else {
							gapOpenj = true;
							++gapCount;
							totalScore += gapOpen;
						}
					}
				} else if (value(mat, j*len + k) != gapChar) {
						if (gapOpeni) {
							++gapExCount;
							totalScore += gap;
						} else {
							++gapCount;
							gapOpeni = true;
							totalScore += gapOpen;
						}
				}
			}
		}
	}
	sscore.alignScore = totalScore;
	sscore.numGap = gapCount;
	sscore.numGapEx = gapExCount;
	sscore.numPairs = pairCount;
	sscore.numIdents = pairIdent;
	sscore.alignLen = len;
	sscore.totalLen = totalLen;
	sscore.colCount = numIdentCols;
	sscore.seqCount = nseq;
	assert(length(numIdentCols)==nseq+1);
	sscore.pairCount = numPairs;
	//return totalScore;
	return sscore;
}

//readBlockFile()
//Read set of LCBs 
//File format is
//I seq1 start-end orient seq2 ....
//V v1 v2 v3 ;
//List of vertices in each LCB on a line ending with ;
//start-end are ignored. Boundaries are determined by extent of the member anchors


template<typename TVertexDescriptor,
	 typename TNames,
	 typename TVertexOrientMap,
	 typename TVertexSeqMap,
	 typename TGraph>
void doReadBlockFile(const std::string & filename,
		     std::map<unsigned int, std::set<TVertexDescriptor> > & block2fragMap,
		     std::vector<std::vector<TVertexDescriptor> > & lcbs,
		     TNames &sequenceNames,
		     TVertexOrientMap & vertexOrientMap,
		     TVertexSeqMap & vertexSeqMap,
		     TGraph & g,
		     bool checkbounds){
  std::ifstream file;
  file.open(filename.c_str(), std::ios_base::in | std::ios_base::binary);
  typedef std::ifstream TFile;
  typedef Value<TFile>::Type TValue;
  std::string line;
  std::vector<unsigned int> currblock;
  std::vector<unsigned int> currlcb;
  //Map sequence name -> char
  std::map<std::string,char> sequenceOrientMap;
  std::map<std::string,std::pair<unsigned int, unsigned int> > sequenceCoordsMap;
  std::map<String<char>,int> seqNamesIdxMap;
  for(int i=0;i<length(sequenceNames);++i){
    seqNamesIdxMap[sequenceNames[i]]=i;
  }

  unsigned int vertexcount = 0;
  while(file){
    getline(file,line); 
    std::istringstream in(line);
    std::string c;
    in >> c;
    if(c == "V"){
      while(in>>c){
	if(c == ";"){
	  for(std::vector<unsigned int>::iterator it=currblock.begin();it!=currblock.end();it++){
	    assert(block2fragMap.find(*it)!=block2fragMap.end());
	    //currlcb.insert(currlcb.end(),block2fragMap[*it].begin(),block2fragMap[*it].end());

	    for(typename std::set<TVertexDescriptor>::iterator vit=block2fragMap[*it].begin();vit!=block2fragMap[*it].end();++vit){
	      assert(vertexSeqMap.find(*vit)!=vertexSeqMap.end());
	      std::string sname(toCString(sequenceNames[vertexSeqMap[*vit]]));
	      //Check if segment is reported as part of the LCB
	      //The LCB identification step may not report all sequences that are part of the anchor
	      if(sequenceOrientMap.find(sname) != sequenceOrientMap.end()){
		if(!checkbounds || fragmentBegin(g,*vit)>=sequenceCoordsMap[sname].first && fragmentBegin(g,*vit)<=sequenceCoordsMap[sname].second){
		  vertexOrientMap[*vit] = sequenceOrientMap[sname];
		  //add node to block
		  currlcb.push_back(*vit);
#ifdef DEBUGGING
		  std::cout << "Adding segment V:" << *vit << " from anchor:" << *it << std::endl;
#endif
		}
		else{
#ifdef DEBUGGING
		  std::cout << "Skipping out-of-bounds anchor segment " << *vit << " len:" << fragmentLength(g,*vit) << " from anchor " << *it << " on sequence " << sname << " " << sequenceOrientMap[sname]
			    << " fragmentBegin:" << fragmentBegin(g,*vit) << " bounds:" << sequenceCoordsMap[sname].first << "-" << sequenceCoordsMap[sname].second << " " << vertexOrientMap[*vit] << std::endl;
#endif
		}
	      }
	      else{
#ifdef DEBUGGING
		std::cout << "Skipping anchor segment " << *vit << " len:" << fragmentLength(g,*vit) << " from anchor " << *it << " on sequence " << sname << std::endl;
	        for(std::map<std::string,char>::iterator sit=sequenceOrientMap.begin();sit!=sequenceOrientMap.end();++sit){
	            std::cerr << sit->first << " " << sit->second << std::endl;
                }	
#endif 
		//currlcb.push_back(*vit);
		//assert(false);
	      }
	    }
	  }
	  vertexcount = vertexcount + currlcb.size();
	  lcbs.push_back(currlcb);
	  currblock.clear();
	  currlcb.clear();
	}
	else{
	  currblock.push_back(atoi(c.c_str()));
	  //Update vertex orientation
	}
      }
      currblock.clear();
    }
    else{
      sequenceOrientMap.clear();
      sequenceCoordsMap.clear();
      if(c == "I"){
	while(in>>c){
	  if(c != ";"){
	    //read sequence name
	    std::string seqname=c;
	    char orient;
	    std::string coords;
	    assert(seqname != ";");
	    //read orient
	    in >> orient;
	    assert(orient != ';');
	    assert(orient == '+' || orient == '-');
	    sequenceOrientMap[seqname] = orient;
	    //read coords
	    in >> coords;
	    std::string start;
	    std::string end;
	    std::istringstream coordsin(coords);
	    getline(coordsin, start, '-');
	    getline(coordsin, end, '-');
	    unsigned int startcoord = atoi(start.c_str());
	    unsigned int endcoord = atoi(end.c_str());
	    /*
	    if(orient == '-'){
	      int slen = length(getValueById(stringSet(g), seqNamesIdxMap[seqname]));
	      unsigned int tmpstartcoord=startcoord;
	      startcoord = slen-endcoord;
	      endcoord = slen-startcoord;
	    }
	    */
	    sequenceCoordsMap[seqname] = make_pair(startcoord,endcoord);
	  }
	}
      }
    }
  }
  std::cerr << "Read " << lcbs.size() << " LCBs containing " << vertexcount << " segments " << std::endl;
  //Sort based on length
  

}
template<typename TVertexDescriptor,
	 typename TNames,
	 typename TVertexOrientMap,
	 typename TGraph>
void readBlockFile(const std::string & filename,
		   std::map<unsigned int, std::set<TVertexDescriptor> > & block2fragMap,
		   std::vector<std::vector<unsigned int> > & lcbs,
		   TNames &sequenceNames,
		   TVertexOrientMap & vertexOrientMap,
		   TGraph & g,
		   bool checkbounds=false){
  typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
  std::map<TVertexDescriptor,unsigned> vertexSeqMap;
  TVertexIterator it(g);
  for(;!atEnd(it);goNext(it)) {
    vertexSeqMap[*it] = sequenceId(g,*it);
  }
  doReadBlockFile(filename,block2fragMap,lcbs,sequenceNames,vertexOrientMap,vertexSeqMap,g,checkbounds); 
  std::cerr << "Sorting LCBs by length" << std::endl; 
  //Sort LCBs in decreasing order by length
  std::map<int,int> lcbidxlenmap;
  std::vector<int> lcbsidx;
  for(unsigned int i=0;i<lcbs.size();i++){
    int totalsize=0;
    std::vector<unsigned int>::const_iterator vit;
    for(vit = lcbs[i].begin();vit!=lcbs[i].end();vit++){
      totalsize = totalsize + fragmentLength(g,*vit);
    }
    assert(totalsize!=0);
    lcbidxlenmap[i] = totalsize;
    lcbsidx.push_back(i);
  }

  assert(lcbidxlenmap.size()==lcbs.size());
  assert(lcbsidx.size()==lcbs.size());
  sort(lcbsidx.begin(),lcbsidx.end(),lcblencmp<std::map<int,int> >(lcbidxlenmap));

  std::vector<std::vector<unsigned int> > newlcbs;
  for(std::vector<int>::iterator lit = lcbsidx.begin();lit != lcbsidx.end();++lit){
    assert((unsigned int)*lit<lcbs.size());
    newlcbs.push_back(lcbs[*lit]);
    if(lit != lcbsidx.begin()){
#ifdef DEBUGGING
      std::cout << "LCB: " << *lit << " " << lcbidxlenmap[*lit] << " <= " << lcbidxlenmap[*(lit-1)] << std::endl;
#endif
      assert(lcbidxlenmap[*lit]<=lcbidxlenmap[*(lit-1)]);
    }
  }
  lcbs=newlcbs;
}
/*
template<typename TSeqs>
void do_segmentation_MERCATOR(){
  std::fstream strm;
  //(2)Support for Mercator
  //G.chroms
  std::string genomestr;
  for(unsigned int i=0;i<nSeq;i++){
    std::ostringstream currfilename;
    currfilename << "G" << i;
    genomestr += currfilename.str();
    genomestr += " ";
    currfilename << ".chroms";
    strm.open(currfilename.str().c_str(), std::ios_base::out | std::ios_base::trunc);
    strm << "S" << i << "\t" //seqname
	 << length(seqSet[i]) << std::endl; //chromLength
    strm.close();
  }
  for(unsigned int i=0;i<nSeq;i++){
    std::ostringstream currfilename;
    currfilename << "G" << i << ".anchors";
    strm.open(currfilename.str().c_str(), std::ios_base::out | std::ios_base::trunc);
    typename std::vector<TBlock>::const_iterator bit = blocks.begin();
    for(bit = blocks.begin();
	bit!=blocks.end();
	bit++){
      if(bit->currentSeq==i){
	//mercator format
	strm << bit->currentSeq << "." << bit->c << "\t" //anchorname
	     << "S" << bit->currentSeq << "\t" //seqname
	     << bit->orient << "\t" //strand
	     <<bit->begCoord << " " << bit->endCoord << "\t" //startCoord endCoord 0-based half open interval [start, end)
	     << 1 << std::endl; //isCoding
      }
    }
    strm.close();
  }  
  //G1-G2.hits
  std::map<std::pair<unsigned int,unsigned int>,std::set<int> > hitMap;
  std::map<std::pair<unsigned int,unsigned int>,std::set<int> >::iterator hmit;
  bool inserted;
  
  bit2 = blocksbycomponent.begin();
  for(; bit2!= blocksbycomponent.end();bit2++){//all cc
    std::vector<TBlock> currblocks = bit2->second;
    typename std::vector<TBlock>::iterator it1,it2;
    for(it1=currblocks.begin();it1!=currblocks.end();it1++){
      for(it2=currblocks.begin();it2!=currblocks.end();it2++){
	if(it1!=it2){
	  std::pair<unsigned int,unsigned int> key;
	  if(it1->currentSeq<it2->currentSeq){
	    key = std::make_pair(it1->currentSeq,it2->currentSeq);
	  }
	  else{
	    key = std::make_pair(it2->currentSeq,it1->currentSeq);
	  }
	  
	  std::pair<std::map<std::pair<unsigned int,unsigned int>,std::set<int> >::iterator,bool> s 
	    = hitMap.insert(std::make_pair(key,std::set<int>()));
	  hmit = s.first;
	  inserted = s.second;
	  hmit->second.insert(it1->c);
	}
      }
    }
  }
  for(hmit=hitMap.begin();hmit!=hitMap.end();hmit++){
    std::ostringstream currfilename;
    currfilename << "G" << hmit->first.first << "-" << "G" << hmit->first.second << ".hits";
    strm.open(currfilename.str().c_str(), std::ios_base::out | std::ios_base::trunc);
    for(std::set<int>::iterator it=hmit->second.begin();it!=hmit->second.end();it++){
      strm << hmit->first.first  << "." << *it << "\t"  //anchorName1
	   << hmit->first.second << "." << *it << "\t" //anchorName2
	   << 1 << "\t"
	   << 1 << std::endl;
    }
    strm.close();
  }
}
*/

template<typename TBlock,
	 typename TVertexDescriptor, 
	 typename TMSAOptions, 
	 typename TSeqs, 
	 typename TNames,
	 typename TGenomeNames,
	 typename TVertexOrientMap,
	 typename TGraph>
void do_segmentation_ENREDO(std::vector<TBlock> & blocks,
			    std::vector<std::vector<unsigned int> > & lcbs,
			    std::map<unsigned int, std::set<TVertexDescriptor> > & block2fragMap, 
			    std::string distance,
			    std::string minlen,
			    TMSAOptions const &msaOpt,
			    TSeqs & seqSet,
			    TNames & sequenceNames,
			    TGenomeNames & genomeNames,
			    TVertexOrientMap & vertexOrientMap,
			    TGraph & g){
  /*
    From Enredo README
  The input file contains the result of mapping a set of anchors onto several
  genomes. Anchors are expected to be sorted by organism, chromosome and
  position. Each line should correspond to an anchor and each line contains 6
  values separated by tabs. The six values are: the anchor name (a string),
  the species name (a string), the chromosome name (a string), the start
  position (an integer value), the end position (an integer value), the strand
  (either + or -) and the score (a real value). Here is an example:
  
  A1      Spcs1   X       53      85      +       123
  B1      Spcs1   X       458     498     +       11
  C1      Spcs1   X       3601    3639    +       434
  B1      Spcs1   X       5480    5520    +       1
  D1      Spcs1   X       6479    6510    +       41
  A       Spcs1   Y       1379    4410    +       1567
  E       Spcs1   Y       5879    5910    +       311
  E       Spcs1   Y       6479    6510    +       217
  D       Spcs1   Y       6567    6593    +       135
  */
  std::fstream strm;
  std::fstream strm2;
  String<char> pf = msaOpt.outfile;
  char * pfilename = toCString(pf);
  std::string projfilename(pfilename);
  projfilename = projfilename + "enredo.anchors";
  std::string idxfilename(pfilename);
  idxfilename = idxfilename + "enredo.idx";
  std::cerr << "Writing ENREDO anchors to " << projfilename.c_str() << std::endl;
  strm.open(idxfilename.c_str(), std::ios_base::out | std::ios_base::trunc);
  strm2.open(projfilename.c_str(), std::ios_base::out | std::ios_base::trunc);
  //enredo anchors file
  //G.anchors
  unsigned int nseq = length(seqSet);
  for(unsigned int i=0;i<nseq;i++){
    std::vector<TBlock *> seqblocks;
    int idx=0;
    strm << i << " " << sequenceNames[i] << std::endl;
    for(typename std::vector<TBlock>::const_iterator bit = blocks.begin(); bit!=blocks.end();bit++){
      if(bit->currentSeq==i){
	seqblocks.push_back(&(blocks[idx]));
      }
      idx++;
    }
    
    std::sort(seqblocks.begin(),seqblocks.end(),blockorder<TBlock,unsigned int>(i));	  
    unsigned int blen = seqblocks.size();
    //if(blen>1){
      for(unsigned int j=0;j<blen;j++){
	strm2 << seqblocks[j]->c << "\t" 
	      << i << "\t"
	      << genomeNames[i] << "\t" 
	  //1-start base coordinates
	      << seqblocks[j]->begCoord+1 << "\t" << seqblocks[j]->endCoord << "\t" 
	      << seqblocks[j]->orient << "\t" 
	      << seqblocks[j]->endCoord - seqblocks[j]->begCoord 
	      << std::endl;
	//}
    }
  }	
  /*  
  unsigned int nSeq = length(seqSet);
  for(unsigned int i=0;i<nSeq;i++){
    strm << i << " " << sequenceNames[i] << std::endl;
    typename std::vector<TBlock>::const_iterator bit = blocks.begin();
    for(bit = blocks.begin();
	bit!=blocks.end();
	bit++){
      if(bit->currentSeq==i){
	//enredo format
	strm2 << bit->c << "\t" 
	      << bit->currentSeq << "\t" 
	      << genomeNames[i] << "\t";
	strm2 << bit->begCoord << "\t" 
	      << bit->endCoord << "\t";
	strm2 << bit->orient << "\t"
	      << bit->endCoord-bit->begCoord<< std::endl;
      }
    }
  }
  */
  strm.close();
  strm2.close();
  //(2)Sort anchors
  std::string sortedprojfilename(projfilename+".sorted");
  std::string sortcmd = "sort -k 2,3 -k 4n,4n < " + projfilename + " > " + sortedprojfilename;
  int res = system(sortcmd.c_str());
  if(res!=0){
    perror("Could not run system command: ");
    std::cerr << sortcmd.c_str() << std::endl 
	      << "SYSTEM:" << res << std::endl;
    exit(1);
  }
  //(3)Run Enredo
  std::string mugsyinstall = std::string(std::getenv("MUGSY_INSTALL"));
  assert(mugsyinstall.length()>0);
  std::string cmd = mugsyinstall+"/enredo ";
  std::string stdoutfilename(boost::lexical_cast<std::string>(getpid())+"lcbs.out");
  std::string stderrfilename(boost::lexical_cast<std::string>(getpid())+"synchain-mugsy.out");
  char * enredoenvopts = std::getenv("ENREDO_OPTS");
  std::string enredoopts;
  if(enredoenvopts==NULL || strlen(enredoenvopts)==0){
       enredoopts = std::string(" --min-score 0 --max-ratio 0 ") + std::string(" --min-length ") + minlen + std::string(" --max-gap-length ") + distance + std::string(" --min-anchors 1 ");
  }
  else{
       enredoopts = std::string(enredoenvopts);
  }
  cmd = cmd + enredoopts + " " + sortedprojfilename
    + " | "+mugsyinstall+"/enredo2mugsy.pl "+idxfilename+" > "+stdoutfilename+" 2> "+stderrfilename;
  assert(cmd.length()>0);
  //#ifdef DEBUGGING
  std::cerr << "Running " << cmd.c_str() << std::endl;
  //#endif
  res = system(cmd.c_str());
  if(res!=0){
    perror("Could not run system command: ");
    std::cerr << cmd.c_str() << std::endl 
	      << "SYSTEM:" << res << std::endl;
    exit(1);
  }
  assert(res==0);
  //(3) Read output file to obtain list of LCBs
  readBlockFile(stdoutfilename,
		block2fragMap,
		lcbs,
		sequenceNames,
		vertexOrientMap,
		g,
		true); //must check bounds
}

template<typename TBlock,
	 typename TNames,
	 typename TGenomeNames>
void writeProjectionFile(std::string projfilename,
			 std::vector<TBlock> & blocks,
			 TNames & sequenceNames,
			 TGenomeNames & genomeNames){
  std::fstream strm;
  strm.open(projfilename.c_str(), std::ios_base::out | std::ios_base::trunc);
  unsigned int nseq = length(sequenceNames);
  assert(nseq==length(sequenceNames));
  assert(nseq==length(genomeNames));
  unsigned int cdist=0;
  //(1)Project blocks onto each sequence and write to a file
  //Topological sort of the blocks over each sequence
  //projecting block onto sequence and printing neighbors (n->n+1)
  for(unsigned int i=0;i<nseq;i++){
    std::vector<TBlock *> seqblocks;
    int idx=0;
    for(typename std::vector<TBlock>::const_iterator bit = blocks.begin(); bit!=blocks.end();bit++){
      if(bit->currentSeq==i){
	seqblocks.push_back(&(blocks[idx]));
      }
      idx++;
    }
    
    std::sort(seqblocks.begin(),seqblocks.end(),blockorder<TBlock,unsigned int>(i));	  
    unsigned int blen = seqblocks.size();
    if(blen>1){
      for(unsigned int j=0;j<blen;j++){
	if(j<blen-1){
	  assert(seqblocks[j]->currentSeq==seqblocks[j+1]->currentSeq);
	  cdist = std::abs((long int)(seqblocks[j]->endCoord - seqblocks[j+1]->begCoord));
	  strm << seqblocks[j]->c << " " << seqblocks[j+1]->c << " " << sequenceNames[i] << " " << cdist << " " << genomeNames[i] << " " 
	       << seqblocks[j]->orient << " " << seqblocks[j+1]->orient << " " 
	       << seqblocks[j]->begCoord << " " << seqblocks[j]->endCoord << " " 
	       << seqblocks[j+1]->begCoord << " " << seqblocks[j+1]->endCoord
	       << std::endl;
	}
      }
    }
    else{
      if(blen==1){
	strm << seqblocks[0]->c << " " << seqblocks[0]->c << " " << sequenceNames[i] << " " << 0 << " " << genomeNames[i] << " " 
	     << seqblocks[0]->orient << " " << seqblocks[0]->orient << " " 
	     << seqblocks[0]->begCoord << " " << seqblocks[0]->endCoord << " " 
	     << seqblocks[0]->begCoord << " " << seqblocks[0]->endCoord
	     << std::endl;
      }
    }
  }
  strm.close();
    

}

template<typename TBlock,
	 typename TVertexDescriptor, 
	 typename TMSAOptions, 
	 typename TNames,
	 typename TGenomeNames,
	 typename TVertexOrientMap,
	 typename TGraph>
void do_segmentation_MUGSY( std::vector<TBlock> & blocks,
			    std::vector<std::vector<unsigned int> > & lcbs,
			    std::map<unsigned int, std::set<TVertexDescriptor> > & block2fragMap, 
			    std::string distance,
			    std::string minlen,
			    TMSAOptions const &msaOpt,
			    TNames & sequenceNames,
			    TGenomeNames & genomeNames,
			    TVertexOrientMap & vertexOrientMap,
			    TGraph & g){

  std::fstream strm;
  String<char> pf = msaOpt.outfile;
  char * pfilename = toCString(pf);

  //(1)Write projection file

  std::string projfilename(pfilename);
  projfilename = projfilename + "projections.out";
  writeProjectionFile(projfilename,blocks,sequenceNames,genomeNames);
  //(2)Run synchain-mugsy using the projection
  std::string mugsyinstall = std::string(std::getenv("MUGSY_INSTALL"));
  assert(mugsyinstall.length()>0);
  std::string cmd = "cat "+projfilename+" | "+mugsyinstall+"/synchain-mugsy ";
  std::string stdoutfilename(boost::lexical_cast<std::string>(getpid())+"lcbs.out");
  std::string stderrfilename(boost::lexical_cast<std::string>(getpid())+"synchain-mugsy.out");
  cmd = cmd + distance + " " + minlen + " "+ minlen
    + " > "+stdoutfilename+" 2> "+stderrfilename;
  assert(cmd.length()>0);
  //#ifdef DEBUGGING
  std::cerr << "Running " << cmd.c_str() << std::endl;
  //#endif
  #ifdef TIMING 
  time(&now);
  std::cerr << "TIME PRE-SYNCHAIN:" << lasttime << " " << now << " " << now-lasttime << std::endl;
  lasttime=now;
  #endif 
  int res = system(cmd.c_str());
  #ifdef TIMING 
  time(&now);
  std::cerr << "TIME SYNCHAIN:" << lasttime << " " << now << " " << now-lasttime << std::endl;
  lasttime=now;
  #endif
  if(res!=0){
    perror("Could not run system command: ");
    std::cerr << cmd.c_str() << std::endl 
	      << "SYSTEM:" << res << std::endl;
    exit(1);
  }
  assert(res==0);
  //(3) Read output file to obtain list of LCBs
  readBlockFile(stdoutfilename,
		block2fragMap,
		lcbs,
		sequenceNames,
		vertexOrientMap,
		g,
		false); //No bounds check
  if(res==0){
#ifdef KEEPCHAINTMP
    ;
#else
    unlink(stdoutfilename.c_str());
    unlink(stderrfilename.c_str());
    unlink(projfilename.c_str());
#endif

  }
  #ifdef TIMING 
  time(&now);
  std::cerr << "TIME POST-SYNCHAIN:" << lasttime << " " << now << " " << now-lasttime << std::endl;
  lasttime=now;
  #endif
}


//Build ungapped profiles from the segment graph
//Output is array of TBlocks needed for mugsy-chaining
template<typename TGraph,
	 typename TComponentMap,
	 typename TComponent,
	 typename TBlock,
	 typename TName,
	 typename TLoc,
	 typename TNames
	 >
void convertCC2Blocks(TGraph &g, 
		      TComponentMap& component,
		      std::map<std::pair<TComponent,TComponent>,TBlock *> & componentVertexMap,
		      std::vector<std::vector<TBlock> > & blocksbycomponent,
		      std::map<TName,std::vector<TLoc> >&aintervals,
		      TNames & sequenceNames){
  
  typedef typename Id<TGraph>::Type TIdType;
  typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
  //typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
  typedef typename Size<TGraph>::Type TSize;
  typedef std::pair<TIdType, TSize> TKey;
  typedef std::map<TKey, TVertexDescriptor> TPosToVertexMap;
  typedef FragmentInfo<TIdType, TSize> TFragmentInfo;

  // data_pvMap is an STL Map to retrieve a vertex given SeqId, Position
  // first.first == seqId
  // first.second == pos
  // second == VertexDescriptor
  typename TPosToVertexMap::const_iterator it1 = g.data_pvMap.begin();
  typename TPosToVertexMap::const_iterator it1End = g.data_pvMap.end();
  typedef typename Position<TGraph>::Type TPos;

  std::map<TComponent,int> seqsPerComponent;
  
  TPos begCoord,endCoord;
  char orient='?';

  int lostbp=0;
  int numlostv=0;
  
  //Track number of sequences per component
  for(;it1!=it1End;++it1) {
    TVertexDescriptor currV = it1->second;
    if(currV != getNil<TVertexDescriptor>()){
      assert(getProperty(component,currV)==component[currV]);
      TComponent c = getProperty(component, currV);
      assert(c < blocksbycomponent.size());
      if(seqsPerComponent.find(c)!=seqsPerComponent.end()){
	seqsPerComponent[c]++;
      }
      else{
	seqsPerComponent.insert(std::make_pair(c,1));
      }
    }
  }
  for(it1 = g.data_pvMap.begin();it1!=it1End;++it1) {
    TVertexDescriptor currV = it1->second;
    if(currV != getNil<TVertexDescriptor>()){
      assert(getProperty(component,currV)==component[currV]);
      TComponent c = getProperty(component, currV);
      assert(c < blocksbycomponent.size());
      TSize currentSeq = sequenceId(g,currV);
      if(seqsPerComponent[c] > 1){
#ifdef DEBUGGING
	std::cout << "Component " << c << " V:" << currV << " seq:" << currentSeq << " degree:" << degree(g,currV) << " coord:" << fragmentBegin(g,currV) << std::endl;
#endif    
	//First block for currentseq
	typename std::map<std::pair<TComponent,TComponent>,TBlock *>::iterator fit = componentVertexMap.find(std::make_pair(c,currentSeq));
	if(fit==componentVertexMap.end()){
	  begCoord = fragmentBegin(g,currV);
	  assert((int)begCoord>=0);
	  endCoord = begCoord+fragmentLength(g,currV);
	  orient='?';
	  typename std::vector<TBlock>::iterator bit = blocksbycomponent[c].insert(blocksbycomponent[c].end(),
										   TBlock(c,currentSeq,begCoord,endCoord,orient,currV));
	  componentVertexMap[std::make_pair(c,currentSeq)] = &(*bit);
	  //blocksbycomponent[c].push_back(TBlock(c,currentSeq,begCoord,endCoord,orient,currV));
	  //unsigned int idx = blocksbycomponent[c].size()-1;
	  //componentVertexMap[std::make_pair(c,currentSeq)] = &(blocksbycomponent[c][idx]);
#ifdef DEBUGGING
	  std::cout << "Adding component " << c << " seq:" << currentSeq << " coords" 
		    << begCoord << "-" << endCoord << " o:" << orient << " V:" << currV << std::endl;
#endif
	}
	else{
	  //Block already inserted
	  TBlock * blk = fit->second;
	  blk->addVertex(g,currV);
#ifdef DEBUGGING
	  std::cout << "Adding vertex to component " << c << " seq:" << currentSeq 
		    << " coords:" << begCoord << "-" << endCoord << " V:" << currV << std::endl;
#endif
	}
      }
      else{
	//Repetitive sequence
	/* Remove to improve reporting of unique sequences
	if(degree(g,currV)>0){
	  lostbp = lostbp + fragmentLength(g,currV);
	  numlostv++;
	  typename std::map<TName,std::vector<TLoc > >::iterator ait = aintervals.find(sequenceNames[currentSeq]);
	  if(ait==aintervals.end()){      
	    aintervals.insert(std::make_pair(sequenceNames[currentSeq],std::vector<TLoc >()));
	  }
	  ait = aintervals.find(sequenceNames[currentSeq]);
	  assert(ait!=aintervals.end());
	  TLoc t1,t2;
	  t1.first = fragmentBegin(g,currV);
	  t1.second = 1;
	  t1.blocknum = 0;
	  ait->second.push_back(t1);
	  t2.first = t1.first+fragmentLength(g,currV);
	  t2.second = -1;
	  t2.blocknum = 0;
	  ait->second.push_back(t2);
	}
	*/
      }
    }
  }
  std::cerr << "Disconnected " << numlostv << " vertices marking " 
	    << lostbp << " aligned bp" << std::endl;
}


//Assign orientation using greedy approach
//Start assignment of edges with best consistency score.
//Break ties with positional score
//?TODO?: Break inconsistent edges
//Output:Blocks:          std::vector<TBlock>
//       vertexOrientMap: map vertex->orient
template<typename TGraph, 
	 typename TBlock, 
	 typename TVertexOrientMap,
	 typename TEdgeDescriptor>
void assignBlockOrientation(TGraph &g,
			    std::vector<std::vector<TBlock> > &blocksbycomponent, 
			    std::vector<TBlock> &blocks,
			    TVertexOrientMap &vertexOrientMap,
			    std::map<TEdgeDescriptor,float> &posScores){
  typedef unsigned int TSize;
  typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
  typename std::vector<std::vector<TBlock> >::iterator bit2 = blocksbycomponent.begin();

  //Track number of conflicting orientation assignments
  int conflicts=0;
  std::set<TEdgeDescriptor> conflictEdges;
  bool removeConflictingEdges=false;

  for(; bit2!= blocksbycomponent.end();bit2++){//all cc
    if(bit2->size()>0){
      std::vector<TEdgeDescriptor> ccedges;
      //std::vector<TBlock> * currblocks = bit2;
      int unorientedSegments=bit2->size();
#ifdef DEBUGGING
      std::cout << "Examining block with " << unorientedSegments << " sequences" << std::endl;
#endif
      std::map<TSize,TBlock *> seqBlockMap;
      //capture all edges in component
      //For all segments on seqs i+1->numseqs{
      for(unsigned int i=0;i<bit2->size();i++){
	for(unsigned int j=i+1;j<bit2->size();j++){
	  //all members of block
	  for(typename std::vector<TVertexDescriptor>::iterator vit=bit2->at(i).currV.begin();vit!=bit2->at(i).currV.end();vit++){
	    TVertexDescriptor currV0 = *vit;
	    //assert(degree(g,currV0)>0);
	    for(typename std::vector<TVertexDescriptor>::iterator vit1=bit2->at(j).currV.begin();vit1!=bit2->at(j).currV.end();vit1++){
	      TVertexDescriptor currV = *vit1;
	      //assert(degree(g,currV)>0);
	      TEdgeDescriptor ed = findEdge(g,currV0,currV);
	      if(ed!=0){
		ccedges.push_back(ed);
	      }
	      else{
		assert(findEdge(g,currV,currV0)==0);
	      }
	    }
	  }
	}
	seqBlockMap.insert(std::make_pair(bit2->at(i).currentSeq,&bit2->at(i)));
#ifdef DEBUGGING
	std::cout << "Block " << i << " seq:" << bit2->at(i).currentSeq  << " vertices:" << bit2->at(i).currV.size() << std::endl;
#endif
#ifdef DEBUGGING
	for(typename std::vector<TVertexDescriptor>::iterator vit=bit2->at(i).currV.begin();vit!=bit2->at(i).currV.end();vit++){
	  TVertexDescriptor currV0 = *vit;
	  std::cout << "V:"<<currV0 << std::endl;
	}
#endif
      }
      assert(seqBlockMap.size()==bit2->size());
#ifdef DEBUGGING
      std::cout << "Number of edges " << ccedges.size() << std::endl;
#endif
      //Sort edges in order of 
      //(1) consistency
      //(2) posscore
      //This way the most consistent and syntenic edges should determine the 
      //relative orientation of segments in the block
      sort(ccedges.begin(),ccedges.end(),edgescorecmp<std::map<TEdgeDescriptor,float> >(&posScores));
      //traverse edges in decreasing order ranked by consistency, posScores
      typename std::vector<TEdgeDescriptor>::reverse_iterator eit=ccedges.rbegin();
      TEdgeDescriptor ed = *eit;
      TVertexDescriptor v1 = getSource(ed);
      TVertexDescriptor v2 = getTarget(ed);
      assert(seqBlockMap.find(sequenceId(g,v1))!=seqBlockMap.end());
      assert(seqBlockMap.find(sequenceId(g,v2))!=seqBlockMap.end());
      TBlock * blockv1 = seqBlockMap[sequenceId(g,v1)];
      TBlock * blockv2 = seqBlockMap[sequenceId(g,v2)];
      assert(blockv2->orient == '?');
      assert(blockv1->orient == '?');
      blockv1->orient = '+';
      assert(cargo(ed)!=0);
      if(cargo(ed)>0){
#ifdef DEBUGGING
	std::cout << " SAME ORIENT " << blockv1->orient << std::endl;
#endif
	blockv2->orient = blockv1->orient;
      }
      else{
#ifdef DEBUGGING
	std::cout << " OPPOSITE ORIENT OF " << blockv1->orient << std::endl;
#endif
	blockv2->orient = (blockv1->orient=='+'?'-':'+');
      }
#ifdef DEBUGGING
      std::cout << "Examining edge " << " " << v1 << "-" << v2 
		<< blockv1->orient << " " << blockv2->orient
		<< std::endl;
#endif
      //check all currV in blockv1,blockv2 for consistency with this assignment
      //First two blocks are oriented relative to each other
      unorientedSegments-=2;
      eit++;
      //Propogate relative orientation through the component graph
      while(unorientedSegments>0){
#ifdef DEBUGGING
	std::cout << "Num blocks unoriented " << unorientedSegments << std::endl;
#endif
	for(;eit!=ccedges.rend();eit++){
	  ed = *eit;
	  v1 = getSource(ed);
	  v2 = getTarget(ed);
	  assert(v1!=v2);
	  assert(seqBlockMap.find(sequenceId(g,v1))!=seqBlockMap.end());
	  assert(seqBlockMap.find(sequenceId(g,v2))!=seqBlockMap.end());
	  blockv1 = seqBlockMap[sequenceId(g,v1)];
	  blockv2 = seqBlockMap[sequenceId(g,v2)];
#ifdef DEBUGGING
	  std::cout << "Examining edge:" << " " << v1 << "-" << v2 
		    << blockv1->orient << " " << blockv2->orient
		    << std::endl;
#endif
	  if(blockv1->orient == '?'){
	    if(blockv2->orient != '?'){
	      //assignment
	      unorientedSegments--;
	      if(cargo(ed)>0){
#ifdef DEBUGGING
		std::cout << " SAME ORIENT " << blockv2->orient << std::endl;
#endif
		blockv1->orient = blockv2->orient;
	      }
	      else{
#ifdef DEBUGGING
		std::cout << " OPPOSITE ORIENT OF " << blockv2->orient << std::endl;
#endif
		blockv1->orient = (blockv2->orient=='+'?'-':'+');
	      }
	    }
	    else{
	      //no assignment
	    }
	  }
	  else{
	    if(blockv2->orient == '?'){
	      if(blockv1->orient != '?'){
		//assignment
		unorientedSegments--;
		if(cargo(ed)>0){
#ifdef DEBUGGING
		  std::cout << " SAME ORIENT " << blockv1->orient << std::endl;
#endif
		  blockv2->orient = blockv1->orient;
		}
		else{
#ifdef DEBUGGING
		  std::cout << " OPPOSITE ORIENT OF " << blockv1->orient << std::endl;
#endif
		  blockv2->orient = (blockv1->orient=='+'?'-':'+');
		}
	      }
	      else{
		//no assignment
	      }
	    }
	    else{
	      //already assigned
	      //check
	      assert(cargo(ed)!=0);
	      if(cargo(ed)>0){
		if(blockv1->orient!=blockv2->orient){
#ifdef DEBUGGING
		  std::cout << "Conflicting orientation. Edge:" << cargo(ed) 
			    << " for vertices V1:" << v1 << "," << blockv1->orient 
			    << " V2:" << v2 << "," << blockv2->orient << std::endl;
#endif
		  //TODO break edge?
		  conflicts++;
		  conflictEdges.insert(ed);
		}
	      }
	      else{
		if(blockv1->orient==blockv2->orient){
#ifdef DEBUGGING
		  std::cout << "Conflicting orientation. Edge:" << cargo(ed) 
			    << " for vertices V1:" << v1 << "," << blockv1->orient 
			    << " V2:" << v2 << "," << blockv2->orient << std::endl;
#endif
		  //TODO break edge?
		  conflicts++;
		  conflictEdges.insert(ed);
		}
	      }
	    }
	  }
	}
	//Start search again at beginning of list of edges
	eit=ccedges.rbegin();
      }
    } //if 
    //copy final output
    blocks.insert(blocks.end(),bit2->begin(),bit2->end());
    for(unsigned int i=0;i<bit2->size();i++){
      //all members of block
      for(typename std::vector<TVertexDescriptor>::iterator vit=bit2->at(i).currV.begin();vit!=bit2->at(i).currV.end();vit++){
	vertexOrientMap[*vit] = bit2->at(i).orient;
      }
    }
  } //for all CC

  //Resolve conflicts if necessary
  if(removeConflictingEdges){
    for(typename std::set<TEdgeDescriptor>::iterator eit=conflictEdges.begin();eit!=conflictEdges.end();++eit){
      removeEdge(g,*eit);
    }
  }
  std::cerr << "Num conflicts: " << conflicts << " when assigning orientation" << std::endl;
}


//Fragments can be oriented
template<typename TFragment,
	 typename TGraph,
	 typename TVertexDescriptor,
	 typename TSize>
void buildFrag(TFragment & frag,
	       TGraph & g,
	       TVertexDescriptor vd1,
	       TVertexDescriptor vd2,
	       TSize id1,
	       int vd1len,
	       unsigned offset1,
	       char orient1,
	       TSize id2,
	       int vd2len,
	       unsigned offset2,
	       char orient2){
  if(orient1 == '-'){
    if(orient2 == '+'){
      //id1:- id2:+
      frag = TFragment(id1,
		       vd1len-(fragmentBegin(g,vd1)+fragmentLength(g,vd1))-offset1,
		       id2,
		       fragmentBegin(g,vd2)-offset2, 
		       fragmentLength(g,vd1),
		       true);
    }
    else{
      //id1:- id2:-
      assert(orient1 == '-');
      assert(orient2 == '-');
      frag = TFragment(id1,
		       vd1len-(fragmentBegin(g,vd1)+fragmentLength(g,vd1))-offset1,
		       id2,
		       vd2len-(fragmentBegin(g,vd2)+fragmentLength(g,vd2))-offset2, 
		       fragmentLength(g,vd1),
		       false);
    }
  }
  else{
    if(orient2 == '-'){
      //id1:+ id2:-
      assert(orient1 == '+');
      assert(orient2 == '-');
      frag = TFragment(id1,
		       fragmentBegin(g,vd1)-offset1,
		       id2,
		       vd2len-(fragmentBegin(g,vd2)+fragmentLength(g,vd1))-offset2, 
		       fragmentLength(g,vd1),
		       true);
    }
    else{
      //id1:+ id2:+
      assert(orient1 == '+');
      assert(orient2 == '+');
      frag = TFragment(id1,
		       fragmentBegin(g,vd1)-offset1,
		       id2,
		       fragmentBegin(g,vd2)-offset2, 
		       fragmentLength(g,vd1),
		       false);
    }
  }
  assert(frag.begin1 >=0);
  assert(frag.begin2 >=0);
  assert(frag.begin1+frag.len <=vd1len);
  assert(frag.begin2+frag.len <=vd2len);
}
	    
//Matches and fragment graph both have coordinates on the leading strand
template<typename TGraph,
	 typename TString,
	 typename TSpec,
	 typename TFragmentString,
	 typename TScoreValues>
void buildMatchesFromGraph(TGraph &g,
			   StringSet<TString, TSpec> &seqSet,
			   TFragmentString &currmatches,
			   TScoreValues &currscores){
  typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
  typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor; 
  typedef typename Size<StringSet<TString, TSpec> >::Type TSize;
  typedef Fragment<> TFragment;

  typedef typename Id<TGraph>::Type TId;
  TVertexDescriptor nilVertex = getNil<TVertexDescriptor>();
  typedef typename Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
  TEdgeIterator itE(g);
  TVertexDescriptor vd1,vd2;
  TSize vd1seq,vd2seq;
  for(;!atEnd(itE);goNext(itE)){
    vd1 = getSource(*itE);
    vd2 = getTarget(*itE);
    vd1seq = sequenceId(g,vd1);
    vd2seq = sequenceId(g,vd2);
    assert(vd1!=nilVertex);
    assert(vd2!=nilVertex);
    assert(vd1!=vd2);
    //findEdge() implemented in graph_impl_undirected.h 
    //traced from data_align object in graph_impl_align.h
    TEdgeDescriptor ed = findEdge(g,vd1,vd2);
    assert(ed);
    //There is an alignment between vd1 and vd2
    assert(vd1seq!=vd2seq);
    int vd1len = length(getValueById(stringSet(g), vd1seq));
    int vd2len = length(getValueById(stringSet(g), vd2seq));
    TFragment currfrag;	
    buildFrag(currfrag,g,vd1,vd2,
	      vd1seq,vd1len,0,'+',
	      vd2seq,vd2len,0,(int)(cargo(ed)<0) ? '-' : '+');
    assert(currfrag.begin1 >=0);
    assert(currfrag.begin2 >=0);
    assert(currfrag.begin1+currfrag.len <=length(getValueById(stringSet(g), vd1seq)));
    assert(currfrag.begin2+currfrag.len <=length(getValueById(stringSet(g), vd2seq)));
    assert(currfrag.len==fragmentLength(g,vd1));
    assert(currfrag.len==fragmentLength(g,vd2));
    if(currfrag.reversed){
      currfrag.begin2 = length(seqSet[currfrag.seqId2]) - (currfrag.begin2+currfrag.len);
    }
    appendValue(currmatches, currfrag);
    appendValue(currscores, fragmentLength(g,vd1));
  }
}

//buildMatchesFromGraph()
//Populate a TFragmentString(currmatches)
//based on edges in TGraph(g) that connect vertices in TLCB(lit)
//Only vertices with sequenceId in seqIdMap are considered

//TODO: Current impl iterates over all pairs of vertices. Performance
//can be improved using a BFS
template<typename TGraph,
	 typename TString,
	 typename TSpec,
	 typename TSeqLenMap,
	 typename TLCB,
	 typename TSize,
	 typename TOffsets,
	 typename TFragmentString,
	 typename TScoreValues>
void buildMatchesFromGraph(TGraph &g,
			   StringSet<TString, TSpec> &seqSet,
			   TSeqLenMap &seqLenMap,
			   TLCB &lit,
			   std::map<TSize,TOffsets> &offsets,
			   std::map<TSize,TSize> &seqIdMap,
			   TFragmentString &currmatches,
			   TScoreValues &currscores){
  typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
  typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
  typename TLCB::const_iterator vit;
  typename TLCB::const_iterator vit2;
  typename TLCB::const_iterator vit_end;
  typename TLCB::const_iterator vit2_end;

  typedef Fragment<> TFragment;

  typedef typename Id<TGraph>::Type TId;
#ifdef NDEBUG
  ;
#else
  TVertexDescriptor nilVertex = getNil<TVertexDescriptor>();
#endif
  //if edge(vit,vit2) present in input graph
  //then mark as present in output graph
  //Add vertex from graph g to LCB graph
  vit_end = lit.end();
  vit2_end = lit.end();
#ifdef DEBUGGING
  std::cout << "Building matches from graph for LCB" << std::endl;
#endif
  TVertexDescriptor vd1,vd2;
  TSize vd1seq,vd2seq;
  TEdgeDescriptor ed;
  for(vit = lit.begin();vit!=vit_end;vit++){
    vd1 = *vit;
    vd1seq = sequenceId(g,vd1);
    for(vit2 = lit.begin();vit2!=vit2_end;vit2++){
      //TODO see if shortcircuit on id1!=id2 improves performance here
      if(vit != vit2){
	assert(vd1!=nilVertex);
	//assert(degree(g,vd1)>0);
	vd2 = *vit2;
	vd2seq = sequenceId(g,vd2);
	assert(vd2!=nilVertex);
	//assert(degree(g,vd2)>0);
	assert(vd1!=vd2);
	//findEdge() implemented in graph_impl_undirected.h 
	//traced from data_align object in graph_impl_align.h
	ed = findEdge(g,vd1,vd2);
	if(!ed){

	}
	//There is an alignment between vd1 and vd2
	if(ed){
	  assert(vd1seq!=vd2seq);
	  TId id1 = idToPosition(seqSet, vd1seq);	  
	  TId id2 = idToPosition(seqSet, vd2seq);
	  assert(id1!=id2);
	  if(seqIdMap.find(id1)!=seqIdMap.end() &&
	     seqIdMap.find(id2)!=seqIdMap.end()){ //Check edge and if sequence was not trimmed out of LCB due to length
	    assert(seqIdMap.find(id1)!=seqIdMap.end());
	    assert(seqIdMap.find(id2)!=seqIdMap.end());
	    
	    int vd1len = seqLenMap[vd1seq];
	    int vd2len = seqLenMap[vd2seq];
	    //assert(vd1len == length(getValueById(stringSet(g), vd1seq)));
	    //assert(vd2len == length(getValueById(stringSet(g), vd2seq)));
	    //length(getValueById(stringSet(g), idToPosition(seqSet,vd1seq)))	    //length(getValueById(stringSet(g),id1));
	    //length(getValueById(stringSet(g), idToPosition(seqSet,vd2seq))) 	    //length(getValueById(stringSet(g),id2));

#ifdef DEBUGGING
	    std::cout << " seqs:" 
		      << seqIdMap[vd1seq] << ":" << offsets[id1].orient
		      << " " 
		      << seqIdMap[vd2seq] << ":" << offsets[id2].orient
		      << " lengths:" << vd1len << " "  << vd2len
		      << " coords:" 
		      << fragmentBegin(g,vd1) << "-"  << fragmentBegin(g,vd1) + fragmentLength(g,vd1)
		      << "," 
		      << fragmentBegin(g,vd2) << "-"  << fragmentBegin(g,vd2) + fragmentLength(g,vd1)
		      << " offset1 " << offsets[id1].offset
		      << " offset2 " << offsets[id2].offset
		      << " edge weight:" << cargo(ed)
		      << std::endl;
#endif
	    

	    if(offsets[id1].orient == '-'){
	      if(offsets[id2].orient=='+'){
		//id1:- id2:+
		assert(offsets[id1].orient == '-');
		assert(offsets[id2].orient == '+');
		appendValue(currmatches, TFragment(seqIdMap[id1],
						   vd1len-(fragmentBegin(g,vd1)+fragmentLength(g,vd1))-offsets[id1].offset,
						   seqIdMap[id2],
						   fragmentBegin(g,vd2)-offsets[id2].offset, 
						   fragmentLength(g,vd1),
						   false));
	      }
	      else{
		//id1:- id2:-
		assert(offsets[id1].orient == '-');
		assert(offsets[id2].orient == '-');
		appendValue(currmatches, TFragment(seqIdMap[id1],
						   vd1len-(fragmentBegin(g,vd1)+fragmentLength(g,vd1))-offsets[id1].offset,
						   seqIdMap[id2],
						   vd2len-(fragmentBegin(g,vd2)+fragmentLength(g,vd2))-offsets[id2].offset, 
						   fragmentLength(g,vd1),
						   false));
	      }
	      appendValue(currscores, fragmentLength(g,vd1));
	    }
	    else{
	      if(offsets[id2].orient == '-'){
		//id1:+ id2:-
		assert(offsets[id1].orient == '+');
		assert(offsets[id2].orient == '-');
		appendValue(currmatches, TFragment(seqIdMap[id1],
						   fragmentBegin(g,vd1)-offsets[id1].offset,
						   seqIdMap[id2],
						   vd2len-(fragmentBegin(g,vd2)+fragmentLength(g,vd1))-offsets[id2].offset, 
						   fragmentLength(g,vd1),
						   false));
	      }
	      else{
		//id1:+ id2:+
		assert(offsets[id1].orient == '+');
		assert(offsets[id2].orient == '+');
		appendValue(currmatches, TFragment(seqIdMap[id1],
						   fragmentBegin(g,vd1)-offsets[id1].offset,
						   seqIdMap[id2],
						   fragmentBegin(g,vd2)-offsets[id2].offset, 
						   fragmentLength(g,vd1),
						   false));
	      }
	      appendValue(currscores, fragmentLength(g,vd1));
	    }	    
	  }
	  else{
	    //Ignore matches, one of the sequences was trimmed from the LCB, probably due to length
#ifdef DEBUGGING
	    if(seqIdMap.find(id1)==seqIdMap.end()){
	      std::cout << "Ignoring match. Trimmed seq V" << vd1 << "-V" << vd2 << " S1" <<id1 << std::endl;
	    }
	    if(seqIdMap.find(id2)==seqIdMap.end()){
	      std::cout << "Ignoring match. Trimmed seq V" << vd1 << "-V" << vd2 << " S2" <<id2 << std::endl;
	    }
#endif
	  }
	}
	else{
	  //Ignore matches
	  //Segments vd1,vd2 are connected in the component
	  //but not directly via an alignment edge
#ifdef DEBUGGING
	  TId id1 = idToPosition(seqSet, vd1seq);	  
	  TId id2 = idToPosition(seqSet, vd2seq);
	  if(id1 != id2
	     && seqIdMap.find(id1)!=seqIdMap.end() 
	     && seqIdMap.find(id2)!=seqIdMap.end()){
	    assert(seqIdMap.find(id1)!=seqIdMap.end());
	    assert(seqIdMap.find(id2)!=seqIdMap.end());
	    /*
	    std::cout << "Ignoring match. Indirect connections V" << vd1 << "-V" << vd2 << " S1" <<id1 << "-" << " S2" << id2 << std::endl;
	    std::cout << "Ignored seqs:" 
		      << seqIdMap[vd1seq] << ":" << offsets[id1].orient
		      << " " 
		      << seqIdMap[vd2seq] << ":" << offsets[id2].orient
		      << " coords:" 
		      << fragmentBegin(g,vd1) << "-"  << fragmentBegin(g,vd1) + fragmentLength(g,vd1)
		      << "," 
		      << fragmentBegin(g,vd2) << "-"  << fragmentBegin(g,vd2) + fragmentLength(g,vd1)
		      << " offset1 " << offsets[id1].offset
		      << " offset2 " << offsets[id2].offset
		      << std::endl;
	    */
	  }
#endif
	}
      }
    }
  }
}


//Set currseqs,offsets
template<typename TSeqID,
	 typename TString,
	 typename TSpec,
	 typename TMap,
	 typename TOffsets>
void setLCBProps(TSeqID i,
		 char currorient,
		 unsigned int min,
		 unsigned int max,
		 TString & str,
		 StringSet<TString, TSpec> &seqSet,
		 TMap & seqIdMap,
		 TOffsets &offsets){
  //Check if mapped ids
  TSeqID seqidx = i;
  if(seqIdMap.find(i)!=seqIdMap.end()){
    seqidx = seqIdMap[i];
  }
  //Set substring in currseqs
  //and offsets
  //First, check orientation
  if(currorient == '+'){
    str = infix(seqSet[seqidx],min,max);
    //appendValue(currseqs,infix(seqSet[i],min,max));
    //orients[i] = '+';
    offsets[i].orient = '+';
    
    //offsets[i] = min;
    offsets[i].offset = min;
    
    //spanlens[i] = max-min;
    offsets[i].spanlen = max-min;
    
    //seqlens[i] = length(seqSet[i]);
    offsets[i].seqlen = length(seqSet[seqidx]);
  }
  else{
    assert(currorient == '-');
    //Handle reverse orientation
#ifdef DEBUGGING
    std::cout << "REVERSING SEQUENCE " << i << " of length " << length(seqSet[seqidx]) 
	      << " " << min << " - " << max << " " << std::endl;
#endif
    assert(min<max);
    assert(int(min)>=0);
    assert(max<=length(seqSet[seqidx]));
    //TString str = DnaStringReverseComplement(infix(seqSet[i],min,max));
    //TString str = infix(seqSet[i],min,max);
    str = infix(seqSet[seqidx],min,max);
    //Reverse complement
    convertInPlace(str, FunctorComplement<Dna5>());
    reverseInPlace(str);
    //appendValue(currseqs,str);
    //Now relative to - strand
    //Offsets assumed relative to matching strand
    //Also MAF stores coordinates relative to matching strand
    int tmpmin = min;
    min = length(seqSet[seqidx]) - max;
    max = length(seqSet[seqidx]) - tmpmin;
    offsets[i].orient = '-';
    offsets[i].offset = min;
    offsets[i].spanlen = max-min;
    offsets[i].seqlen = length(seqSet[seqidx]);
  }

}


//Determines orienation for sequence i in LCB lit
//Uses the vertexOrientMap that was previously built using relative orientation
//of most consistent matches in the original alignment graph
template<typename TGraph,
	 typename TString,
	 typename TSpec,
	 typename TLCB,
	 typename SeqID,
	 typename TVertexOrientMap>
bool getLCBProps(TGraph &g,
		 StringSet<TString, TSpec> &seqSet,
		 TLCB &lit,
		 SeqID i,
		 char & currorient, 
		 unsigned int &min,
		 unsigned int &max,
		 unsigned int &alnlen,
		 TVertexOrientMap &vertexOrientMap){
  typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
  typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
  std::vector<unsigned int>::const_iterator vit;
  std::vector<unsigned int>::const_iterator vit_end;
  bool seqPresent=false;
  bool resetOrientMajorityRule=false;
  if(resetOrientMajorityRule){
    //Resolve conflicts in orientation, use a majority rule to assign block orientation
    //TODO, determine what best to do with misoriented, conflicting blocks
    int plusorient=0;
    int minusorient=0;
    vit_end = lit->end();
    for(vit = lit->begin();vit!=vit_end;++vit){
      TVertexDescriptor vd1 = *vit;
      if(idToPosition(seqSet, sequenceId(g,vd1))==i){
	seqPresent=true;
	if(vertexOrientMap[vd1] == '+'){
	  plusorient++;
	}
	else{
	  assert(vertexOrientMap[vd1] == '-');
	  minusorient++;
	}
      }
    }
    if(plusorient>=minusorient){
      currorient = '+';
    }
    else{
      currorient = '-';
    }
  }
  else{
    //Orient already set in vertexOrientMap
  }
  

  //Currently, this method will use the first encountered
  //orientation for sequence $i in LCB $lit
  vit_end = lit->end();
  for(vit = lit->begin();vit!=vit_end;++vit){
    TVertexDescriptor vd1 = *vit;
    if(idToPosition(seqSet, sequenceId(g,vd1))==i){
      seqPresent=true;
      //Determine orientation for the block
      assert(vertexOrientMap.find(vd1) != vertexOrientMap.end());
      assert(vertexOrientMap[vd1] != '?');
      if(currorient != '?'){
	//All vertices in a block should have the same orientation
	if(vertexOrientMap[vd1] != currorient){
	  //There is a conflict
#ifdef DEBUGGING
	  std::cout << "Conflicting orientation on seq:" << i << " currorient:" << currorient << " V:" << vd1 << " expecting:" << vertexOrientMap[vd1] << std::endl;
#endif
	  if(resetOrientMajorityRule){
	    vertexOrientMap[vd1] = currorient;
	  }
	  else{
	    //if(msaOpt.refine == "colinear"){
	    //assert(false);
	    //}
	  }
	}
      }
      else{
	if(resetOrientMajorityRule){
	  assert(false);//should be using majority rule code now
	}
	else{
	  currorient = vertexOrientMap[vd1];
	}
      }
      assert(currorient != '?');
#ifdef DEBUGGING
      std::cout << "Determining orient for LCB " << " seq:" << i << " V:" << vd1 << " orient:" << currorient 
		<< " block:" << std::endl;
      std::cout << "Coords:" << fragmentBegin(g,vd1) << "-" << fragmentBegin(g,vd1)+fragmentLength(g,vd1) 
		<< " min:" << min << " max:" << max << std::endl;
#endif
      //Min max are always on the leading strand here
      assert((int)fragmentBegin(g,vd1)>=0);
      alnlen = alnlen + fragmentLength(g,vd1);
      min = (fragmentBegin(g,vd1)<min) ? fragmentBegin(g,vd1) : min;
      max = (fragmentBegin(g,vd1)+fragmentLength(g,vd1)>max) ? fragmentBegin(g,vd1)+ fragmentLength(g,vd1) : max;
    }
  }
    
  return seqPresent;
}


//retrieveLCBSegments()
//Populate TFragmentString(currmatches) for TLCB(currlcb) using TGraph(g)

//TODO:retrieve LCBSegments from initial match set (TFragmentString matches)
//rather than the alignment graph
template<typename TGraph, 
	 typename TString,
	 typename TSpec,
	 typename TString2,
	 typename TSpec2,
	 typename TMap,
	 typename TVertexOrientMap, 
	 typename TLCB, 
	 typename TNames, 
	 typename TSequence,
	 typename TFragmentString, 
	 typename TScoreValues,
	 typename TSize,
	 typename TOffsets,
	 typename TCoveredSet,
	 typename TSortedV>
void retrieveLCBSegments(TGraph & g, 
			 StringSet<TString, TSpec> &seqSetv,    //can include placeholder seqs, ie. no seq strings
			 StringSet<TString2, TSpec2> &seqSetReal, //can include virtual seqs, must have strings
			 TMap & seqIdxMap,                      //mapping between seqSetv->seqSetReal
			 TVertexOrientMap &vertexOrientMap,
			 TLCB &currlcb, 
			 unsigned int lcbid,
			 TNames & sequenceNames,
			 StringSet<TSequence, Owner<> > & currseqs,
			 TFragmentString &currmatches, 
			 TScoreValues &currscores,
			 TNames & currnameSet,
			 std::map<TSize,TOffsets> &offsets,
			 TCoveredSet & coveredSet,
			 TSortedV & vseqs,
			 unsigned int MIN_FRAGMENT_SIZE){
  assert(MIN_FRAGMENT_SIZE>=1);
  if(lcbid>0){};
  //Vars
  typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
  typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
  
  //Need to trim seqSet to members and subsequences present in the current LCB
  //Map to track old seqid to newseqid
  std::map<TSize,TSize> seqIdMap;
  std::map<TSize,TSize> seqLenMap;

  std::set<TVertexDescriptor> currlcbset;

  typename std::set<TVertexDescriptor>::iterator pos;

  typedef typename Id<TGraph>::Type TId;

  typename std::vector<TVertexDescriptor>::const_iterator vit;
  typename std::vector<TVertexDescriptor>::const_iterator vit_end;
  vit_end = currlcb->end();
  for(vit = currlcb->begin();vit!=vit_end;vit++){
    currlcbset.insert(*vit);
  }
  
  //Determine sequences present in LCB
  //and calculate spanning coords min,max
  for(TSize i = 0; i<length(seqSetv); ++i) {
    //Change to intmax
    unsigned int min=std::numeric_limits<unsigned int>::max();
    unsigned int max=0;
    unsigned int alnlen=0;
    char currorient = '?';
    //
    //Filter LCBs so that we only include sequences that span >
    //MIN_FRAGMENT_SIZE
    if(getLCBProps(g,seqSetv,currlcb,
		   i,currorient,min,max,alnlen,
		   vertexOrientMap)){
      if(alnlen>=MIN_FRAGMENT_SIZE &&
	 max-min>=MIN_FRAGMENT_SIZE){
#ifdef DEBUGGING
	assert(currorient!='?');
	assert(max>0); 
	assert((int)min>=0);
	assert(min<max);
	std::cout << "LCB:" << lcbid << " seq:" << i << " " 
		  << min << "-" << max << " spanlen:" << max-min 
		  << " alnlen:" << alnlen
		  << " orient: " << currorient << std::endl;
#endif
#ifdef NDEBUG
	;
#else
	unsigned int nSeq = length(currseqs);
#endif
	//Subsequence of lcb on seq $i 
	TString lcbseqstr;
	setLCBProps(i,currorient,min,max,lcbseqstr,
		    seqSetReal,
		    seqIdxMap,
		    offsets);
	//Save association between current seq $i
	//and position in $currseqs
	assert(length(lcbseqstr)==max-min);
	appendValue(currseqs,lcbseqstr);
	appendValue(currnameSet,sequenceNames[i]);
	assert(length(currnameSet)==length(currseqs));
	
	seqIdMap.insert(std::make_pair(i,length(currseqs)-1));

	if(seqIdxMap.find(i)!=seqIdxMap.end()){
	  seqLenMap.insert(std::make_pair(i,length(getValueById(stringSet(g),idToPosition(seqSetReal, seqIdxMap[i])))));
	}
	else{
	  seqLenMap.insert(std::make_pair(i,length(getValueById(stringSet(g),idToPosition(seqSetReal, i)))));
	}

	//Make sure that we have added one seq
	assert(length(currseqs)==nSeq+1);

	//Using sort vertices on seqs, capture any missing vertices that are spanned by the LCB
	//TODO, use findVertex and index to avoid looking through all vertices
	//
	//TVertexDescriptor act_knot = findVertex(ali_g,seq_id,begin_pos);
	if(i<vseqs.size()){
	  for(vit = vseqs[i].begin();vit!=vseqs[i].end();++vit){
	    //Check that vertex is not already aligned
	    assert(sequenceId(g,*vit)==i);
	    if(fragmentBegin(g,*vit)>=min){
	      if(fragmentBegin(g,*vit)<max){
		  currlcbset.insert(*vit);
		  coveredSet.insert(*vit);
	      }
	      else{
		//past max, we can stop looking
		break;
	      }
	    }
	  }
	}
	else{
	  //Vseqs not populated
	}
      }
      else{
	//Seq fragment is too short to include in LCB
      }
    }
    else{
      //Seq not present in LCB
      //This is ok
    }
    
  }
  assert(seqIdMap.size()==length(currseqs));
#ifdef SEQAN_PROFILE2
  std::cerr << "LCB Init of segments done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) 
	    << " seconds" << std::endl;
#endif
  
  
#ifdef SEQAN_PROFILE2
  std::cerr << "LCB Building segments" << std::endl;
#endif
  
  //All coordinates for fragments
  //must be relative to the orientation determined previously
  //offsets array is always relative to the matching strand
  typedef Fragment<> TFragment;
  
  buildMatchesFromGraph(g,
			seqSetv,
			seqLenMap,
			currlcbset,//currlcb,
			offsets,
			seqIdMap,
			currmatches,
			currscores);

  //TODO double check that matches all match expected values in lcbseqstr here

#ifdef SEQAN_PROFILE2
    std::cerr << "LCB Building segments done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) 
	      << " seconds" << std::endl;
#endif
}

//TODO retrieve LCBSegments from initial match set (TFragmentString matches)
//Rather than the alignment graph
template<typename TGraph, 
	 typename TString,
	 typename TSpec,
	 typename TVertexOrientMap, 
	 typename TLCB, 
	 typename TNames, 
	 typename TSequence,
	 typename TFragmentString, 
	 typename TScoreValues,
	 typename TSize,
	 typename TOffsets,
	 typename TCoveredSet,
	 typename TV>
void retrieveLCBSegments(TGraph & g, 
			 StringSet<TString, TSpec> &seqSet,
			 TVertexOrientMap &vertexOrientMap,
			 TLCB &currlcb, 
			 unsigned int lcbid,
			 TNames & sequenceNames,
			 StringSet<TSequence, Owner<> > & currseqs,
			 TFragmentString &currmatches, 
			 TScoreValues &currscores,
			 TNames & currnameSet,
			 std::map<TSize,TOffsets> &offsets,
			 TCoveredSet &coveredset,
			 TV & vseqs,
			 unsigned int MIN_FRAGMENT_SIZE){
  //dummy empty map
  std::map<TSize,TSize> seqIdxMap;
  retrieveLCBSegments(g,
		      seqSet,
		      seqSet,
		      seqIdxMap,
		      vertexOrientMap,
		      currlcb,lcbid,
		      sequenceNames,
		      currseqs,
		      currmatches,
		      currscores,
		      currnameSet,
		      offsets,
		      coveredset,
		      vseqs,
		      MIN_FRAGMENT_SIZE);
}

void transformMAF(const char * maffile, 
		  FILE * outstrm, 
		  std::map<std::string,std::string> &currnameSet, 
		  std::map<std::string,unsigned int> &offsets, 
		  std::map<std::string,char> &orients, 
		  std::map<std::string,unsigned int> & seqlens){
  struct mafFile *mf;
#ifdef DEBUGGING
  std::cout << "Transforming maf file " << maffile << std::endl;
#endif
  mf = mafOpen(maffile, 0);
  struct mafAli *a, *A, *last_a;
  struct mafComp *c;
  A = last_a = NULL;
  while ((a = mafNext(mf)) != NULL) {
    if ((c = a->components) == NULL)
      assert(false);//fatal("empty maf entry");
    if (last_a == NULL)
      A = a;
    else
      last_a->next = a;
    last_a = a;
  }
  if(A==NULL){
#ifdef DEBUGGING
    std::cout << "can't find any alignments" << std::endl;
#endif
  }
  else{
    //Do transform
    char chrName[200], species_name[200];
    for (a = A; a != NULL; a = a->next) {
      int i=0;
      for(c=a->components; c!=NULL; c=c->next) {
	//Update coordinates
	parseSrcName(c->src, species_name, chrName);
	assert(currnameSet.find(std::string(chrName))!=currnameSet.end());
	//From UCSC FAQ The start of the aligning
	//region in the source sequence. This is a
	//zero-based number. If the strand field is
	//'-' then this is the start relative to the
	//reverse-complemented source sequence.
	
	//TODO 
	//Confirm reverse alignments during refine are not handled properly
	if(c->strand == '+'){
	  c->start = c->start+offsets[std::string(chrName)];
	}
	else{
	  //Must convert relative to matching strand from original match
	  
	  assert(c->strand == '-');
	  c->start = c->start+offsets[std::string(chrName)];
	}
	c->strand = orients[std::string(chrName)];
	c->src = (char *)currnameSet[std::string(chrName)].c_str();
	c->srcSize = seqlens[std::string(chrName)];
	i++; 
      }
      mafWrite(outstrm, a);
    }
    mafFileFree(&mf);
  }
}

template<typename TScore>
void runIterativeMUGSY(std::string & outputdir,
		       const std::string & fastafiles,
		       std::string & prefix,
		       const std::string & outprefix,
		       MsaOptions<Dna5 , TScore> const& msaOpt){

  char * mugsyinstall = std::getenv("MUGSY_INSTALL");
  std::ostringstream refinecmd;
  refinecmd << mugsyinstall
	    << "/mugsy " 
#ifdef DEBUGGING
	    << "--debug 5 --log refine.log "
#endif
	    << " --distance " << msaOpt.distance
	    << " --minlength 15"
	    << " --nucmeropts \"-l 10 -c 15\"" //relax matchlen
    //TODO consider removing this option during refine to allow for refinement of short blocks
	    << " --skipunique --directory " << outputdir 
	    << " --skiprefine --colinear "
	    << " --prefix " << prefix
	    << " " << fastafiles
	    << " 1>"
	    << outprefix << "mugsyrefine.stdout"
	    << " 2>" << outprefix << "mugsyrefine.stderr";
#ifdef DEBUGGING
    std::cout << refinecmd.str() << std::endl;
#endif
      int ret = system(refinecmd.str().c_str());
      if(ret!=0){
	std::cerr << refinecmd.str() << std::endl 
		  << "SYSTEM:" << ret << std::endl;
      }
      else{
#ifdef DEBUGGING
	;
#else
	std::string stdout(outprefix+"mugsyrefine.stdout");
	std::string stderr(outprefix+"mugsyrefine.stderr");
	std::string log(prefix+".mugsy.log");
	unlink(stdout.c_str());
	unlink(stderr.c_str());
	unlink(log.c_str());
#endif
      }
      assert(ret==0);
}

//Refinement using Mugsy
//Also support for fsa,pecan,lagan aligners. They must be in your path
//TODO, save label,dups from original MAF
template<typename TScore>
void refineMSA(const char * maffile,
	       MsaOptions<Dna5 , TScore> const& msaOpt){
  std::fstream strmmaf;
  FILE * strmmafrefined;

  struct mafFile *mf;
  mf = mafOpen(maffile, 0);
  struct mafAli *a, *A, *last_a;
  struct mafComp *c;
  A = last_a = NULL;
  while ((a = mafNext(mf)) != NULL) {
    if ((c = a->components) == NULL)
      assert(false);//fatal("empty maf entry");
    if (last_a == NULL)
      A = a;
    else
      last_a->next = a;
    last_a = a;
  }
  if(A==NULL){
#ifdef DEBUGGING
    std::cout << "can't find any alignments" << std::endl;
#endif
  }
  else{
    std::string outfile(msaOpt.outfile);

    std::vector<char> writable(outfile.size() + 1);
    std::copy(outfile.begin(), outfile.end(), writable.begin());
    std::string outputdir(dirname(&writable[0]));
    if(outputdir[outputdir.length()-1] != '/'){
      outputdir = outputdir + '/';
    }

    strmmafrefined = fopen(std::string(outfile+".maf.refined").c_str(),"w");//, std::ios_base::out | std::ios_base::trunc);
    mafWriteStart(strmmafrefined, "mugsy_refined");
    //Do transform
    char chrName[200], species_name[200];
    int lcbid=0;
    int COL_WIDTH=60;
    for (a = A; a != NULL; a = a->next) {
      std::map<std::string,unsigned int> curroffsets;
      std::map<std::string,unsigned int> currspanlens;
      std::map<std::string,unsigned int> currseqlens;
      std::map<std::string,char> currorients;
      std::map<std::string,std::string> currnameSetv;
      
      int ncol = a->textSize;
      int i=0;
      std::ostringstream tmpgraph;
      tmpgraph << "MUGTMP" << getpid() << "_" << lcbid;

      std::vector<std::string> fnames;
      for(c=a->components; c!=NULL; c=c->next) {
	std::fstream strmfsa;
	std::string fname(tmpgraph.str());
	fname = outputdir+fname + "_S"+boost::lexical_cast<std::string>(i) + ".fsa";
	strmfsa.open(fname.c_str(), std::ios_base::out | std::ios_base::trunc);
	fnames.push_back(fname);
	parseSrcName(c->src, species_name, chrName);
	//Write FASTA
	if(msaOpt.refine=="fsa"){
	  //Write XMFA style
	  strmfsa << ">" << tmpgraph.str() << "_S" <<boost::lexical_cast<std::string>(i) 
		  << "." << c->src << ":" << 1 << "-" << c->size << " " << c->strand << " " << c->size << std::endl;
	}
	else{
	  strmfsa << ">" << c->src << std::endl ;
	}
	int col=0;
	int j=0;
	for (col = j = 0; j < ncol; ++j) {
	  if(c->text[j]=='-'){

	  }
	  else{
	    strmfsa << c->text[j];
	    ++col;
	    if (col == COL_WIDTH) {
	      strmfsa << std::endl;
	      col = 0;
	    }
	  }
	}
	if (col != 0){
	  strmfsa << std::endl;
	}
	std::string sname(c->src);
	curroffsets[sname] = c->start;
	currspanlens[sname] = c->size;
	currseqlens[sname] = c->srcSize;
	currorients[sname] = c->strand;
	currnameSetv[sname] = sname;
	++i;
	strmfsa.close();
      }
      
      //Require more than one sequence
      if(fnames.size()>1){
	//
	std::ostringstream fastafiles;
	for(int k=0;k<(int)fnames.size();k++){
	  fastafiles << fnames[k] << " " ;
	}
	
	
	//Output MAF file
	std::string prefix("MGREF");
	if(a->label>=0){
	  prefix = prefix + boost::lexical_cast<std::string>(a->label);
	}
	std::string maffile(outputdir+"/"+prefix+".maf");
	//Clean up old maf with same name
	unlink(maffile.c_str());
	
	//Run refinement 
	//Support for other aligners is provided for evaluation
	if(msaOpt.refine=="pecan"){
	  //Support for pecan aligner
	  std::ostringstream treecmd;
	  //treecmd << "cat " << fastafiles.str() << " | /usr/local/projects/angiuoli/developer/sangiuoli/muscle/trunk/muscle -clusteronly -in - -tree1 /tmp/pecan.tree 1> /dev/null 2> /dev/null";
	  treecmd << "cat " << fastafiles.str() << " | muscle -clusteronly -in - -tree1 /tmp/pecan.tree 1> /dev/null 2> /dev/null";
	  int ret = system(treecmd.str().c_str());
	  if(ret!=0){
	    std::cerr << treecmd.str() << std::endl 
		      << "SYSTEM:" << ret << std::endl;
	  }
	  
	  std::ostringstream refinecmd;
	  //refinecmd << "java -cp /usr/local/projects/angiuoli/developer/sangiuoli/pecan_v0.8/pecan_v0.8.jar bp.pecan.Pecan -J /usr/local/projects/angiuoli/developer/sangiuoli/exonerate-2.2.0-x86_64/bin/exonerate -E `cat /tmp/pecan.tree | perl -ne 'chomp;print'` -F " << fastafiles.str() << " >> pecan." << getpid() << ".mfa";
	  refinecmd << "java -cp pecan_v0.8.jar bp.pecan.Pecan -J exonerate -E `cat /tmp/pecan.tree | perl -ne 'chomp;print'` -F " << fastafiles.str() << " >> pecan." << getpid() << ".mfa";
	  ret = system(refinecmd.str().c_str());
	  if(ret!=0){
	    std::cerr << refinecmd.str() << std::endl 
		      << "SYSTEM:" << ret << std::endl;
	  }
	  else{
	  }
	}
	else if(msaOpt.refine == "mlagan"){
	  std::ostringstream refinecmd;
	  //refinecmd << "/usr/local/projects/angiuoli/developer/sangiuoli/lagan20/mlagan.sh " << fastafiles.str() << " >> lagan."<< getpid() << ".mfa 2> /dev/null";
	  refinecmd << "mlagan.sh " << fastafiles.str() << " >> lagan."<< getpid() << ".mfa 2> /dev/null";
	  int ret = system(refinecmd.str().c_str());
	  if(ret!=0){
	    std::cerr << refinecmd.str() << std::endl 
		      << "SYSTEM:" << ret << std::endl;
	  }
	}
	else if(msaOpt.refine == "fsa"){
	  ostringstream refinecmd;
	  //refinecmd << "/usr/local/projects/angiuoli/developer/sangiuoli/fsa-1.15.3/src/main/fsa --fast --noindel2 --refinement 0 " << fastafiles.str() << " > fsa." << getpid() << ".mfa 2>> test.fsa.stderr";
	  //refinecmd << "fsa --anchored --maxram 15000 --fast --noindel2 --refinement 0 " << fastafiles.str() << " > fsa." << getpid() << ".mfa 2>> test.fsa.stderr";
	  refinecmd << "fsa --fast --noindel2 --refinement 0 " << fastafiles.str() << " > fsa." << getpid() << ".mfa 2>> test.fsa.stderr";
	  int ret = system(refinecmd.str().c_str());
	  if(ret!=0){
	    std::cerr << refinecmd.str() << std::endl 
		      << "SYSTEM:" << ret << std::endl;
	  }
	  ostringstream convertcmd;
	  convertcmd << "echo '=' >> fsa." << getpid() << ".mfa;" 
		     << std::getenv("MUGSY_INSTALL") << "/xmfa2maf.pl < fsa." << getpid() << ".mfa > " << maffile.c_str() << " 2>> test.maf.stderr";
	  ret = system(convertcmd.str().c_str());
	  if(ret!=0){
	    std::cerr << convertcmd.str() << std::endl 
		      << "SYSTEM:" << ret << std::endl;
	  }
	}
	else{
	  runIterativeMUGSY(outputdir,fastafiles.str(),prefix,tmpgraph.str(),msaOpt);
	}
	std::cerr << ".";
	//Need to clean up here to prevent huge proliferation of files
#ifdef DEBUGGING
	;
#else
	for(int k=0;k<(int)fnames.size();k++){
	  unlink(fnames[k].c_str());
	}
#endif
	//Library call added to multiz for parsing
	FILE* intFileDescriptor ;
	struct stat stat_FileStatistics ;
	
	intFileDescriptor = fopen(maffile.c_str(), "r");
	if(intFileDescriptor != NULL){
	  fstat(fileno(intFileDescriptor), &stat_FileStatistics) ;
	  unsigned long size = stat_FileStatistics.st_size ;
	  fclose(intFileDescriptor);
	  if(size>0){
	    assert(currnameSetv.size()==fnames.size()); 
	    transformMAF(maffile.c_str(),
			 strmmafrefined,
			 currnameSetv,
			 curroffsets,
			 currorients,
			 currseqlens);
#ifdef DEBUGGING
	    ;
#else
	    unlink(maffile.c_str());
#endif
	  }
	  else{
	    //Refined MAF file has zero length
	  }
	}
	++lcbid;
      }
      else{
	mafWrite(strmmafrefined, a);      
      }
      //mafAliFree(&a); 
    }
  }
}
	
template<typename TStringSet,
	 typename TCargo,
	 typename TSpec,
	 typename TLCB,
	 typename TStringSet1,
	 typename TNames,
	 typename TGenomeNames,
	 typename TVertexOrientMap,
	 typename TIntervals,
	 typename TScore>
void generateLCBs(Graph<Alignment<TStringSet, TCargo, TSpec> > &g,
		  TLCB &LCBs,
		  TStringSet1 &seqSet,
		  TNames &sequenceNames,
		  TGenomeNames &genomeNames,
		  TVertexOrientMap &vertexOrientMap,
		  TIntervals & aintervals,
		  MsaOptions<Dna5 , TScore> const& msaOpt){
  //Configurable options
  bool useadjscores=false; //Generate and use adjacency scores  
  bool bpanalysis=(msaOpt.segmentation=="none") ? false : true;    //Set to false to skip breakpoint analysis entirely, each CC in segment graph will be an LCB
  //*******
  //Retrieve initial set of alignment blocks(LCBs) from segment graph
  //A block is a set of segments that are connected in the segment graph
#ifdef SEQAN_PROFILE
  std::cerr << "Converting segments to multi-genome anchors " 
	    << length(seqSet) << " " 
	    << length(sequenceNames) << " " 
	    << length(genomeNames) << " " 
	    << vertexOrientMap.size() << " " 
	    << LCBs.size() << " " 
	    << numVertices(g) << " " 
	    << numEdges(g) << std::endl;
#endif
  typedef Dna5 TAlphabet;
  typedef typename Value<TScore>::Type TScoreValue;
  typedef typename Size<TStringSet>::Type TSize;
  typedef typename Value<TStringSet1>::Type TString;
  typedef typename Value<TNames>::Type TName;
  //typedef Graph<Alignment<TStringSet, TSize> > TGraph;
  //Using int to support negative edge scores
  typedef Graph<Alignment<TStringSet, int> > TGraph;
  typedef typename Id<TGraph>::Type TId; 
  typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
  typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
  typedef typename EdgeType<TGraph>::Type TEdgeStump;
  typedef typename Iterator<String<TEdgeStump*> const, Rooted>::Type TIterConst;
  typedef typename Iterator<String<TEdgeStump*>, Rooted>::Type TIter;
  
  //
  typedef std::map<unsigned int, unsigned int> TComponentLength;
  
  // Strongly Connected Components, topological sort, and length of each component
  typedef String<unsigned int> TComponentMap;
  typedef typename Value<TComponentMap>::Type TComponent;
  typedef typename Position<TGraph>::Type TPos;
  typedef SVABlock<TComponent,TSize,TVertexDescriptor,TPos> TBlock;
  
  TComponentMap component;
  typedef typename Value<TComponentMap>::Type TComponent;
  //Hold input blocks that will be used to generate LCBs
  
  std::map<std::pair<TComponent,TComponent>,TBlock *> componentVertexMap;
  std::vector<std::vector<TBlock> > blocksbycomponent; 

  TSize numComponents;

  //Greedy algorithm for resolving conflicts in connecting segments
  //into blocks. Conflicts arise when there are more than 2 segments
  //connected from the same genome seperated by < msaOpt.poscombinewindow

  //Considering two methods
  //(1)Connect using best positional score first.
  //Derive positional score from an intial clustering
  //Break gaps that violate constraints using a mincut
  //(2)Connect using best consistency score
  //Start new cluster whenever a repeat/dup is to be added

  std::map<TEdgeDescriptor,float> posScores;
  if(useadjscores){
    //Adjacency scoring is optional, off by default
    //(1) Using adjacency and consistency score
    std::cerr << "Not implemented" << std::endl;
    exit(1);
    /*
    numComponents = convertSegments2BlocksAdjacency(g,
						    component,
						    componentVertexMap,
						    blocksbycomponent,
						    seqSet,
						    genomeNames,
						    posScores,
						    msaOpt,
						    cuts);
    convertCC2Blocks(g,
		     component,
		     componentVertexMap,
		     blocksbycomponent,
		     aintervals,
		     sequenceNames);
    */

  }
  else{
    //(2) Using consistency score
    std::cerr << "Greedy CC on consistency score " << std::endl;
    //numComponents = connected_components_by_genome_ranked_RECURSIVE(g, component, genomeNames, 100);
    //std::cerr << "numc ranked recur" << numComponents << std::endl;
    //numComponents = connected_components(g,component);
    //std::cerr << "numc reg" << numComponents << std::endl;
 
    //Convert segment graph (V=genome segments on one genome) into
    //anchor graph (V=genome segments on multiple genomes)
    numComponents = connected_components_by_genome_ranked(g, component, genomeNames, msaOpt.anchorwin);
    std::cerr << "Num components:" << numComponents << std::endl;
    blocksbycomponent.resize(numComponents);
    //Collapse CC into blocks
    convertCC2Blocks(g,
		     component,
		     componentVertexMap,
		     blocksbycomponent,
		     aintervals,
		     sequenceNames);
#ifdef DEBUGGING
    //Ensure there are not blocks with 2 seqs from the same genome
    typename std::vector<std::vector<TBlock> >::iterator bit2 = blocksbycomponent.begin();
    for(; bit2!= blocksbycomponent.end();bit2++){//all cc
      if(bit2->size()>0){
	for(unsigned int i=0;i<bit2->size();i++){
	  for(unsigned int j=i+1;j<bit2->size();j++){
	    for(typename std::vector<TVertexDescriptor>::iterator vit=bit2->at(i).currV.begin();vit!=bit2->at(i).currV.end();vit++){
	      TVertexDescriptor currV0 = *vit;
	      //assert(degree(g,currV0)>0);
	      for(typename std::vector<TVertexDescriptor>::iterator vit1=bit2->at(j).currV.begin();vit1!=bit2->at(j).currV.end();vit1++){
		TVertexDescriptor currV = *vit1;
		if(sequenceId(g,currV0)!=sequenceId(g,currV)){
		  std::cout << "V1:" << currV0 << " " << " V2:" << currV 
			    << " seq: " << sequenceId(g,currV0) << " " << sequenceId(g,currV) 
			    << " genome:" << genomeNames[sequenceId(g,currV0)] << "-" << genomeNames[sequenceId(g,currV)] << std::endl;
		  assert(genomeNames[sequenceId(g,currV0)]!=genomeNames[sequenceId(g,currV)]);
		}
	      }
	    }
	  }
	}
      }
    }
#endif
  }
  
  
  //
  //
  //Report some stats after building blocks
  unsigned totalavaillen=0;
  //unsigned totalmatchingbp=0;
  unsigned totalseqlen=0;
  typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
  TVertexIterator itV2(g);
  for(;!atEnd(itV2);goNext(itV2)){
    if(degree(g,*itV2)>0){
      totalavaillen+=fragmentLength(g,*itV2);
    }
  }
  //TVertexIterator itV(g);
  //for(;!atEnd(itV);goNext(itV)){
  //if(degree(g,*itV)>0){
  //totalmatchingbp+=fragmentLength(g,*itV);
  //}
  //}

  TSize seqSetLen = length(seqSet);
  for(unsigned int i=0;i<seqSetLen;i++){
   totalseqlen+=length(seqSet[i]);
  }

  //std::cerr << "Excluded unique, repeat/duplicated bp:" << totalmatchingbp-totalavaillen 
  //<< " " 
  //<< "=" << (float)(totalmatchingbp-totalavaillen)/totalseqlen << std::endl;
  std::cerr << "Percentage matching bp (not including matching repeats/dups):" 
	    << totalavaillen << "/" << totalseqlen
	    << "=" << (float)totalavaillen/totalseqlen << std::endl;
  //assert(totalavaillen<=totalmatchingbp);

#ifdef SEQAN_PROFILE
  std::cerr << "Anchor conversion done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;
  std::cerr << "Num anchors: " << numComponents << std::endl;
#endif	  
  //*********
  //Need to rescore because edges may have been removed
  //invalidating edge pointers
  posScores.clear();
  if(useadjscores){
    //Score for positional conservation
    std::cerr << "Rescoring for positional conservation" << std::endl;
    std::cerr << "Not implemented" << std::endl;
    exit(1);
    //scorePosCons(g,
    //	 component,
    //	 numComponents,
    //	 posScores,
    //	 msaOpt.posscorewindow);
  }
  //*********
  //Assign orientation to LCBs
#ifdef SEQAN_PROFILE
  std::cerr << "Assigning orientation to " << blocksbycomponent.size() << " anchors" << std::endl;
#endif
#ifdef DEBUGGING_GRAPH
  std::fstream rawstrm;
  rawstrm.open("refinegraphpreorient.out", std::ios_base::out | std::ios_base::trunc);
  write(rawstrm,g,seqSet,Raw());
  rawstrm.close();
#endif

  //Array of blocks
  std::vector<TBlock> blocks;
  //Assign orientation to blocks
  //+ reversed==false
  //- reversed==true
  //Save the orientation using vertexOrientMap 
  //Need a map so we can lookup orientation for each vertex in a block
  //std::map<TVertexDescriptor,char> vertexOrientMap;
  assert(vertexOrientMap.size()==0); //expecting empty map to start
  assignBlockOrientation(g,
			 blocksbycomponent,
			 blocks,
			 vertexOrientMap,
			 posScores); 
  posScores.clear();
  blocksbycomponent.clear();
  componentVertexMap.clear();
  clear(component);
  
#ifdef SEQAN_PROFILE
  std::cerr << "Orientation done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;
#endif
  
#ifdef SEQAN_PROFILE
  std::cerr << "Building an orthology map using " << blocks.size() << " anchors" << std::endl;
#endif
  //Determine collinear runs, ie LCBs
  //Save runs in the LCBs vector
  //Each LCB is a set of TVertexDescriptors
  //std::vector<std::vector<TVertexDescriptor> > LCBs;
  
  std::map<unsigned int, std::set<TVertexDescriptor> > block2fragMap;
  typename std::vector<TBlock>::const_iterator bit = blocks.begin();
  for(bit = blocks.begin();
      bit!=blocks.end();
      bit++){
#ifdef DEBUGGING
    typename std::vector<TVertexDescriptor>::const_iterator dvit;
    for(dvit = bit->currV.begin();dvit!=bit->currV.end();++dvit){ 
      assert((*dvit)!=getNil<TVertexDescriptor>());
    }
#endif
    if(block2fragMap.find(bit->c) == block2fragMap.end()){
      block2fragMap.insert(std::make_pair(bit->c,std::set<TVertexDescriptor>()));
      block2fragMap[bit->c].insert(bit->currV.begin(),bit->currV.end());
    }
    else{
      block2fragMap[bit->c].insert(bit->currV.begin(),bit->currV.end());
    }
  }

  if(!bpanalysis){
    //Assign each block as an LCB
    //Useful for some simple testing 
    typename std::map<unsigned int, std::set<TVertexDescriptor> >::iterator it;
    for(it = block2fragMap.begin();it!=block2fragMap.end();it++){
      std::vector<unsigned int> currlcb;
      currlcb.insert(currlcb.end(),it->second.begin(),it->second.end());
      LCBs.push_back(currlcb);
    }
  }
  else{ 
    //Read LCBs/blocks from input file if specified
    if(msaOpt.blockfile.length()>0){
      readBlockFile(msaOpt.blockfile,
		    block2fragMap,
		    LCBs,
		    sequenceNames,
		    vertexOrientMap,
		    g,
		    false);
    }
    else{
      //Calculate LCBs using segmentation method
      //Enredo, Mercator, or MUGSY
      std::string diststr(msaOpt.distance);
      std::string minlenstr(msaOpt.minlength);
      if(msaOpt.segmentation == "enredo"){
	do_segmentation_ENREDO(blocks,
			       LCBs,
			       block2fragMap,
			       diststr,
			       minlenstr,
			       msaOpt,
			       seqSet,
			       sequenceNames,
			       genomeNames,
			       vertexOrientMap,
			       g);
      }
      else{
	if(msaOpt.segmentation == "mercator"){
	  //TODO, fix this tested but now broken do_segmentation_MERCATOR()
	  std::cerr << "Mercator segmentation not implemented" << std::endl;
	}
	else{
	  if(msaOpt.refine == "colinear"){
	    //For refinement, assume all within a single LCB
	    typename std::map<unsigned int, std::set<TVertexDescriptor> >::iterator it;
	    std::vector<unsigned int> currlcb;
	    for(it = block2fragMap.begin();it!=block2fragMap.end();it++){
	      currlcb.insert(currlcb.end(),it->second.begin(),it->second.end());
	    }
	    LCBs.push_back(currlcb);
	  }
	  else{
	    //Default method is our own
	    do_segmentation_MUGSY(blocks,
				  LCBs,
				  block2fragMap,
				  diststr,
				  minlenstr,
				  msaOpt,
				  sequenceNames,
				  genomeNames,
				  vertexOrientMap,
				  g);
	  }
	}
      }
    }
  }
  blocks.clear();
#ifdef DEBUGGING
  //Print out LCBs
  for(typename std::vector<std::vector<unsigned int> >::const_iterator lit = LCBs.begin();lit!=LCBs.end();lit++){
    std::vector<unsigned int>::const_iterator vit;
    for(vit = lit->begin();vit!=lit->end();vit++){
      std::cout << *vit << ",";
      assert(*vit!=getNil<TVertexDescriptor>());//Null vertex
    }
    std::cout << std::endl;
  }	
#endif

#ifdef SEQAN_PROFILE
  std::cerr << "LCB identification done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;
#endif

}


template<typename TSequenceSet,
	 typename TIds,
	 typename TTreeMap,
	 typename TDistanceValue>
void getGuideTree(TSequenceSet &seqSet,
		  TIds &curridset,
		  TTreeMap &seqguideTrees,
		  Graph<Tree<TDistanceValue> > &currguideTree){
  typedef String<TDistanceValue> TDistanceMatrix;
  TDistanceMatrix distanceMatrix;

  typedef Dna5 TAlphabet;
  typedef unsigned TSize;
  typedef String<TAlphabet> TSequence;
  //Parse subtree
  TDistanceMatrix currdistanceMatrix;
  typedef typename Value<TDistanceMatrix>::Type TValue;
  typedef typename Iterator<TDistanceMatrix>::Type TMatrixIterator;
  
  //This assumes the best estimate is obtained by kmers across whole genome
  //Alternatively, build guide tree using subsequence of the LCBs for LCBs of length 
  //greater than some cutoff to avoid needlessly building trees for extremely short LCBs
  std::ostringstream curridsetstr;
  for(std::set<unsigned int>::iterator it = curridset.begin();it!=curridset.end();it++){
    curridsetstr << *it << ":";
  }
  std::string curridsetstring = curridsetstr.str();
  if(seqguideTrees.find(curridsetstring) == seqguideTrees.end()){
#ifdef DEBUGGING
    std::cout << "Generating a new guide tree for " << curridsetstring << std::endl;
#endif
    //Copy genome string for sequences present in the current LCB
    StringSet<TSequence, Owner<> > currStringSet;
    for(TSize i = 0; i<length(seqSet); ++i) {
      if(curridset.find(i) != curridset.end()){
	appendValue(currStringSet,seqSet[i]);
      }
    }
    
    //if (empty(distanceMatrix)) getDistanceMatrix(currG, currdistanceMatrix, KmerDistance());
    getKmerSimilarityMatrix(currStringSet, currdistanceMatrix, 3, TAlphabet());
    // Similarity to distance conversion
    TMatrixIterator matIt = begin(currdistanceMatrix);
    TMatrixIterator endMatIt = end(currdistanceMatrix);
    for(;matIt != endMatIt;++matIt) value(matIt) = (1 - (*matIt)) * 100;
    upgmaTree(currdistanceMatrix, currguideTree, UpgmaMin());
    //njTree(currdistanceMatrix, currguideTree);
    //Save tree for the combination specified by curridset
    seqguideTrees.insert(std::make_pair(curridsetstring,currguideTree));
    clear(currdistanceMatrix);
  }
  else{
#ifdef DEBUGGING
      std::cout << "Using existing guide tree for " << curridsetstring << std::endl;
#endif
      currguideTree=seqguideTrees[curridsetstring];
  }

  //  else{
  //use guidetree that covers all sequences
  //assert(numVertices(seqguideTree)>0);
  //currguideTree = seqguideTree;
  //assert(false);
  //}
#ifdef SEQAN_PROFILE2
    std::cout << "LCB Guide tree done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;
#endif
}

template<typename TGraph,
	 typename TGraph2,
	 typename TSize,
	 typename TStringSet,
	 typename TGuideTree,
	 typename TScore>
s_score alignSingleLCB(TGraph &currG,
			TGraph2 &currgOut,
			TSize lcbid,
			TStringSet &currseqs,
			TGuideTree &currguideTree,
			MsaOptions<Dna5 , TScore> const& msaOpt){
  typedef typename Value<TScore>::Type TScoreValue;
  typedef Dna5 TAlphabet;
  //Currently disabled
  bool inlinerefine= (msaOpt.refine=="true") ? true : false; //Compute iterative refinement inline
  if(lcbid>0){}
  if(inlinerefine){}
  TSize nSeq = length(currseqs);
  TSize threshold = 100;
#ifdef SEQAN_PROFILE2
  std::cout << "LCB Performing triplet extension " << lcbid << std::endl;
#endif
  if (nSeq < threshold) tripletLibraryExtension(currG);
  else tripletLibraryExtension(currG, currguideTree, threshold / 2);
#ifdef SEQAN_PROFILE2
  std::cout << "LCB Triplet extension done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;
#endif
  //
  //*******
  //Alignment
  //*******  
#ifdef SEQAN_PROFILE2
  std::cout << "LCB Performing progressive alignment" << std::endl;
#endif	
#ifdef DEBUGGING_GRAPH
  std::fstream rawstrm2;
  std::string lcbgname = "lcbgraph"+boost::lexical_cast<std::string>(lcbid)+".out";
  rawstrm2.open(lcbgname.c_str(), std::ios_base::out | std::ios_base::trunc);
  write(rawstrm2,currG,currseqs,Raw());
  rawstrm2.close();
#endif


  //Perform the alignment
  progressiveAlignment(currG, currguideTree, currgOut);


#ifdef SEQAN_PROFILE2
  std::cout << "LCB progressive alignment done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;
#endif  
  
  
#ifdef DEBUGGING_GRAPH
  std::fstream rawstrm3;
  std::string alcbgname = "alignedlcbgraph"+boost::lexical_cast<std::string>(lcbid)+".out";
  rawstrm3.open(alcbgname.c_str(), std::ios_base::out | std::ios_base::trunc);
  write(rawstrm3,currgOut,currseqs,Raw());
  rawstrm3.close();
#endif

  s_score sscore;
#ifdef SCORING
  // Print the scoring information
  //TScoreValue gop = msaOpt.sc.data_gap_open;
  //TScoreValue gex = msaOpt.sc.data_gap_extend;
  //TSize alphSize = ValueSize<TAlphabet>::VALUE;
  
  // Print the alignment information
  TSize numGapEx = 0;
  TSize numGap = 0;
  TSize numPairs = 0;
  TSize numIdents = 0;
  TSize alignLen = 0;
  TSize totalLen = 0;
  String<TSize> pairCount;
  String<TSize> colCount;
  //TScoreValue alignScore;
  sscore = alignmentEvaluationCustom(currgOut, 
			       msaOpt.sc, 
			       numGapEx, 
			       numGap, 
			       numPairs, 
			       numIdents,
			       pairCount, 
			       colCount,
			       alignLen, 
			       totalLen);
  /*
  sscore.alignScore = alignScore;
  sscore.numGap = numGap;
  sscore.numGapEx = numGapEx;
  sscore.numPairs = numPairs;
  sscore.numIdents = numIdents;
  sscore.alignLen = alignLen;
  sscore.totalLen = totalLen;
  sscore.colCount = colCount;
  sscore.seqCount = nSeq;
  assert(length(colCount)==nSeq+1);
  sscore.pairCount = pairCount;
  */
#endif
#ifdef DEBUGGING2
  TSize alphSize = ValueSize<TAlphabet>::VALUE;
  TScoreValue gop = msaOpt.sc.data_gap_open;
  TScoreValue gex = msaOpt.sc.data_gap_extend;
  std::cout << "LCBID:" << lcbid << std::endl;
  std::cout << "Scoring parameters:" << std::endl;
  std::cout << "*Gap opening: " << gop << std::endl;
  std::cout << "*Gap extension: " << gex << std::endl;
  std::cout << "*Scoring matrix: " << std::endl;
  std::cout << "   ";
  for(TSize col = 0; col<alphSize; ++col) std::cout << TAlphabet(col) << ',';
  std::cout << std::endl;
  for(TSize row = 0; row<alphSize; ++row) {
    for(TSize col = 0; col<alphSize; ++col) {
      if (col == 0) std::cout << TAlphabet(row) << ": ";
      //std::cout << score(scType, TAlphabet(row), TAlphabet(col));
      if (col < alphSize - 1) std::cout << ',';
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
  std::cout << "Alignment Score: " << sscore.alignScore << std::endl;
  std::cout << "Alignment Length: " << alignLen << std::endl;
  std::cout << "#Match-Mismatch pairs: " <<numPairs << std::endl;
  std::cout << "Score contribution by match-mismatch pairs: " << (sscore.alignScore - (((TScoreValue) numGap * gop) + ((TScoreValue) numGapEx * gex))) << std::endl;
  std::cout << "#Gap extensions: " << numGapEx << std::endl;
  std::cout << "Score contribution by gap extensions: " << ((TScoreValue) numGapEx * gex) << std::endl;
  std::cout << "#Gap openings: " << numGap << std::endl;
  std::cout << "Score contribution by gap openings: " << ((TScoreValue) numGap * gop) << std::endl;
  std::cout << std::endl;
  std::cout << "#Pairs: " << std::endl;
  std::cout << "   ";
  for(TSize col = 0; col<alphSize; ++col) std::cout << TAlphabet(col) << ',';
  std::cout << std::endl;
  for(TSize row = 0; row<alphSize; ++row) {
    for(TSize col = 0; col<alphSize; ++col) {
      if (col == 0) std::cout << TAlphabet(row) << ": ";
      std::cout << value(pairCount, row * alphSize + col);
      if (col < alphSize - 1) std::cout << ',';
    }
    std::cout << std::endl;
  }
#endif
  return sscore;
}



template<typename TName,
	 typename TLoc,
	 typename TNames>
void saveInterval(std::map<TName,std::vector<TLoc> >&aintervals,
		  TNames &currnameSet,
		  std::vector<unsigned int> &curroffsets,
		  std::vector<unsigned int> &currspanlens,
		  std::vector<unsigned int> &currseqlens,
		  std::vector<char> &currorients,
		  int lcbid,
		  bool dup=false){
  for(int i=0;i<(int)length(currnameSet);i++){
    TName n = currnameSet[i];
#ifdef DEBUGGING
    std::cout << "Saving aligned intervals for " << n << std::endl;
#endif
    typename std::map<TName,std::vector<TLoc > >::iterator ait = aintervals.find(n);
    if(ait==aintervals.end()){      
      aintervals.insert(std::make_pair(n,std::vector<TLoc >()));
    }
    ait = aintervals.find(n);
    assert(ait!=aintervals.end());
    int fmin,fmax;
    //assert(curroffsets[i]>=0);
#ifdef DEBUGGING
    std::cout << "Interval orient " << currorients[i] << " offset" << curroffsets[i] << " spans:" << currspanlens[i] <<  " len: " << currseqlens[i] << std::endl;
#endif
    if(currorients[i] == '+'){
      fmin=curroffsets[i];
      fmax=curroffsets[i]+currspanlens[i];
    }
    else{
      fmax=currseqlens[i]-curroffsets[i];
      fmin=currseqlens[i]-(curroffsets[i]+currspanlens[i]);
    }
#ifdef DEBUGGING
    std::cout << "Intervals for " << fmin << "-" << fmax <<  " " << currseqlens[i] << std::endl;
#endif
    assert(fmin>=0);
    assert(fmax<=(int)currseqlens[i]);
    //-1 - duplication end
    //0 align end
    //1 duplication start
    //2 align start
    if(dup){
      //ait->second.push_back(make_pair(fmin,1));//align start
      //ait->second.push_back(make_pair(fmax,-1));//align end 
      TLoc t1,t2;
      t1.first = fmin;
      t1.second = 1;
      t1.blocknum = lcbid;
      ait->second.push_back(t1);
      t2.first = fmax;
      t2.second = -1;
      t2.blocknum = lcbid;
      ait->second.push_back(t2);
    }
    else{
      //ait->second.push_back(make_pair(fmin,2));//align start
      //ait->second.push_back(make_pair(fmax,0));//align end
      TLoc t1,t2;
      t1.first = fmin;
      t1.second = 2;
      t1.blocknum = lcbid;
      ait->second.push_back(t1);
      t2.first = fmax;
      t2.second = 0;
      t2.blocknum = lcbid;
      ait->second.push_back(t2);
    }
  }
}

template<typename TGraph,
	 typename TLCBs,
	 typename TStringSet1,
	 typename TStringSet2,
	 typename TNames,
	 typename TGenomeNames,
	 typename TVertexOrientMap,
	 typename TStream1,
	 typename TName,
	 typename TLoc,
	 typename TScore>
std::vector<s_score> alignLCBs(TGraph &g,
			       TLCBs &LCBs,
			       TStringSet1 &seqSet,
			       TStringSet2 &genomeSeqSet,
			       TNames &sequenceNames,
			       TGenomeNames &genomeNames,
			       TVertexOrientMap &vertexOrientMap,
			       TStream1 &strmmaf,
			       std::map<TName,std::vector<TLoc> > &aintervals,
			       MsaOptions<Dna5 , TScore> const &msaOpt){

  typedef double TDistanceValue;
  typedef unsigned TSize;
  typedef Dna5 TAlphabet;
  typedef String<TDistanceValue> TDistanceMatrix;
  typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
  typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
  typedef typename Value<TScore>::Type TScoreValue;

  std::vector<s_score> allscores;
  TDistanceMatrix distanceMatrix;

  Graph<Tree<TDistanceValue> > seqguideTree;
  std::map<std::string, Graph<Tree<TDistanceValue> > > seqguideTrees;
  
  TSize nSeq = length(seqSet);
  
  TSize nGenomes=0;
  for(TSize i=0;i<length(genomeNames);i++){
    nGenomes = (genomeNames[i] > nGenomes) ? genomeNames[i] : nGenomes;
  }
  nGenomes = nGenomes+1;

#ifdef SEQAN_PROFILE
    std::cerr << "Saving interval tree marking location of duplications" << std::endl;
#endif   
  //Save interval tree of duplications
  //Save interval trees
  typedef IntervalAndCargo<int,TSize> TInterval;
  typedef Graph<Directed<void,WithoutEdgeId> > TIGraph;
  typedef IntervalTreeNode<TInterval> TNode;
  typedef String<TNode> TPropertyMap;
  String<String<TInterval> > dintervals;
  resize(dintervals,length(seqSet));
  String<TIGraph> dgs;
  String<TPropertyMap> dpms;
  for(int i=0;i<(int)length(seqSet);i++){
    std::map<int,std::pair<int,int> > tmpintervals;
    std::map<int,std::pair<int,int> >::iterator pos;
    bool inserted=false;
    typename std::map<TName,std::vector<TLoc> >::iterator ait=aintervals.find(sequenceNames[i]);
    if(ait!=aintervals.end()){
      for(typename std::vector<TLoc>::iterator pit = ait->second.begin();pit!=ait->second.end();pit++){
	boost::tie(pos, inserted) = tmpintervals.insert(std::make_pair(pit->blocknum,std::make_pair(0,0)));
	if(pit->second==1){
	  pos->second.first = pit->first;
	}
	else{
	  if(pit->second==-1){
	    pos->second.second = pit->first;
	  }
	}
      }
    }
    for(std::map<int,std::pair<int,int> >::iterator it = tmpintervals.begin();it!=tmpintervals.end();++it){
      //std::cout << it->second.first << " " << it->second.second << " bnum:" << it->first << std::endl;
      appendValue(dintervals[i],IntervalAndCargo<int,unsigned int>(it->second.first,it->second.second,it->first)); 
    }
  }
  
  resize(dgs,nSeq);
  resize(dpms,nSeq);
  for(unsigned int i=0;i<nSeq;i++){
    unsigned center = length(seqSet[i])/2;
    createIntervalTree(dgs[i],dpms[i],dintervals[i],center);
    //intervals for sequence i are not needed anymore
    clear(dintervals[i]);
  }
#ifdef SEQAN_PROFILE
    std::cerr << "Saving interval tree done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;
#endif
#ifdef SEQAN_PROFILE
    std::cerr << "Sorting vertices on each seq" << std::endl;
#endif   
  //Sort vertices on seq 
  std::vector<std::vector<TVertexDescriptor> > vseqs;
  vseqs.resize(nSeq);
  typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
  TVertexIterator it(g);
  for(;!atEnd(it);goNext(it)) {
    TVertexDescriptor v = *it;
    if(degree(g,v)>0){
      assert(sequenceId(g,*it)<vseqs.size());
      vseqs[sequenceId(g,v)].push_back(v);
    }
  }
  for(unsigned int i=0;i<nSeq;i++){
    //Sort vertices on seq 
    sort(vseqs[i].begin(),vseqs[i].end(),vertexposcmp<TGraph>(g));
  }
#ifdef SEQAN_PROFILE
    std::cerr << "Sorting vertices done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;
#endif
  unsigned int lcbid=0;

  //
  //Loop over each LCB and align
#ifdef SEQAN_PROFILE
  std::cerr << "Aligning " << LCBs.size() << " LCBs" << std::endl;
#endif

  //For tracking if an anchor has been aligned in an LCB
  std::set<TVertexDescriptor> coveredSet;

  // TODO, for parallel mugsy, refactor into a sub-process and parallel loop
  for(typename std::vector<std::vector<TVertexDescriptor> >::iterator lit = LCBs.begin();lit!=LCBs.end();lit++){
#ifdef DEBUGGING
    std::cout << "LCB:" << lcbid << " num_anchors:" << lit->size() << std::endl;	  
#endif
#ifdef SEQAN_PROFILE2
    std::cout << "LCB of size " << lit->size() << std::endl;	  
    std::cout << "LCB Initializing segments for LCB" << std::endl;
#endif	
#ifdef SEQAN_PROFILE
    if(lcbid%1000==0){
      std::cout << ".";
    }
#endif
    
    //Matches, scores, seqs, ids for current LCB
    typedef String<Fragment<> > TFragmentString;
    typedef typename Value<TScore>::Type TScoreValue;
    typedef String<TScoreValue> TScoreValues;
    typedef String<TAlphabet> TSequence;
    TFragmentString currmatches;
    TScoreValues currscores;
    StringSet<TSequence, Owner<> > currseqs;
    TNames currnameSet;
    
    //For tracking substrings 
    //std::map<TSize,unsigned int> offsets;
    //std::map<TSize,unsigned int> spanlens;
    //std::map<TSize,unsigned int> seqlens;
    //std::map<TSize,char> orients;
    std::map<TSize,s_offset> offsets;
    
    //Copy links between set of vertices in LCB $lit
    //from Graph $g and store in $currmatches,$currscores,$currseqs

    std::vector<TVertexDescriptor> currlcb;
    std::set<TVertexDescriptor> currlcbset;

    //Unless allownestedlcbs, each anchor vertex can contribute to
    //exactly one LCB; the longest LCB spanning the anchor
    //If the vertex is reported in subsequent LCBs it will be skipped
    if(msaOpt.allownestedlcbs == "true"){
      //skip checks if anchor has already been aligned
      //default for allownestedlcbs is false
    }
    else{
      for(std::vector<unsigned int>::const_iterator vit = lit->begin();vit!=lit->end();++vit){
	if(coveredSet.find(*vit)==coveredSet.end()){
	  currlcb.push_back(*vit);
	}
      }
      *lit=currlcb;
    }
    
    retrieveLCBSegments(g,
			seqSet,
			vertexOrientMap,
			lit,
			lcbid,
			sequenceNames,
			currseqs,
			currmatches,
			currscores,
			currnameSet,
			offsets,
			coveredSet,
			vseqs,
			boost::lexical_cast<int>(msaOpt.minlength));

    
    //TODO, Add matches for duplications and overlapping matches



#ifdef DEBUGGING
    std::cout << "LCB: " << lcbid 
	      << " vertices:" << lit->size()
	      << " seqset:" << length(seqSet)
	      << " offsets:" << offsets.size() 
	      << " currseqs:" << length(currseqs)
	      << " sequenceNames:" << length(sequenceNames)
	      << " currnameset:" << length(currnameSet)
	      << std::endl;
#endif
    assert(length(sequenceNames)==length(seqSet));
    assert(length(currnameSet)==length(currseqs));

    if(length(currseqs)>1 && length(currmatches)>0){
      assert(length(currmatches)>0);
#ifdef SEQAN_PROFILE2
      std::cout << "LCB Building alignment graph" << std::endl;
#endif
      //Build new graph from matches
      
      //Since LCBs contain no reversals
      //All matches should be relative to the forward/leading strand only
      //(non-reversed here)
      TGraph currG(currseqs);
      buildAlignmentGraph(currmatches, currscores, currG, FractionalScore());
#ifdef DEBUGGING
      std::cout << "Graph built V:"  << numVertices(currG) << " E:" << numEdges(currG) 
		<< " number of seqs:" << length(currseqs)
		<< " number of matches:" << length(currmatches) 
		<< " number of scores:" << length(currscores)
		<< std::endl;
#endif

      //Double check edge weights
      typedef typename Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
      TEdgeIterator itE(currG);
      //Undo Hack that stores reverse complement matches using
      //negative edge weights
      for(;!atEnd(itE);goNext(itE)){
	if(cargo(value(itE))<0){
	  cargo(value(itE)) = cargo(value(itE))*-1;
	}
      }

      //TESTING Test code to detect additional matches between
      //disconnected vertices in the LCB
      bool hashnonmatches=false;
      if(hashnonmatches){
	std::map<std::string,std::vector<TVertexDescriptor> > vhash;
	typename std::map<std::string,std::vector<TVertexDescriptor> >::iterator vhashpos;
	bool inserted;
	//Attempt to combine non-matching vertices
	typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
	std::cout << "Attempting to add addl edges" << std::endl;
	int newedges=0;
	TVertexIterator itV(currG);
	for(;!atEnd(itV);goNext(itV)){
	  if(degree(currG,*itV)==0){
	    String<char> sseq = infix(currseqs[sequenceId(currG,*itV)],fragmentBegin(currG,*itV),fragmentBegin(currG,*itV)+fragmentLength(currG,*itV));
	    char * c = toCString(sseq);
	    std::string cstr(c);
	    std::cout << sequenceId(currG,*itV) << " " << cstr << std::endl;
	    boost::tie(vhashpos, inserted) = vhash.insert(std::make_pair(cstr,std::vector<TVertexDescriptor>()));
	    vhashpos->second.push_back(*itV);
	  }
	}
	for(typename std::map<std::string,std::vector<TVertexDescriptor> >::iterator hit=vhash.begin();hit!=vhash.end();++hit){
	  for(typename std::vector<TVertexDescriptor>::iterator vit1=hit->second.begin();vit1!=hit->second.end();++vit1){
	    for(typename std::vector<TVertexDescriptor>::iterator vit2=vit1+1;vit2!=hit->second.end();++vit2){
	      assert(vit1!=vit2);
	      addEdge(currG,*vit1,*vit2,1);
	      newedges++;
	    }
	  }
	}
	std::cout << "Added " << newedges << " new edges" << std::endl;
      }

#ifdef DEBUGGING
      TEdgeIterator itE2(currG);
      for(;!atEnd(itE2);goNext(itE2)){
	assert(cargo(value(itE2))>0);
      }
#endif
      
      //Map between LCB array (curr*) and segment graph array
      std::vector<unsigned int> curroffsets;
      std::vector<unsigned int> currspanlens;
      std::vector<unsigned int> currseqlens;
      std::vector<char> currorients;
      assert(length(currseqs)==length(currnameSet));
      currorients.resize(length(currnameSet));
      currseqlens.resize(length(currnameSet));
      currspanlens.resize(length(currnameSet));
      curroffsets.resize(length(currnameSet));
      String<int> relevant_segments;
      std::set<unsigned int> currgenomes; //List of sequence ids
#ifdef DEBUGGING_GRAPH
      std::fstream strminfo;
      std::string lcbgname = "lcbgraph"+boost::lexical_cast<std::string>(lcbid)+".info";
      strminfo.open(lcbgname.c_str(), std::ios_base::out | std::ios_base::trunc);
#endif
      //TODO, refactor using a id map
      for(TSize currrow = 0; currrow<length(currnameSet); ++currrow) {
	for(TSize row = 0; row<length(sequenceNames); ++row) {
	  if(currnameSet[currrow]==sequenceNames[row]){
#ifdef DEBUGGING_GRAPH
	    strminfo << "SeqId:" << currrow << " " << currnameSet[currrow] 
		     << " origId:" << row 
		     << " offset:" << offsets[row].offset 
		     << " orient:" << offsets[row].orient << std::endl;
#endif
	    curroffsets[currrow] = offsets[row].offset;
	    currspanlens[currrow] = offsets[row].spanlen;
	    currseqlens[currrow] = offsets[row].seqlen;
	    currorients[currrow] = offsets[row].orient;
	    //See if we overlap a duplication
	    findIntervals(dgs[row],dpms[row],offsets[row].offset,relevant_segments);
	    findIntervals(dgs[row],dpms[row],offsets[row].offset+offsets[row].spanlen,relevant_segments);
	    currgenomes.insert(genomeNames[row]);
#ifdef DEBUGGGING
	    std::cout << currnameSet[currrow] << " " << currorients[currrow] << std::endl;
#endif
	  }
	}
      }
#ifdef DEBUGGING_GRAPH
      strminfo.close();
#endif
      
#ifdef SEQAN_PROFILE2
      std::cout << "LCB Building alignment graph done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;
      std::cout << "LCB Scoring matches" << std::endl;
#endif
      //Use or build guide tree based on genome sequences
      typedef double TDistanceValue;
      Graph<Tree<TDistanceValue> > currguideTree;
      assert(currgenomes.size()>0);
      getGuideTree(genomeSeqSet,currgenomes,seqguideTrees,currguideTree);

      TGraph currgOut(currseqs);
      s_score sscores = alignSingleLCB(currG,
				       currgOut,
				       ++lcbid,
				       currseqs,
				       currguideTree,
				       msaOpt);
#ifdef SCORING
      //Map to orig seq ids
      //TODO, refactor using a id map
      String<TSize> colCount;
      resize(colCount,nSeq+1);
      for(TSize i=0;i<nSeq;++i){
	for(TSize j=0;j<length(currseqs);++j){
	  if(sequenceNames[i]==currnameSet[j]){
	    colCount[i+1] += sscores.colCount[j+1];
	  }
	}
      }
      sscores.colCount=colCount;
      allscores.push_back(sscores);
#endif

      if(msaOpt.unique == "true"){
	//save interval 
	saveInterval(aintervals,
		     currnameSet,
		     curroffsets,
		     currspanlens,
		     currseqlens,
		     currorients,
		     lcbid);
      }
      if(strmmaf.is_open()){
	//Optionally write output
	//mafformat defined in refinement/graph_impl_align.h
	std::ostringstream lcblabel;
	lcblabel << " label=" << lcbid;

	if(msaOpt.duplications == "true"){
	  std::set<int> dblocks;
	  for(unsigned int i=0;i<length(relevant_segments);i++){
	    dblocks.insert(relevant_segments[i]);
	  }
	  for(std::set<int>::iterator dit=dblocks.begin();dit!=dblocks.end();++dit){
	    //std::cout << "LCB:" <<lcbid << " overlaps duplicated block d" << *dit << std::endl;
	    lcblabel << " dup=d" << *dit;
	  }
	}
	write(strmmaf,currgOut,currnameSet,MafFormat(),curroffsets,currspanlens,currseqlens,currorients,lcblabel.str());
	strmmaf.flush();
      }
    }
  }
#ifdef SEQAN_PROFILE
    std::cerr << "Aligning LCBs done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;
#endif
  return allscores;
}

template<typename TGraph, 
	 typename TStringSet,
	 typename TStringSet2,
	 typename TNames, 
	 typename TGenomeNames, 
	 typename TScore, 
	 typename TLCBs,
	 typename TVMap,
	 typename TStream,
	 typename TIMap>
void wholeGenomeAlignment(TGraph &g,
			  TStringSet &seqSet,
			  TStringSet2 &genomeSeqSet,
			  TNames &sequenceNames,
			  TGenomeNames &genomeNames,
			  MsaOptions<Dna5, TScore> const& msaOpt,
			  TLCBs &LCBs,
			  TVMap &vertexOrientMap,
			  TStream &strmmaf,
			  TIMap &aintervals){
  typedef Dna5 TAlphabet;
  typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
  typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
  //Perform consistency scoring
  //Do not add any edges to the graph
  //Simple score existing match edges for consistency
  
  //extension preserves directionality of matches
  //Edge weight > 0 same strand
  //Edge weight < 0 opposite strand
#ifdef SEQAN_PROFILE
  std::cerr << "Performing consistency scoring for connected edges only" << std::endl;
#endif
#ifdef NDEBUG
  ;
#else
  unsigned int nEdges = numEdges(g);
#endif
  //tripletLibraryExtensionCond(g,false);
  //tripletLibraryExtension(g, genomeguideTree, threshold / 2);
  graphBasedTripletLibraryExtension(g,false);
  std::cerr << "Num edges after consistency scoring " << numEdges(g) << std::endl;

  assert(nEdges==numEdges(g));
#ifdef SEQAN_PROFILE
  std::cerr << "Consistency scoring done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;
#endif

#ifdef DEBUGGING
  double vm, rss;
  process_mem_usage(vm, rss);
  cout << "VM: " << vm << "; RSS: " << rss << endl;
#endif
  
  /*
  if(msaOpt.filter){
#ifdef SEQAN_PROFILE
    std::cerr << "Filtering segment graph" << std::endl;
#endif
    filterSegmentGraph(g,seqSet,genomeNames,genomeguideTree);
#ifdef SEQAN_PROFILE
    std::cerr << "Filtering segment graph done: " 
	      << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;
#endif	
  }
  */
  //*******
  //Generate LCBs
  //
  //std::vector<std::vector<TVertexDescriptor> > LCBs;
  //std::map<TVertexDescriptor,char> vertexOrientMap;
  #ifdef TIMING 
  time(&now);
  std::cerr << "TIME ALIGNMENT_GRAPH:" << lasttime << " " << now << " " << now-lasttime << std::endl;
  lasttime=now;
  #endif 
  generateLCBs(g,
	       LCBs,
	       seqSet,
	       sequenceNames,
	       genomeNames,
	       vertexOrientMap,
	       aintervals,
	       msaOpt);


#ifdef SEQAN_PROFILE
  std::cerr << "Generating " << LCBs.size() << " alignments " << std::endl;
#endif
#ifdef DEBUGGING
  process_mem_usage(vm, rss);
  cout << "VM: " << vm << "; RSS: " << rss << endl;
#endif
  typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
  typedef typename EdgeType<TGraph>::Type TEdgeStump;
  typename std::vector<std::vector<unsigned int> >::const_iterator lit;
  typedef Fragment<> TFragment;
  typedef String<TAlphabet> TSequence;
  //Retrieve LCBs from the complete alignment graph $g
  std::vector<s_score> allscores = alignLCBs(g,
					     LCBs,
					     seqSet,
					     genomeSeqSet,
					     sequenceNames,
					     genomeNames,
					     vertexOrientMap,
					     strmmaf,
					     aintervals,
					     msaOpt);
  #ifdef TIMING 
  time(&now);
  std::cerr << "TIME ALIGN_LCB:" << lasttime << " " << now << " " << now-lasttime << std::endl;
  lasttime=now;
  #endif 

#ifdef SEQAN_PROFILE
  std::cerr << std::endl 
	    << "Finished aligning LCBs: " 
	    << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;
#endif 
  
#ifdef SCORING

#ifdef SEQAN_PROFILE
  std::cerr << "Calculating scores "<< std::endl;
#endif

  typedef unsigned int TSize;
  typedef typename Value<TScore>::Type TScoreValue;
  TSize numGapEx = 0;
  TSize numGap = 0;
  TSize numPairs = 0;
  TSize numIdents = 0;
  TSize alignLen = 0;
  TSize totalLen = 0;
  TSize lcbCount = 0;
  unsigned int minLen = std::numeric_limits<unsigned int>::max();
  unsigned int maxLen = 0;
  String<TSize> colCount;
  String<TSize> seqCount;
  TSize alignScore=0;
  TSize nSeq = length(seqSet);
  TSize nGen = length(genomeNames);
  fill(colCount,nSeq+1,0);
  fill(seqCount,nSeq+1,0);
  for(TSize i=0;i<nSeq;i++){
    assert(colCount[i]==0);
    assert(seqCount[i]==0);
  }
  for(std::vector<s_score>::iterator sit=allscores.begin();sit!=allscores.end();++sit){
    TSize nSeqn = sit->seqCount;
    seqCount[nSeqn]++;
    lcbCount++;
    alignScore += sit->alignScore;
    numGap += sit->numGap;
    numGapEx += sit->numGapEx;
    numPairs += sit->numPairs;
    numIdents += sit->numIdents;
    minLen = (sit->alignLen < minLen) ? sit->alignLen : minLen;
    maxLen = (sit->alignLen > maxLen) ? sit->alignLen : maxLen;
    alignLen += sit->alignLen;
    totalLen += sit->totalLen;
    for(TSize i=0;i<nSeqn;++i){
      colCount[i+1] += sit->colCount[i+1];
    }
  }
  
  std::string outfile(msaOpt.outfile);
  std::fstream strmscore;
  strmscore.open(std::string(outfile+".scores").c_str(), std::ios_base::out | std::ios_base::trunc);
  TScoreValue gop = msaOpt.sc.data_gap_open;
  TScoreValue gex = msaOpt.sc.data_gap_extend;
  strmscore << "Num LCBs: " << lcbCount << std::endl;
  strmscore << "Avg length: " << (float)alignLen/(float)lcbCount << "bp Range:" << minLen << "-" << maxLen << "bp " << std::endl;
  strmscore << "Total scoring summary over all LCBs" << std::endl;
  strmscore << "SP alignment Score: " << alignScore << std::endl;
  strmscore << "Alignment Length: " << alignLen << std::endl;
  strmscore << "Sum of sequence length: " << totalLen << std::endl;
  strmscore << "#Match-Mismatch pairs: " << numPairs << std::endl;
  strmscore << "#Match pairs: " << numIdents << std::endl;
  strmscore << "Score contribution by match-mismatch pairs: " << (alignScore - (((TScoreValue) numGap * gop) + ((TScoreValue) numGapEx * gex))) << std::endl;
  strmscore << "#Gap extensions: " << numGapEx << std::endl;
  strmscore << "Score contribution by gap extensions: " << ((TScoreValue) numGapEx * gex) << std::endl;
  strmscore << "#Gap openings: " << numGap << std::endl;
  strmscore << "Score contribution by gap openings: " << ((TScoreValue) numGap * gop) << std::endl;
  //strmscore << "Average percent identity: " << << std::endl;
  //strmscore << "Average percent aligned: " << << std::endl;
  if(nSeq!=nGen){
    strmscore << "WARNING: Some of the following stats are not calculated correctly for incomplete genomes" << std::endl;
  }
  strmscore << "Count of columns with identical characters " << std::endl;
  for(TSize i=0;i<=nSeq;++i){
    if(i!=0)
      strmscore << " " << i << ":" << colCount[i];
  }
  strmscore << std::endl;
  strmscore << "Counts of seqs per LCB " << std::endl;
  for(TSize i=1;i<=nSeq;++i){
    strmscore << " " << i << ":" << seqCount[i];
  }
  strmscore << std::endl;
  strmscore << "Lengths per seq" << std::endl;
  int minseq = std::numeric_limits<unsigned int>::max();
  int maxseq = 0;
  int totalunaligned=0;
  for(TSize i=0;i<nSeq;++i){
    minseq = (minseq < length(seqSet[i])) ? minseq : length(seqSet[i]);
    maxseq = (maxseq > length(seqSet[i])) ? maxseq : length(seqSet[i]);
    strmscore << sequenceNames[i] << " len: " << length(seqSet[i])  << " unaligned lcb,bp,%: " << " aligned lcb,bp,%: " << endl;
    //totalunaligned+=
  }
  //TODO these stats only work for completed genomes
#ifdef SCORING_NEW
  if(length(sequenceNames)==length(genomeNames)){
    //bionomial coffecient
    int nfac=1;
    int n1fac=1;
    for(int i=1;i<=nSeq;++i){nfac *= i;}
    for(int i=1;i<=(nSeq-2);++i){n1fac *= i;}
    assert(n1fac>0);
    int possiblematchpairsmin = (nfac/(2*n1fac))*minseq;
    int possiblematchpairsaln = (nfac/(2*n1fac))*alignLen;
    std::cout << nfac << " " << n1fac << std::endl;
    assert(possiblematchpairsmin>0);
    assert(possiblematchpairsaln>0);
    strmscore << "Estimate of average %id (using min seq len) " << (float)numIdents / (float)possiblematchpairsmin << std::endl;
    strmscore << "Estimate of average %id (using aln len) " << (float)numIdents / (float)possiblematchpairsaln << std::endl;
    strmscore << "Estimate of average %aln (using min seq len) " << (float)numPairs / (float)possiblematchpairsmin << std::endl;
    strmscore << "Estimate of average %aln (using aln len) " << (float)numPairs / (float)possiblematchpairsaln << std::endl;
    
    strmscore << "Estimate of overall %id (maxseq/alignLen+unaligned) " << (float)maxseq / (float)(alignLen+totalunaligned) << std::endl;
    strmscore << std::endl;
  }
#endif
  strmscore.close();
#ifdef SEQAN_PROFILE
  std::cout << std::endl 
	    << "Finished calculating scores " 
	    << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;
#endif 
  
#endif 

}  

//-1 - duplication end
//0 align end
//1 duplication start
//2 align start

template<typename TStringSet,
	 typename TNames,
	 typename TName,
	 typename TLoc,
	 typename TStream>
void printUniques(TStringSet &seqSet,
		  TNames &sequenceNames,
		  std::map<TName,std::vector<TLoc> >&aintervals,
		  TStream &strmmaf){ 
  int icount=0;
  for(int i=0;i<(int)length(seqSet);i++){
    typename std::map<TName,std::vector<TLoc> >::iterator ait=aintervals.find(sequenceNames[i]);
    if(ait!=aintervals.end()){
      sort(ait->second.begin(),ait->second.end(),poscmp<TLoc>());
      int last=-1;
      int open=0;
      int indup=0;
      std::vector<int> currdups;
      //pair<coord,start_end>
      strmmaf << std::endl;
      assert(ait->second[ait->second.size()-1].first<=(int)length(seqSet[i]));
      for(typename std::vector<TLoc>::iterator pit = ait->second.begin();pit!=ait->second.end();pit++){
	//Start of alignment
	assert(pit->first>=0);
	assert((unsigned)pit->first<=length(seqSet[i]));
	if(pit->second>0){
	  if(pit->second==1){//duplication start
	    indup++;
	    currdups.push_back(pit->blocknum);
	  }
	  else{
#ifdef DEBUGGING
	    std::cout << "OPEN:" << open << " type:" << pit->second << " seq: " << ait->first << " coord:" << pit->first << " last_close:" << last 
		      << " spanlen: " << length(infix(seqSet[i],last,pit->first)) 
		      << " == " <<  pit->first - last 
		      << " indup:" << indup <<std::endl;
#endif
	    if(open==0){
	      if(last>0 && pit->first-last>0){
		if(indup){
		  //print as dup
		  //std::cout << "DUP" << std::endl;
		  strmmaf << "a score=0 label=u" << icount++ << " mult=1 dup=";
		  for(std::vector<int>::iterator it =currdups.begin();it!=currdups.end();++it){
		    strmmaf << "d" << *it;
		    if(it+1!=currdups.end()){
		      strmmaf << ",";
		    }
		  }
		  strmmaf << std::endl;
		}
		else{
		  strmmaf << "a score=0 label=u" << icount++ << " mult=1" << std::endl;
		}
		assert((int)length(infix(seqSet[i],last,pit->first)) == pit->first - last);
		strmmaf << "s " << ait->first << " " << last << " " << pit->first - last << " + " 
			<< length(seqSet[i]) << " " << infix(seqSet[i],last,pit->first)
			<< std::endl
			<< std::endl;
	      }
	    }
	    if(pit->second==2){
	      open++;
	    }
	  }
	}
	else{
	  if(pit->second==-1){//duplication stop
	    indup--;
	    currdups.pop_back();
	  }
	  if(pit->second==0){//End of alignment
	    open--;
	  }
	  last=pit->first;
#ifdef DEBUGGING
	  std:: cout << "CLOSE: " << open << " type:" << pit->second << " coord:" << last << std::endl;
#endif
	}
	if(last<0){
	  last=0;
	}
	assert(open>=0);
	assert(last<=(int)length(seqSet[i]));
      }
      //assert(open==0);
      if(last && length(seqSet[i])-last>0 && open==0){
#ifdef DEBUGGING
	std::cout << "Printing to end of sequence " << last << "-" << length(seqSet[i]) << std::endl;
#endif
	if(indup){
	  //print as dup
	  //std::cout << "DUP" << std::endl;
	  strmmaf << "a score=0 label=u" << icount++ << " mult=1 dup=";
	  for(std::vector<int>::iterator it =currdups.begin();it!=currdups.end();++it){
	    strmmaf << "d" << *it;
	    if(it+1!=currdups.end()){
	      strmmaf << ",";
	    }
	  }
	  strmmaf << std::endl;
	}
	else{
	  strmmaf << "a score=0 label=u" << icount++ << " mult=1" << std::endl;
	}
	assert(length(infix(seqSet[i],last,length(seqSet[i]))) == length(seqSet[i]) - last);
	strmmaf << "s " << ait->first << " " << last << " " << length(seqSet[i]) - last << " + " 
		<< length(seqSet[i]) << " " << infix(seqSet[i],last,length(seqSet[i])) 
		<< std::endl
		<< std::endl;
      }
    }
    else{
#ifdef DEBUGGING
      std::cout << "No alignment on sequence " << sequenceNames[i] << std::endl;
      std::cout << "Printing entire sequence" << std::endl;
#endif
      strmmaf << "a score=0 label=0 mult=1" << std::endl;
      strmmaf << "s " << sequenceNames[i] << " 0 " << length(seqSet[i]) << " + " 
	      << length(seqSet[i]) << " " << seqSet[i]
	      << std::endl
	      << std::endl;
    }
  }
}

template<typename TStringSet, 
    typename TCargo, 
    typename TSpec, 
    typename TStringSet1, 
    typename TNames, 
    typename TGenomeNames,
    typename TIntervals,
    typename TScore>
void
singlepass_wholeGenomeAlignment(Graph<Alignment<TStringSet, TCargo, TSpec> >& gAlign, 
				TStringSet1& sequenceSet,
				TNames& sequenceNames,
				TGenomeNames& genomeNames,
				TIntervals &aintervals,
				MsaOptions<Dna5 , TScore> const& msaOpt)
{
  typedef Dna5 TAlphabet;
  typedef typename Value<TScore>::Type TScoreValue;
  typedef typename Size<TStringSet>::Type TSize;
  typedef typename Value<TStringSet1>::Type TString;
  typedef typename Value<TNames>::Type TName;
  //typedef Graph<Alignment<TStringSet, TSize> > TGraph;
  //Using int to support negative edge scores
  typedef Graph<Alignment<TStringSet, int> > TGraph;
  typedef typename Id<TGraph>::Type TId; 
  typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
  typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
  
  typedef double TDistanceValue;

#ifdef SEQAN_PROFILE
  std::cerr << "Mugsy WGA" << std::endl;
  std::cerr << "Reading sequences and alignments " << std::endl;
  std::cerr << "--distance=" << msaOpt.distance << std::endl;
  std::cerr << "--minlength=" << msaOpt.minlength << std::endl;
#endif
  // Initialize alignment object
  clear(gAlign);
  assignStringSet(gAlign, sequenceSet);
  // Some alignment constants
  TStringSet& seqSet = stringSet(gAlign);
  TSize nSeq = length(seqSet);
  TSize nGenomes=0;
  for(TSize i=0;i<length(genomeNames);i++){
    nGenomes = (genomeNames[i] > nGenomes) ? genomeNames[i] : nGenomes;
  }
  nGenomes = nGenomes+1;
  std::cerr << "Number of genomes:" << nGenomes << std::endl;
  std::cerr << "Number of sequences:" << nSeq << std::endl;
      
  // Containers for segment matches and corresponding scores 
  typedef String<Fragment<> > TFragmentString;
  TFragmentString matches;
  typedef String<TScoreValue> TScoreValues;
  TScoreValues scores;
  
  // Include segment matches from subalignments
  if (!empty(msaOpt.alnfiles)) {
    typedef typename Iterator<String<std::string>, Standard>::Type TIter;
    TIter begIt = begin(msaOpt.alnfiles, Standard() );
    //TIter begItEnd = end(msaOpt.alnfiles, Standard() );
    //Only read first alignment file
    //for(;begIt != begItEnd; goNext(begIt)) {
#ifdef SEQAN_PROFILE
    std::cerr << "*Alignment file XMFA format: " << (*begIt).c_str() << std::endl;
#endif
    std::ifstream strm_lib;
    strm_lib.open((*begIt).c_str(), ::std::ios_base::in | ::std::ios_base::binary);
    read(strm_lib, matches, scores, sequenceSet, sequenceNames, MultiFastaAlign());
    strm_lib.close();
    //  }
  }
  /*
  //TODO, read mummer for defining MUMi
  // Include MUMmer segment matches
  if (!empty(msaOpt.mummerfiles)){
#ifdef SEQAN_PROFILE
    std::cout << "Parsing MUMmer segment matches:" << std::endl;
#endif
    String<char> mummerFiles = value(msaOpt.mummerfiles);
    String<char> currentMumFile;
    for(TSize i = 0; i<=length(mummerFiles); ++i) {
      if ((i == length(mummerFiles) || (value(mummerFiles, i) == ','))) {		
#ifdef SEQAN_PROFILE
	std::cout << "*MUMmer file: " << currentMumFile << std::endl;
#endif
	std::stringstream input;
	input << currentMumFile;
	std::ifstream strm_lib;
	strm_lib.open(input.str().c_str(), std::ios_base::in | std::ios_base::binary);
	read(strm_lib, matches, scores, seqSet, sequenceNames, MummerLib());		
	strm_lib.close();
	clear(currentMumFile);
      } else {
	if ((value(mummerFiles, i) != ' ') && (value(mummerFiles, i) != '\t')) appendValue(currentMumFile, value(mummerFiles, i));
      }
    }
#ifdef SEQAN_PROFILE
    std::cout << "Parsing done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;
#endif
  }
  */
#ifdef SEQAN_PROFILE
  std::cerr << "Reading alignments done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;
#endif
#ifdef DEBUGGING
  double vm, rss;
  process_mem_usage(vm, rss);
  cout << "VM: " << vm << "; RSS: " << rss << endl;
#endif

  //Build StringSet for each genome
  TStringSet1 genomeSeqSet;
  TSize seqSetLen = length(seqSet);
  std::vector<std::vector<TSize> > genomeMap;
  std::map<TSize,TSize> genomeLenMap;
  genomeMap.resize(nGenomes);
  for(TSize i=0;i<seqSetLen;++i){
    genomeMap[genomeNames[i]].push_back(i);
  }
  TSize mapLen=genomeMap.size();
  resize(genomeSeqSet,mapLen,Exact());
  for(TSize i=0;i<mapLen;++i){
    std::stringstream ss;
    //Concatenate all sequences for the genome
    for(typename std::vector<TSize>::iterator sit = genomeMap[i].begin();sit != genomeMap[i].end();++sit){
      append(genomeSeqSet[i],seqSet[*sit]);
    }
  }
  //Save sum of genome lengths for later
  for(unsigned int i=0;i<length(genomeSeqSet);++i){
    genomeLenMap[i]=length(genomeSeqSet[i]);
  }
#ifdef SEQAN_PROFILE
  std::cerr << "Building guide trees" << std::endl;
#endif

  /*
  // Set-up a distance matrix
  typedef String<TDistanceValue> TDistanceMatrix;
  TDistanceMatrix distanceMatrix;
 
  clear(distanceMatrix);
  //Calculate initial
  //Guide tree over all genomes
  typedef Graph<Tree<TDistanceValue> > TGuideTree;
  TGuideTree genomeguideTree;
  TSize ktup=3; //3mers
  getKmerSimilarityMatrix(genomeSeqSet, distanceMatrix, ktup, TAlphabet());
  // Similarity to distance conversion
  typedef typename Value<TDistanceMatrix>::Type TValue;
  typedef typename Iterator<TDistanceMatrix, Standard>::Type TMatrixIterator;
  TMatrixIterator matIt = begin(distanceMatrix, Standard());
  TMatrixIterator endMatIt = end(distanceMatrix, Standard());
  for(;matIt != endMatIt;++matIt) 
    *matIt = SEQAN_DISTANCE_UNITY - (*matIt);
  if (msaOpt.build == 0) njTree(distanceMatrix, genomeguideTree);
  else if (msaOpt.build == 1) upgmaTree(distanceMatrix, genomeguideTree, UpgmaMin());
  else if (msaOpt.build == 2) upgmaTree(distanceMatrix, genomeguideTree, UpgmaMax());
  else if (msaOpt.build == 3) upgmaTree(distanceMatrix, genomeguideTree, UpgmaAvg());
  else if (msaOpt.build == 4) upgmaTree(distanceMatrix, genomeguideTree, UpgmaWeightAvg());
  clear(distanceMatrix);
  */
#ifdef SEQAN_PROFILE
  std::cerr << "Building guide trees done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;
#endif
  //*******
  //Build alignment graph
  //
  // Use these segment matches for the initial alignment graph
#ifdef SEQAN_PROFILE
  std::cerr << "Building alignment graph from " << length(matches) << " matches" << std::endl;
#endif
  TGraph g(seqSet);
  if (!msaOpt.rescore) buildAlignmentGraph(matches, scores, g, FractionalScore() );
  else buildAlignmentGraph(matches, scores, g, msaOpt.sc, ReScore() );
  //clear these here to save memory
  clear(matches);
  clear(scores); 
#ifdef SEQAN_PROFILE
  std::cerr << "Building alignment graph done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;
#endif
#ifdef DEBUGGING
  process_mem_usage(vm, rss);
  cout << "VM: " << vm << "; RSS: " << rss << endl;
#endif

  std::cerr << std::endl << "Refined alignment graph built. E: " << numEdges(g) << " V:" << numVertices(g) << std::endl;
  
  //Stats
  typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
  TVertexIterator itV(g);
  unsigned totalmatchingbp=0;
  unsigned totalseqlen=0;
  for(;!atEnd(itV);goNext(itV)){
    if(degree(g,*itV)>0){
      totalmatchingbp+=fragmentLength(g,*itV);
    }
  }
  for(unsigned int i=0;i<seqSetLen;i++){
    totalseqlen+=length(seqSet[i]);
  }
  std::cerr << "Average fragment length: " 
	    << (float)(totalmatchingbp/numVertices(g)) 
	    << "bp" << std::endl;
  std::cerr << "Percentage matching bp:" 
	    << totalmatchingbp << "/" << totalseqlen  
	    << "=" << (float)totalmatchingbp/totalseqlen 
	    << std::endl;
  //Calculate a distance measure similar to MUMi
  //Print range
  std::map<std::pair<int,int>,int > flengths;
  std::map<std::pair<int,int>,int>::iterator pos;
  typedef typename Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
  TEdgeIterator itE(g);
  TVertexDescriptor source,target;
  TEdgeDescriptor ed;
  bool inserted=false;
  TSize sgen,tgen;
  for(;!atEnd(itE);goNext(itE)){
    ed = *itE;
    source = getSource(ed);
    target = getTarget(ed);
    sgen = genomeNames[sequenceId(g,source)];
    tgen = genomeNames[sequenceId(g,target)];
    if(sgen>tgen){
      TSize tmp=sgen;
      sgen=tgen;
      tgen=tmp;
    }
    assert(fragmentLength(g,source)==fragmentLength(g,target));
    //TODO, fix may count same vertex more than once
    boost::tie(pos, inserted) = flengths.insert(std::make_pair(std::make_pair(sgen,tgen),fragmentLength(g,source)));
    if(!inserted){
      pos->second+=fragmentLength(g,source);
    }
  }
  float minnia=2;
  float maxnia=0;
  float minnim=2;
  float maxnim=0;
  for(std::map<std::pair<int,int>,int>::iterator mit = flengths.begin();mit!=flengths.end();++mit){
    float avgsize = (genomeLenMap[mit->first.first]+genomeLenMap[mit->first.second])/2;
    float minsize = genomeLenMap[mit->first.first] < genomeLenMap[mit->first.second] ? genomeLenMap[mit->first.first] : genomeLenMap[mit->first.second];
    assert(avgsize>0);
    float nia = 1 - (mit->second/avgsize);
    float nim = 1 - (mit->second/minsize);
    minnia = (minnia < nia) ? minnia : nia;
    maxnia = (maxnia > nia) ? maxnia : nia;
    minnim = (minnim < nim) ? minnim : nim;
    maxnim = (maxnim > nim) ? maxnim : nim;
  }
  flengths.clear();
  //clear(genomeSeqSet);
  std::cerr << "D=1-Lseq/Lavg min-max: " << minnia << "-" << maxnia << std::endl;
  std::cerr << "D=1-Lseq/Lmin min-max: " << minnim << "-" << maxnim << std::endl;
  //
#ifdef DEBUGGING_GRAPH
  std::fstream rawstrm;
  rawstrm.open("origrefinegraph.out", std::ios_base::out | std::ios_base::trunc);
  write(rawstrm,g,sequenceSet,Raw());
  rawstrm.close();
#endif
  //*******
  //
  //TODO 
  //testing partitioning
  //This partitioning is just for testing.
  //The refined graph is already partitioned.
  //Ideally matches could be filtered and partitioned prior to building the large alignment graph
  //Convert graph back to matches
  /*
  TFragmentString filtmatches;
  std::vector<TFragmentString> matchSets;
  if(doPartitioning){
    typedef typename Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
    TEdgeIterator itE(g);
    for(;!atEnd(itE);goNext(itE)){
      TEdgeDescriptor ed = *itE;
      TVertexDescriptor vd1 = getSource(ed);
      TVertexDescriptor vd2 = getTarget(ed);
      appendValue(filtmatches, Fragment<>(sequenceId(g,vd1), 
					  fragmentBegin(g,vd1), 
					  sequenceId(g,vd2),
					  fragmentBegin(g,vd2), 
					  fragmentLength(g,vd1), 
					  ((int)(cargo(ed)<0)) ? true : false));
    }
    partitionSegments(seqSet,filtmatches,matchSets,msaOpt.partition);
    for(typename std::vector<TFragmentString>::iterator it = matchSets.begin();
	it!=matchSets.end();it++){
      String<Fragment<> > matchset = *it;
#ifdef DEBUGGING
      std::cout << "Length of matchset " << length(matchset) << std::endl;
#endif
    }
  }
  else{
    //Using full graph
  }
  */
  std::fstream strmmaf;
  std::string outfile(msaOpt.outfile);
  strmmaf.open(std::string(outfile+".maf").c_str(), std::ios_base::out | std::ios_base::trunc);
  _streamWrite(strmmaf,"##maf version=1 scoring=mugsy");
  typedef String<unsigned int> TComponentMap;
  typedef typename Value<TComponentMap>::Type TComponent;
  typedef typename Position<TGraph>::Type TPos;
  typedef SVABlock<TComponent,unsigned,TVertexDescriptor,unsigned> TBlock;
  std::vector<std::vector<TVertexDescriptor> > lcbs;
  std::map<TVertexDescriptor,char> vertexOrientMap;
 
  //LCBs are saved in lcbsp
  //Optionally, can also store profiles and write directly to strmmaf
  wholeGenomeAlignment(g,
		       seqSet,
		       genomeSeqSet,
		       sequenceNames,
		       genomeNames,
		       msaOpt,
		       lcbs,
		       vertexOrientMap,
		       strmmaf,
		       aintervals);
  //Close out MAF
  strmmaf << std::endl;

  //loop over all profiles and print

  //Print all remaining unaligned sequences
  //TODO broken in refactor, fix 
  if(msaOpt.unique == "true"){
    printUniques(seqSet,sequenceNames,aintervals,strmmaf);
  }

  //Close output streams
  strmmaf.close();

#ifdef SEQAN_PROFILE
  std::cerr << "Mugsy all done" << std::endl;
#endif
}


//
//Input all pairwise matches in duplicated regions
//Build refinement graph to reduce into non-overlapping segments
//Sort over each sequence and build runs of regions < DUP_ADJ
//Save runs of length > DUP_CMB
template<typename TStringSet, 
	 typename TCargo, 
	 typename TSpec, 
	 typename TStringSet1, 
	 typename TNames, 
	 typename TGenomeNames,
	 typename TScore,
	 typename TIntervals>
inline void
findDuplications(Graph<Alignment<TStringSet, TCargo, TSpec> >& gAlign, 
		 TStringSet1& sequenceSet,
		 TNames& sequenceNames,
		 TGenomeNames& genomeNames,
		 TIntervals& dupintervals,
		 MsaOptions<Dna5 , TScore> const& msaOpt)
{
  typedef Dna5 TAlphabet;
  typedef typename Value<TScore>::Type TScoreValue;
  typedef typename Size<TStringSet>::Type TSize;
  typedef typename Value<TStringSet1>::Type TString;
  typedef typename Value<TNames>::Type TName;
  //Using int to support negative edge scores
  typedef Graph<Alignment<TStringSet, int> > TGraph;
  typedef typename Id<TGraph>::Type TId; 
  typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
  typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
  //
  typedef std::map<unsigned int, unsigned int> TComponentLength;
  
  // Strongly Connected Components, topological sort, and length of each component
  typedef String<unsigned int> TComponentMap;
  typedef typename Value<TComponentMap>::Type TComponent;
  typedef typename Position<TGraph>::Type TPos;
  typedef SVABlock<TComponent,TSize,TVertexDescriptor,TPos> TBlock;
  
  typedef typename Value<TComponentMap>::Type TComponent;
  typedef std::pair<TId, TSize> TKey;
  typedef std::map<TKey, TVertexDescriptor> TPosToVertexMap;
  typedef FragmentInfo<TId, TSize> TFragmentInfo;

  typedef double TDistanceValue;
#ifdef SEQAN_PROFILE
  std::cerr << "Detecting duplications " << std::endl;
#endif
  // Initialize alignment object
  clear(gAlign);
  assignStringSet(gAlign, sequenceSet);
  // Some alignment constants
  TStringSet& seqSet = stringSet(gAlign);
  TSize nSeq = length(seqSet);
  TSize nGenomes=0;
  for(TSize i=0;i<length(genomeNames);i++){
    nGenomes = (genomeNames[i] > nGenomes) ? genomeNames[i] : nGenomes;
  }
  nGenomes = nGenomes+1;
  std::cerr << "Number of genomes:" << nGenomes << std::endl;
  std::cerr << "Number of sequences:" << nSeq << std::endl;
  
  // Set-up a distance matrix
  typedef String<TDistanceValue> TDistanceMatrix;
  TDistanceMatrix distanceMatrix;
  
  // Containers for segment matches and corresponding scores 
  typedef String<Fragment<> > TFragmentString;
  TFragmentString matches;
  typedef String<TScoreValue> TScoreValues;
  TScoreValues scores;
  	
  // Include segment matches from subalignments
  if (!empty(msaOpt.alnfiles)) {
    typedef typename Iterator<String<std::string>, Standard>::Type TIter;
    TIter begIt = begin(msaOpt.alnfiles, Standard() );
    TIter begItEnd = end(msaOpt.alnfiles, Standard() );
    goNext(begIt);//alignment XMFA is second alignment file passed
    for(;begIt != begItEnd; goNext(begIt)) {
#ifdef SEQAN_PROFILE
      std::cerr << "*Alignment file XMFA format: " << (*begIt).c_str() << std::endl;
#endif
      std::ifstream strm_lib;
      strm_lib.open((*begIt).c_str(), ::std::ios_base::in | ::std::ios_base::binary);
      //defined in graph_align_tcoffee_io.h
      read(strm_lib, matches, scores, sequenceSet, sequenceNames, MultiFastaAlign());
      strm_lib.close();
      //clear(alignmentFile);read(strm_lib, matches, scores, sequenceNames, FastaAlign());
    }
  }
  else{
    assert(false);
  }
#ifdef SEQAN_PROFILE
    std::cerr << "Building alignment graph" << std::endl;
#endif
    TGraph g(seqSet);
    //defined in graph_align_tcoffee_base.h
    if (!msaOpt.rescore) buildAlignmentGraph(matches, scores, g, FractionalScore() );
    else buildAlignmentGraph(matches, scores, g, msaOpt.sc, ReScore() );
    //clear(matches);
    //clear(scores); 
#ifdef SEQAN_PROFILE
    std::cerr << "Building alignment graph done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;
#endif
    std::cerr << std::endl << "Refined alignment graph built. E: " << numEdges(g) << " V:" << numVertices(g) << std::endl;
    //Stats
    typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
    TVertexIterator itV(g);
    unsigned totalmatchingbp=0;
    unsigned totalseqlen=0;
    TSize seqSetLen = length(seqSet);
    for(;!atEnd(itV);goNext(itV)){
      if(degree(g,*itV)>0){
	totalmatchingbp+=fragmentLength(g,*itV);
      }
    }
    for(unsigned int i=0;i<seqSetLen;i++){
      totalseqlen+=length(seqSet[i]);
    }
    std::cerr << "Average fragment length: " 
	      << (float)(totalmatchingbp/numVertices(g)) 
	      << "bp" << std::endl;
    std::cerr << "Percentage matching bp:" 
	      << totalmatchingbp << "/" << totalseqlen  
	      << "=" << (float)totalmatchingbp/totalseqlen 
	      << std::endl;
    
#ifdef DEBUGGING_GRAPH
    std::fstream rawstrm;
    rawstrm.open("origrefinegraph.out", std::ios_base::out | std::ios_base::trunc);
    write(rawstrm,g,sequenceSet,Raw());
    rawstrm.close();
#endif
    
#ifdef DEBUGGING
  std::cout << "Finding connected components " << std::endl;
#endif

#ifdef SEQAN_PROFILE
    std::cerr << "Determining duplicated regions" << std::endl;
#endif
  
  //Sequence set that will capture each copy of the duplication
  TStringSet1 runSeqSet;
    //StringSet<TString, TSpec> runSeqSet;
  Graph<Directed<> > runG;

  // Connected Components
  // Each CC represents an UNGAPPED set of aligned fragments/segments across sequences
  // A CC is an ungapped block 
  // A CC is also an LCB at this point of the algorithm but may be extended
  TComponentMap componentall; 
  std::map<std::pair<TComponent,TComponent>,TBlock *> componentVertexMap;
  std::vector<std::vector<TBlock> > blocksbycomponent; 
  
  //TSize numComponents = connected_components(g, componentall);
  //TSize numComponents = connected_components_by_genome_ranked(g, componentall, genomeNames,std::numeric_limits<unsigned int>::max());
  TSize numComponents = connected_components_ranked(g, componentall);
  std::cerr << "Determined " << numComponents << " component segments in graph of size " << numVertices(g) << std::endl;
  assert(numComponents>0);
  //std::cerr << "Calculating positional scores" <<std::endl;
  //scorePosCons(g,componentall,numComponents,posScores,POS_ADJ);
  //std::cerr << "Set positional scores" << std::endl;

  //Identify runs
  int POS_CMB = 100;

  std::map<TSize,std::vector<TVertexDescriptor> > componentSeqMap; 
  std::map<TComponent,std::vector<TVertexDescriptor> > componentMap;
  
  typename TPosToVertexMap::const_iterator it2 = g.data_pvMap.begin();
  typename TPosToVertexMap::const_iterator it2End = g.data_pvMap.end();
  for(it2 = g.data_pvMap.begin();it2!=it2End;++it2) {
    TVertexDescriptor currV = it2->second;
    assert(getProperty(componentall,currV)==componentall[currV]);

    TSize currentSeq = sequenceId(g,currV);

    typename std::map<TSize,std::vector<TVertexDescriptor> >::iterator fit = componentSeqMap.find(currentSeq); 
    if(fit==componentSeqMap.end()){
      componentSeqMap[currentSeq] = std::vector<TVertexDescriptor>();
    }
    componentSeqMap[currentSeq].push_back(currV);
    TComponent c = getProperty(componentall, currV);
    if(componentMap.find(c)==componentMap.end()){
      componentMap[c] = std::vector<TVertexDescriptor>();
    }
    componentMap[c].push_back(currV);
  }
  

  std::map<TSize,TSize> seqIdxMap;
  std::map<TVertexDescriptor,TSize> runmap;
  std::map<TSize,std::vector<TVertexDescriptor> > vrunmap;
  int runcount=0;
  for(typename std::map<TSize,std::vector<TVertexDescriptor> >::iterator it = componentSeqMap.begin();it!=componentSeqMap.end();++it){
    std::set<std::pair<int,int> > runs;

    TSize currentSeq = it->first;
    //std::cout << "Examining sequence " << currentSeq
    //      << " with num vertices:" << it->second.size() << std::endl;
    //Sort vertices in G on sequence currentSeq
    std::set<TSize> repeatCC;
    std::vector<TVertexDescriptor> vlist;
    int lastcoord=0;
    int runstart=0;
    int runend=0;

    sort(it->second.begin(),it->second.end(),vertexposcmp<TGraph>(g));
    for(typename std::vector<TVertexDescriptor>::iterator vit = it->second.begin();vit!=it->second.end();++vit){
      TVertexDescriptor currV = *vit;
      TComponent c = getProperty(componentall, *vit);
      //Only consider segments that are part of matches
      if(componentMap[c].size()>1){
	if(lastcoord>0){
	  int dist = fragmentBegin(g,*vit)-lastcoord;
	  //it->first.first is CC label
	  if(dist>(int)POS_CMB 
	     || repeatCC.find(c)!=repeatCC.end()){
	    runend = lastcoord;
	    if(runend - runstart > POS_CMB){
	      //std::cout << runstart << " runstart " << runstart << " runend: " << runend << " len:" << runend - runstart << std::endl;
	      runs.insert(std::make_pair(runstart,runend));
	      for(unsigned i=0;i<vlist.size();++i){
		runmap[vlist[i]] = runcount;
		assert(sequenceId(g,vlist[i])==currentSeq);
	      }
	      if(vrunmap.find(runcount)==vrunmap.end()){
		vrunmap[runcount] = std::vector<TVertexDescriptor>();
	      }
	      vrunmap[runcount].insert(vrunmap[runcount].end(),vlist.begin(),vlist.end());
	      addVertex(runG);
	      seqIdxMap[runcount]=currentSeq;
	      runcount++;
	    }
	    repeatCC.clear();
	    vlist.clear();
	    runstart = fragmentBegin(g,*vit);
	  }
	}
	lastcoord = fragmentBegin(g,*vit)+fragmentLength(g,*vit);
	//it->first.first is CC labe
	repeatCC.insert(c); 
	vlist.push_back(currV);
	//std::cout << "Last coord:" << lastcoord << " component:" << c << " size:" << componentMap[c].size() << std::endl;
      }
    }
    //assert(runcount+1==runs.size());
        
    //A run is a list of CC
    //Create a new seqSet that contains all the runs
    
    //for(typename std::set<std::pair<int,int> >::iterator rit = runs.begin();rit != runs.end();++rit){
      //std::cout << "Seq: " << currentSeq 
      //	<< rit->first << "-"  << rit->second << " " 
      //	<< rit->second - rit->first << std::endl;
      //TString newseq=seqSet[i];
      //For each run, create a new seqset
      //appendValue(runSeqSet,seqSet[currentSeq]);
      //addVertex(runG);
    //}
    //std::cout << "Current seq run count " << runcount  << " " << " . Total runs " << runs.size() << std::endl;

  }
  //
  //Each run represents a copy of a duplicated region
  //We will determine the copies that need to be aligned and store them in an LCB.
  //Also stated, an LCB is a set of runs, where each run is a set of vertices in G(0).
  //
  //Determining the LCBs that represent duplications
  //
  //After the runs have been defined
  //an LCB is simply a list of all vertices in G(0) reachable in the run; ie. all members of the CCs in that run
  //
  //
  //To obtain the list of runs that comprise an LCB:
  //Build a graph G(l) where each node is a runidx and an edge connects any two runidx that share a ccidx from G(0)
  //To find the LCBs, simply determine the connected components in this graph G(l).  Each component defines corresponds to an LCB 
  //and the list of vertices is obtained from G(0) using ccidx at each node

  //foreach run in runs
  //  addvertex(run,g.l)
  //
  //foreach v1 in cc
  // foreach v2 in cc
  //  if(runmap[v1]!=runmap[v2])
  //   addedge(runv1,runv2,g.l)
  //
  //numcc = connected_components(rccmap,g.l)
  //lcbs.resize(numcc)
  //
  //foreach runv (g.l)
  // //save all vertices associated with
  // lcbidx <- rccmap[runv]
  // lcbs[lcbidx].push_back(runv)

  
  //Next, the vertices need to be updated to map to unique sequence ids for each run.
  //foreach lcb (lcbs)
  //  numruns <- lcbs[lcb].size()
  //  Built new set set
  //  Build a new graph 
  //  foreach run (lcbs[lcb])
  //     foreach v (vmap[run])
  //      seqIdMap[v] = runmap[v]
  //      lcbv.push_back(v)


  //At this point, runSeqSet should be populated
  //graph G(0) should be updated to reference sequence ids in runSeqSet
  //seqIdMap should be an identity map

  //retrieveLCBSegments()
  
  //Probably not necessary
  //foreach match (currmatches)
  //  fragment(match,0).seq = sequenceid[fragment(match,0).seq]
  //  fragment(match,1).seq = sequenceid[fragment(match,1).seq]

  //buildAlignmentGraph()

  TComponentMap runccmap;


  //Add edges between any 2 runs that are connected in G(0)

  //Edge iterator 
  typedef typename Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
  TEdgeIterator itE1(g);
  for(;!atEnd(itE1);++itE1){
    TEdgeDescriptor ed = *itE1;
    TVertexDescriptor source = getSource(ed);
    TVertexDescriptor target = getTarget(ed);
    if(runmap.find(source)!=runmap.end()
       && runmap.find(target)!=runmap.end()){
      if(runmap[source]!=runmap[target]){
	addEdge(runG,runmap[source],runmap[target]);
      }
      else{
	//std::cout << runmap[source] << " " << runmap[target] << " s:" << sequenceId(g,source) << " " << sequenceId(g,target) <<std::endl;
	assert(runmap[source]==runmap[target]);
	assert(sequenceId(g,source)==sequenceId(g,target));
      }
    }
  }

  int numcc = connected_components(runG,runccmap);
  std::cerr << "Determined " << numcc << " LCBs from a graph of runs: " << numVertices(runG) << std::endl;

  typedef std::vector<TVertexDescriptor> TLCB;
  std::vector<std::vector<TVertexDescriptor> > runs;
  std::vector<std::vector<TVertexDescriptor> > LCBs;
  //List of runs in an LCB
  runs.resize(numcc);
  //List of vertices in an LCB
  LCBs.resize(numcc);

  for(unsigned i=0;i<numVertices(runG);i++){
    int lcbidx = runccmap[i];
    //std::cout << "LCB " << lcbidx << " contains run " << i << std::endl;
    runs[lcbidx].push_back(i);
  }
  std::map<TSize,TSize> seqIdMap;
  //std::cout << "LCBs " << runs.size() << std::endl;
  for(unsigned j=0;j<runs.size();j++){
    int lcbidx = j;
    int numruns = runs[j].size();
    //std::cout << "LCB " << j << " runs " << numruns << std::endl;
    for(int k=0;k<numruns;k++){
      int runidx = runs[j][k];
      LCBs[lcbidx].insert(LCBs[lcbidx].end(),vrunmap[runidx].begin(),vrunmap[runidx].end());
      for(unsigned i=0;i<vrunmap[runidx].size();i++){
	TVertexDescriptor currV = vrunmap[runidx][i];
	//std::cout << "CurrV: " << currV << " runidx: " << runidx << " " << " seq:" << sequenceId(g,currV) << std::endl;
	g.data_fragment[currV].data_seq_id = runidx;
	//std::cout << "CurrV: " << currV << " runidx: " << runidx << " " << " seq:" << sequenceId(g,currV) << std::endl;
	assert(sequenceId(g,currV)==(unsigned)runidx);
      }
      seqIdMap[runidx] = runidx;
    }
  }

  //Set vertex orientation
  std::map<TVertexDescriptor,char> vertexOrientMap;
  typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
  TVertexIterator itV2(g);
  for(;!atEnd(itV2);goNext(itV2)){
    vertexOrientMap[*itV2] = '+';
  }

  //Trim graph and remove vertices that are not part of runs
  typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
  TVertexIterator itv(g);
  std::vector<TVertexDescriptor> removeV;
  for(;!atEnd(itv);goNext(itv)) {
    TVertexDescriptor currV = *itv;
    if(runmap.find(currV)==runmap.end()){
      removeV.push_back(currV);
    }
  }
  for(typename std::vector<TVertexDescriptor>::iterator vit = removeV.begin();vit!=removeV.end();++vit){
    removeVertex(g,*vit);
  }

  //At this point, runSeqSet should be populated
  //graph G(0) should be updated to reference sequence ids in runSeqSet
  //seqIdMap should be an identity map
  
  resize(runSeqSet,numVertices(runG));
  TNames sequenceRunsNames;
  //resize(sequenceRunsNames,numVertices(runG));
  for(unsigned int i=0;i<numVertices(runG);i++){
    std::string name(toCString(sequenceNames[seqIdxMap[i]]));
    std::string count(boost::lexical_cast<std::string>(i));
    appendValue(sequenceRunsNames,name+"_"+count);
  }
  //For tracking substrings 
  //std::map<TSize,unsigned int> offsets;
  //std::map<TSize,unsigned int> spanlens;
  //std::map<TSize,unsigned int> seqlens;
  //std::map<TSize,char> orients;
  typedef unsigned TSize2;
  std::map<TSize2,s_offset> offsets;
  
  //Copy links between set of vertices in LCB $lit
  //from Graph $g and store in $currmatches,$currscores,$currseqs
  
  blocksbycomponent.resize(numComponents);
  convertCC2Blocks(g,
		   componentall,
		   componentVertexMap,
		   blocksbycomponent,
		   dupintervals,
		   sequenceNames);
#ifdef SEQAN_PROFILE
    std::cerr << "Determining duplicated regions done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;
    std::cerr << "Aligning duplicated regions" << std::endl;
#endif
  unsigned int lcbid=0;
  std::fstream strmmaf;
  std::string outfile(msaOpt.outfile);
  strmmaf.open(std::string(outfile+".dups.maf").c_str(), std::ios_base::out | std::ios_base::trunc);
  _streamWrite(strmmaf,"##maf version=1 scoring=mugsy");
  //std::cout << "Iterating over LCBs" << std::endl;
  for(typename std::vector<std::vector<TVertexDescriptor> >::iterator lit = LCBs.begin();lit!=LCBs.end();lit++){
    //Matches, scores, seqs, ids for current LCB
    typedef String<Fragment<> > TFragmentString;
    typedef typename Value<TScore>::Type TScoreValue;
    typedef String<TScoreValue> TScoreValues;
    typedef String<TAlphabet> TSequence;
    TFragmentString currmatches;
    TScoreValues currscores;
    StringSet<TSequence, Owner<> > currseqs;
    std::set<unsigned int> curridset; 
    TNames currnameSet;

    std::set<TVertexDescriptor> coveredSet;
    std::vector<std::vector<TVertexDescriptor> > vseqs;
    retrieveLCBSegments(g,
			runSeqSet,
			seqSet,
			seqIdxMap,
			vertexOrientMap,
			lit,
			++lcbid,
			sequenceRunsNames,
			currseqs,
			currmatches,
			currscores,
			currnameSet,
			offsets,
			coveredSet,
			vseqs,
			boost::lexical_cast<unsigned int>(msaOpt.minlength));

    
    
    //std::cout << "Retrieving LCB segments for LCB " << lcbid << std::endl; 
    if(length(currseqs)>1 && length(currmatches)>0){
      TGraph currG(currseqs);
      buildAlignmentGraph(currmatches, currscores, currG, FractionalScore());
      //Double check edge weights
      typedef typename Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
      TEdgeIterator itE(currG);
      //Undo Hack that stores reverse complement matches using
      //negative edge weights
      for(;!atEnd(itE);goNext(itE)){
	if(cargo(value(itE))<0){
	  cargo(value(itE)) = cargo(value(itE))*-1;
	}
      }
      
      typedef double TDistanceValue;
      Graph<Tree<TDistanceValue> > currguideTree;
      Graph<Tree<TDistanceValue> > seqguideTree;
      std::map<std::string, Graph<Tree<TDistanceValue> > > seqguideTrees;

      //Build guide tree using current list of seqs
      for(unsigned int i=0;i<length(currseqs);i++){
	//std::cout << i << " " << length(currseqs[i])  << " " << ((curridset.find(i)!=curridset.end()) ? 1 : 0) << " " << currnameSet[i]<< std::endl;
	curridset.insert(i);//force inclustion of this seq in building the guide tree
      }
      getGuideTree(currseqs,curridset,seqguideTrees,currguideTree);

      typedef Fragment<> TFragment;
      typedef String<TAlphabet> TSequence;
      
      std::cout << "Aligning LCB " << lcbid << " with " << length(currseqs) << std::endl;
      assert(curridset.size()>0);
      TGraph currgOut(currseqs);
      s_score sscores = alignSingleLCB(currG,
				       currgOut,
				       lcbid,
				       currseqs,
				       currguideTree,
				       msaOpt);
      //Write MAF format
      std::vector<unsigned int> curroffsets;
      std::vector<unsigned int> currspanlens;
      std::vector<unsigned int> currseqlens;
      std::vector<char> currorients;
      assert(length(currseqs)==length(currnameSet));
      currorients.resize(length(currnameSet));
      currseqlens.resize(length(currnameSet));
      currspanlens.resize(length(currnameSet));
      curroffsets.resize(length(currnameSet));
      //TODO, refactor using a id map
      for(TSize currrow = 0; currrow<length(currnameSet); ++currrow) {
	for(TSize row = 0; row<length(sequenceRunsNames); ++row) {
	  if(currnameSet[currrow]==sequenceRunsNames[row]){
	    curroffsets[currrow] = offsets[row].offset;
	    currspanlens[currrow] = offsets[row].spanlen;
	    currseqlens[currrow] = offsets[row].seqlen;
	    currorients[currrow] = offsets[row].orient;
	    //reset name
	    currnameSet[currrow] = sequenceNames[seqIdxMap[row]];
	  }
	}
      } 
      saveInterval(dupintervals,
		   currnameSet,
		   curroffsets,
		   currspanlens,
		   currseqlens,
		   currorients,
		   lcbid,
		   true);
      //mafformat defined in refinement/graph_impl_align.h
      write(strmmaf,currgOut,currnameSet,MafFormat(),curroffsets,currspanlens,currseqlens,currorients,"label=d"+boost::lexical_cast<std::string>(lcbid));
      
      strmmaf.flush();
    }
  }
  strmmaf.close();
  //Alignments of many duplicated regions tend to be fragmented on first pass
  //Consider running refinement by default
  //refineMSA(std::string(outfile+".dups.maf").c_str(),msaOpt);
#ifdef DEBUGGING
  //TODO
  //Resolve repetitive clusters here
  //(1)break edges with weak support from adjacent matches
  //(2)determine mincut on repeatitive clusters
  typedef typename Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;

  std::fstream dotstrm;
  dotstrm.open("refinegraphpos.dot", std::ios_base::out | std::ios_base::trunc);
  dotstrm << "graph g{" << std::endl;
  typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
  TVertexIterator it(g);
  for(;!atEnd(it);goNext(it)) {
    dotstrm << *it << " [label=\"" << *it << " S" << sequenceId(g,*it) << ","<<fragmentBegin(g,*it) << ","<<fragmentLength(g,*it) << "\"];" << std::endl;
  }
#endif
#ifdef SEQAN_PROFILE
    std::cerr << "Aligning duplicated regions done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;
#endif
}


//////////////////////////////////////////////////////////////////////////////////

inline void
_addVersion(CommandLineParser& parser) {
	::std::string rev = "$Revision: 4637 $";
	addVersionLine(parser, "Version 1.00 (10 Oct 2009) Revision: " + rev.substr(11, 4) + "");
}

//////////////////////////////////////////////////////////////////////////////////

template <typename TSeqSet, typename TNameSet>
bool _loadSequences(TSeqSet& sequences, 
		    TNameSet& fastaIDs, 
		    TNameSet& genomes,
		    const char *fileName)
{
  assert(length(genomes)==0);
	MultiFasta multiFasta;
	if (!open(multiFasta.concat, fileName, OPEN_RDONLY)) return false;
	AutoSeqFormat format;
	guessFormat(multiFasta.concat, format);	
	split(multiFasta, format);
	unsigned seqCount = length(multiFasta);
	resize(sequences, seqCount, Exact());
	resize(fastaIDs, seqCount, Exact());
	resize(genomes, seqCount, Exact());
	unsigned skippedseqCount = 0;
	for(unsigned i = 0; i < seqCount; ++i) 
	  {
	    char seqname[100],genomename[100];
	    std::string idline;
	    assignSeqId(idline, multiFasta[i], format);
	    int matches = sscanf(idline.c_str(),"%s %s",seqname,genomename);
	    if(matches==2){
	      fastaIDs[i]=seqname;
	      genomes[i]=genomename;
	    }
	    else{
		assignSeqId(fastaIDs[i], multiFasta[i], format);
		assignSeqId(genomes[i], multiFasta[i],format);
	    }
	    assignSeq(sequences[i], multiFasta[i], format);
	    //SVA check for bad inputs here < kmer size
	    if(length(sequences[i])<3){
	      skippedseqCount++;
	    }
	}
	if(skippedseqCount>0){
	  clear(sequences);
	  clear(fastaIDs);
	  clear(genomes);
	  seqCount = length(multiFasta)-skippedseqCount;
	  //std::cerr << "Updated seqCount " << seqCount << ". Skipping" << skippedseqCount << std::endl;
	  resize(sequences, seqCount, Exact());
	  resize(fastaIDs, seqCount, Exact());
	  resize(genomes, seqCount, Exact());
	  unsigned sidx=0;
	  String<char> testseq;
	  unsigned oseqCount = length(multiFasta);
	  for(unsigned i = 0; i < oseqCount; ++i) 
	  {
	    assignSeq(testseq, multiFasta[i], format);
	    //SVA check for bad inputs here < kmer size
	    if(length(testseq)>=3){
	      char seqname[100],genomename[100];
	      std::string idline;
	      assignSeqId(idline, multiFasta[i], format);
	      int matches = sscanf(idline.c_str(),"%s %s",seqname,genomename);
	      if(matches==2){
		fastaIDs[sidx]=seqname;
		genomes[sidx]=genomename;
	      }
	      else{
		assignSeqId(fastaIDs[sidx], multiFasta[i], format);
		assignSeqId(genomes[sidx], multiFasta[i],format);
	      }
	      assignSeq(sequences[sidx], multiFasta[i], format);
	      sidx++;
	    }
	    else{
	      std::cerr << "Skipping sequence of length " << length(testseq) << std::endl;
	    }
	  }
	  assert(sidx==seqCount);
	}
	return (seqCount > 0);
}

template<typename TAlphabet, typename TScore>
inline void
customizedMsaAlignment(MsaOptions<TAlphabet, TScore> const& msaOpt) {
	typedef String<TAlphabet> TSequence;
	StringSet<TSequence, Owner<> > sequenceSet;
	StringSet<String<char> > sequenceNames;
	StringSet<String<char> > genomeNames;
	_loadSequences(sequenceSet, sequenceNames, genomeNames, msaOpt.seqfile.c_str());
	assert(length(sequenceNames)==length(sequenceSet));
#ifdef DEBUGGING
	for(unsigned int j = 0; j<length(sequenceNames); ++j) {
	  std::cout << j << " " << sequenceNames[j] << std::endl;
	  assert(value(sequenceNames,j)==sequenceNames[j]);
	}
#endif
	// Alignment of the sequences
	Graph<Alignment<StringSet<TSequence, Dependent<> >, void, WithoutEdgeId> > gAlign;
	typedef unsigned int TSize;
	TSize gidx=0;
	std::map<String<char>,TSize> genomeIdx;
	String<TSize> genomeIndices;
	//Convert Names to indicies
	for(TSize i=0;i<length(genomeNames);++i){
	  TSize cidx;
	  if(genomeIdx.find(genomeNames[i])==genomeIdx.end()){
	    genomeIdx[genomeNames[i]]=gidx;
	    cidx=gidx;
	    gidx++;
	  }
	  else{
	    cidx=genomeIdx[genomeNames[i]];
	  }
	  appendValue(genomeIndices,cidx);
	}
	
	// Calc MSA
	 //Aligned intervals, used to determine remaining segments unaligned
	typedef iloc TLoc;
	std::map<String<char>,std::vector<TLoc> > aintervals;
	if(msaOpt.duplications == "true"){
	  findDuplications(gAlign, sequenceSet, sequenceNames, genomeIndices, aintervals, msaOpt);
	}

	//if(msaOpt.partition >= 2){
	  //Testing new code
	  //Prototype only
	  //multipassprog_wholeGenomeAlignment(gAlign, sequenceSet, sequenceNames, genomeIndices, msaOpt);
	  //std::cerr << "Bad partition parameter " << msaOpt.partition << std::endl;
	  //exit(1);
	//}
	//else{
	singlepass_wholeGenomeAlignment(gAlign, sequenceSet, sequenceNames, genomeIndices, aintervals, msaOpt);
	//}
		
	// Alignment output
	if (msaOpt.outputFormat == 0) {
		FILE* strmWrite = fopen(msaOpt.outfile.c_str(), "w");
		write(strmWrite, gAlign, sequenceNames, FastaFormat());
		fclose(strmWrite);
	} else if (msaOpt.outputFormat == 1) {
		FILE* strmWrite = fopen(msaOpt.outfile.c_str(), "w");
		write(strmWrite, gAlign, sequenceNames, MsfFormat());
		fclose(strmWrite);
	}

}

//////////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TScore, typename TSc>
inline void
_setMatchScore(MsaOptions<TAlphabet, TScore>&, TSc) {
	// No operation
}

//////////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TScore, typename TSc>
inline void
_setMismatchScore(MsaOptions<TAlphabet, TScore>&, TSc) {
	// No operation
}

//////////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TSc>
inline void
_setMatchScore(MsaOptions<TAlphabet, Score<int, Simple> >& msaOpt, TSc msc) {
	msaOpt.sc.data_match = msc;
}

//////////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TSc>
inline void
_setMismatchScore(MsaOptions<TAlphabet, Score<int, Simple> >& msaOpt, TSc mmsc) {
	msaOpt.sc.data_mismatch = mmsc;
}


//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
template<typename TConfigOptions, typename TScore>
inline void
evaluateAlignment(TConfigOptions const& cfgOpt, TScore const& scType, Dna5) {
  std::fstream strmmaf;
  //FILE * strmmafrefined;

  struct mafFile *mf;
  mf = mafOpen(cfgOpt.infile.c_str(), 0);
  struct mafAli *a, *A, *last_a;
  struct mafComp *c;
  A = last_a = NULL;
  while ((a = mafNext(mf)) != NULL) {
    if ((c = a->components) == NULL)
      assert(false);//fatal("empty maf entry");
    if (last_a == NULL)
      A = a;
    else
      last_a->next = a;
    last_a = a;
  }
  if(A==NULL){
#ifdef DEBUGGING
    std::cout << "can't find any alignments" << std::endl;
#endif
  }
  else{
    int lcbid=0;
    char chrName[200], species_name[200];
    int COL_WIDTH=60;
    long unsigned int totalscore=0;
    long unsigned int totallen=0;
    for (a = A; a != NULL; a = a->next) {
      int ncol = a->textSize;
      std::ostringstream tmpgraph;
      tmpgraph << "MUGTMP" << getpid() << "_" << ++lcbid;
      std::fstream strmfsa;
      std::string fname(tmpgraph.str());
      fname = "/tmp/"+fname + ".eval.fsa";
      strmfsa.open(fname.c_str(), std::ios_base::out | std::ios_base::trunc);
      for(c=a->components; c!=NULL; c=c->next) {
	parseSrcName(c->src, species_name, chrName);
	//Write FASTA
	strmfsa << ">" << c->src << std::endl ;
	int col=0;
	int j=0;
	for (col = j = 0; j < ncol; ++j) {
	  strmfsa << c->text[j];
	  ++col;
	  if (col == COL_WIDTH) {
	    strmfsa << std::endl;
	    col = 0;
	  }
	}
	if (col != 0){
	  strmfsa << std::endl;
	}
      }
      strmfsa.close();
        
      typedef typename Value<TScore>::Type TScoreValue;
      typedef String<Dna5> TSequence;
      typedef typename Size<TSequence>::Type TSize;
      typedef String<char> TName;
      StringSet<TSequence, Owner<> > origStrSet;
      StringSet<TName> names;
      
      // Read the sequences
      std::fstream strm;
      strm.open(fname.c_str(), std::ios_base::in | std::ios_base::binary);
      read(strm,origStrSet,names,FastaAlign());	
      strm.close();
      
      // Make a dependent StringSet
      typedef StringSet<TSequence, Dependent<> > TDepSequenceSet;
      TDepSequenceSet strSet(origStrSet);
      
      // Read the alignment
      typedef String<Fragment<> > TFragmentString;
      String<TScoreValue> scores;
      TFragmentString matches;
      std::fstream strm_lib;
      strm_lib.open(fname.c_str(), std::ios_base::in | std::ios_base::binary);
      read(strm_lib,matches, scores, names, FastaAlign());	
      strm_lib.close();
      unlink(fname.c_str());
      // Build the alignment graph
      typedef Graph<Alignment<TDepSequenceSet, TSize> > TGraph;
      TGraph g(strSet);
      buildAlignmentGraph(matches, g, FrequencyCounting() );
      
      // Print the scoring information
      TScoreValue gop = scType.data_gap_open;
      TScoreValue gex = scType.data_gap_extend;
      std::cout << "Scoring parameters:" << std::endl;
      std::cout << "*Gap opening: " << gop << std::endl;
      std::cout << "*Gap extension: " << gex << std::endl;
      std::cout << "*Scoring matrix: " << std::endl;
      TSize alphSize = ValueSize<Dna5>::VALUE;
      std::cout << "   ";
      for(TSize col = 0; col<alphSize; ++col) std::cout << Dna5(col) << ',';
      std::cout << std::endl;
      for(TSize row = 0; row<alphSize; ++row) {
	for(TSize col = 0; col<alphSize; ++col) {
	  if (col == 0) std::cout << Dna5(row) << ": ";
	  std::cout << score(scType, Dna5(row), Dna5(col));
	  if (col < alphSize - 1) std::cout << ',';
	}
	std::cout << std::endl;
      }
      std::cout << std::endl;
      
      // Print the alignment information
      TSize numGapEx = 0;
      TSize numGap = 0;
      TSize numPairs = 0;
      TSize alignLen = 0;
      String<TSize> pairCount;
      TScoreValue alignScore = alignmentEvaluation(g, scType, numGapEx, numGap, numPairs, pairCount, alignLen);
      totalscore+=alignScore;
      totallen+=alignLen;
      std::cout << "Alignment Score: " << alignScore << std::endl;
      std::cout << "Alignment Length: " << alignLen << std::endl;
      std::cout << "#Match-Mismatch pairs: " << numPairs << std::endl;
      std::cout << "Score contribution by match-mismatch pairs: " << (alignScore - (((TScoreValue) numGap * gop) + ((TScoreValue) numGapEx * gex))) << std::endl;
      std::cout << "#Gap extensions: " << numGapEx << std::endl;
      std::cout << "Score contribution by gap extensions: " << ((TScoreValue) numGapEx * gex) << std::endl;
      std::cout << "#Gap openings: " << numGap << std::endl;
      std::cout << "Score contribution by gap openings: " << ((TScoreValue) numGap * gop) << std::endl;
      std::cout << std::endl;
      std::cout << "#Pairs: " << std::endl;
      std::cout << "   ";
      for(TSize col = 0; col<alphSize; ++col) std::cout << Dna5(col) << ',';
      std::cout << std::endl;
      for(TSize row = 0; row<alphSize; ++row) {
	for(TSize col = 0; col<alphSize; ++col) {
	  if (col == 0) std::cout << Dna5(row) << ": ";
	  std::cout << value(pairCount, row * alphSize + col);
	  if (col < alphSize - 1) std::cout << ',';
	}
	std::cout << std::endl;
	
      }
      /*
      struct mafAli *nexta;
      for (a = A; a != NULL; a = nexta) {
	nexta=a->next;
	mafAliFree(&a);
      }
      */
    }
    mafFileFree(&mf);
    std::cout << "Total alignment score: " << totalscore << std::endl;
    std::cout << "Total alignment length: " << totallen << std::endl;
  }
}

template<typename TAlphabet, typename TScore>
inline void
_initMsaParams(CommandLineParser& parser, TScore& scMat) {
	
	// Msa configuration
	MsaOptions<TAlphabet, TScore> msaOpt;
	
	// Set main options
	getOptionValueLong(parser, "seq", msaOpt.seqfile);
	getOptionValueLong(parser, "outfile", msaOpt.outfile);
	// MUGSY specific options
	getOptionValueLong(parser, "distance", msaOpt.distance);
	getOptionValueLong(parser, "minlength", msaOpt.minlength);
	getOptionValueLong(parser, "refine", msaOpt.refine);
	getOptionValueLong(parser, "duplications", msaOpt.duplications);
	getOptionValueLong(parser, "unique", msaOpt.unique);
	getOptionValueLong(parser, "allownestedlcbs", msaOpt.allownestedlcbs);
	getOptionValueLong(parser, "anchorwin", msaOpt.anchorwin);
	getOptionValueLong(parser, "blockfile", msaOpt.blockfile);
	getOptionValueLong(parser, "segmentation", msaOpt.segmentation);
	if(msaOpt.segmentation != "none" && msaOpt.segmentation != "enredo" && msaOpt.segmentation != "mercator"){
	  msaOpt.segmentation = "mugsy";
	}
	String<char> optionVal;
	getOptionValueLong(parser, "format", optionVal);
	if (optionVal == "maf") msaOpt.outputFormat = 0;
	else if (optionVal == "msf") msaOpt.outputFormat = 1;

	unsigned int beg = 0;
	::std::string tmpVal;

	if (beg != tmpVal.length())
		appendValue(msaOpt.libfiles, tmpVal.substr(beg, tmpVal.length() - beg));	
	getOptionValueLong(parser, "aln", tmpVal);
	beg = 0;
	for(unsigned int i = 0; i<tmpVal.length(); ++i) {
		if (tmpVal[i] == ',') {
			appendValue(msaOpt.alnfiles, tmpVal.substr(beg, i - beg));
			beg = i + 1;
		}
	}
	if (beg != tmpVal.length())
		appendValue(msaOpt.alnfiles, tmpVal.substr(beg, tmpVal.length() - beg));

	// Set scoring options
	msaOpt.sc = scMat;
	getOptionValueLong(parser, "gop", msaOpt.sc.data_gap_open);
	getOptionValueLong(parser, "gex", msaOpt.sc.data_gap_extend);
	int msc = 0;
	getOptionValueLong(parser, "msc", msc);
	_setMatchScore(msaOpt, msc);
	int mmsc = 0;
	getOptionValueLong(parser, "mmsc", mmsc);
	_setMismatchScore(msaOpt, mmsc);

	// Set guide tree options

	if (optionVal == "nj") msaOpt.build = 0;
	else if (optionVal == "min") msaOpt.build = 1;
	else if (optionVal == "max") msaOpt.build = 2;
	else if (optionVal == "avg") msaOpt.build = 3;
	else if (optionVal == "wavg") msaOpt.build = 4;

	// Set alignment evaluation	options
	getOptionValueLong(parser, "infile", msaOpt.infile);

	// Check if any segment-match generation procedure is selected, otherwise set the default
	if ((empty(msaOpt.alnfiles)) && (empty(msaOpt.method))) {
		appendValue(msaOpt.method, 0);
		appendValue(msaOpt.method, 1);
	}

	// Evaluation mode?
	if (isSetLong(parser, "infile")) {
	  if(length(msaOpt.refine) > 0 && msaOpt.refine != "colinear"){ //Refinement mode
	    refineMSA(msaOpt.infile.c_str(),msaOpt);
	  }
	  else {
	    //typedef typename Value<TScore>::Type TScoreValue;
	    //TScore scType(boost::lexical_cast<int>(value(msaOpt, "msc")),
	    //	  boost::lexical_cast<int>(value(msaOpt, "mmsc")),-1 * boost::lexical_cast<int>(value(msaOpt, "gex")),-1 * boost::lexical_cast<int>(value(msaOpt, "gop")));
	    evaluateAlignment(msaOpt, msaOpt.sc, Dna5() );
	    //evaluateAlignment(msaOpt);
	  }
	} else { // or alignment mode?
	  if (!isSetLong(parser, "seq")) { 
	    shortHelp(parser, std::cerr);	// print short help and exit
	    exit(0);
	  }
	  customizedMsaAlignment(msaOpt);
	}
}



inline void
_initScoreMatrix(CommandLineParser& parser, Dna5 const) {
	String<char> matrix;
	getOptionValueLong(parser, "matrix", matrix);
	if (isSetLong(parser, "matrix")) {
		Score<int, ScoreMatrix<> > sc;
		loadScoreMatrix(sc, matrix);
		_initMsaParams<Dna5>(parser, sc);
	} else {
		Score<int> sc;
		_initMsaParams<Dna5>(parser, sc);
	}
}



int main(int argc, const char *argv[]){
#ifdef TIMING 
  time(&now);
  lasttime=now;
#endif
  //////////////////////////////////////////////////////////////////////////////
  // Command line parsing
  //////////////////////////////////////////////////////////////////////////////
  std::string versionstring = std::string("1.3");  
  // Set the keys
  CommandLineParser parser;
  _addVersion(parser);
  
  addTitleLine(parser, "*************************************************");
  addTitleLine(parser, "* mugsyWGA                                      *");
  addTitleLine(parser, "* v"+versionstring+"                                          *");
  addTitleLine(parser, "* Multiple whole genome aligner                 *");
  addTitleLine(parser, "* using graph based LCB identification          *");
  addTitleLine(parser, "* and Seqan::TCoffee                            *");
  addTitleLine(parser, "*************************************************");

  addUsageLine(parser, "-seq <multi-FASTA sequence file> -aln <Aligned pairwise FASTA library> [-distance <LCB chaining distance>] [-minlength <LCB minimum length>] [Other options]");

  //Many config options lifted from seqan::tcoffee
  addSection(parser, "Main Options:");
  addOption(parser, addArgumentText(CommandLineOption("s", "seq", "multi-FASTA file with all input sequences. For draft genomes, FASTA headers should be in the form >seqname genomename.", OptionType::String), "<FASTA Sequence File>"));
  addOption(parser, addArgumentText(CommandLineOption("al", "aln", "Library of pairwise alignments. Aligned multi-FASTA format (XMFA)", OptionType::String), "<File1>,<File2>,..."));
  addOption(parser, addArgumentText(CommandLineOption("o", "outfile", "output filename prefix", (int)OptionType::String, "outfile"), "<Filename>"));
  addOption(parser, addArgumentText(CommandLineOption("distance", "distance", "LCB chaining distance", (int)OptionType::String,"1000"), "<String>"));
  addOption(parser, addArgumentText(CommandLineOption("minlength", "minlength", "Minimum LCB segment span", (int)OptionType::String,"100"), "<Int>"));
  addOption(parser, addArgumentText(CommandLineOption("unique", "unique", "Report unique regions", OptionType::String,"true"), "[true|false]"));
  addOption(parser, addArgumentText(CommandLineOption("duplications", "duplications", "Report duplications. Requires a second alignment file of pairwise duplications is passed to --aln. ", OptionType::String,"false"), "[true|false]"));
  
  addSection(parser, "Other Options:");

  addOption(parser, addArgumentText(CommandLineOption("f", "format", "output format", (int)OptionType::String, "maf"), "[maf | msf]"));
  addOption(parser, addArgumentText(CommandLineOption("anchorwin", "anchorwin", "bp window to consider for collapsing anchors", (int)OptionType::Int,0), "<Int>"));

  //synchain-mugsy can return overlapping and nested synteny blocks with the extent determined by --distance
  //allownestedlcbs=false ensures each multi-genome anchor contributes to exactly one LCB; the longest LCB spanning the anchor
  //The LCBs are sorted by length in descending order. Each anchor is
  //removed from the anchor graph as soon as it is aligned in an LCB.
  addOption(parser, addArgumentText(CommandLineOption("allownestedlcbs", "allownestedlcbs", "allow anchors to contribute to multiple LCBs. Default=false", OptionType::String,"false"), "[true|false]"));

  addOption(parser, addArgumentText(CommandLineOption("refine", "refine", "refinement method: mugsy,fsa,pecan,mlagan", OptionType::String), "<String>"));
  //addOption(parser, addArgumentText(CommandLineOption("poscorewindow", "psw", "posscorewindow", (int)OptionType::Int,1000), "<Int>"));
  //addOption(parser, addArgumentText(CommandLineOption("possharedcutoff", "pscut", "possharedcutoff", (int)OptionType::Double,(double)0.1), "<Int>"));

  addOption(parser, addArgumentText(CommandLineOption("segmentation", "segmentation", "Segmentation method. mugsy,enredo,mercator", OptionType::String), "<String>"));
  addOption(parser, addArgumentText(CommandLineOption("blockfile", "blockfile", "Bypass segmentation and use this output file from synchain-mugsy", OptionType::String), "<String>"));

  addSection(parser, "Scoring Options:");
  addOption(parser, addArgumentText(CommandLineOption("g", "gop", "gap open penalty", (int)OptionType::Int, -13), "<Int>"));
  addOption(parser, addArgumentText(CommandLineOption("e", "gex", "gap extension penalty", (int)OptionType::Int, -1), "<Int>"));
  addOption(parser, addArgumentText(CommandLineOption("ma", "matrix", "score matrix", (int)OptionType::String, "Blosum62"), "<Matrix file>. Ignored."));
  addOption(parser, addArgumentText(CommandLineOption("ms", "msc", "match score", (int)OptionType::Int, 5), "<Int>"));
  addOption(parser, addArgumentText(CommandLineOption("mm", "mmsc", "mismatch penalty", (int)OptionType::Int, -4), "<Int>"));

  addSection(parser, "Guide Tree Options:");
  //addOption(parser, addArgumentText(CommandLineOption("u", "usetree", "tree filename", OptionType::String), "<Newick guide tree>"));
  addOption(parser, addArgumentText(CommandLineOption("b", "build", "tree building method for progressive aln", (int)OptionType::String, "nj"), "[nj, min, max, avg, wavg]"));
  addHelpLine(parser, "nj = Neighbor-joining");
  addHelpLine(parser, "min = UPGMA single linkage");
  addHelpLine(parser, "max = UPGMA complete linkage");
  addHelpLine(parser, "avg = UPGMA average linkage");
  addHelpLine(parser, "wavg = UPGMA weighted average linkage");
  addHelpLine(parser, "Neighbor-joining creates an");
  addHelpLine(parser, "  unrooted tree. We root that tree");
  addHelpLine(parser, "  at the last joined pair.");
  // Alignment evaluation	
  addSection(parser, "Alignment Evaluation Options:");
  addOption(parser, addArgumentText(CommandLineOption("i", "infile", "alignment file", OptionType::String), "<FASTA alignment file>"));
  
  if (argc == 1)
    {
      shortHelp(parser, std::cerr);	// print short help and exit
      return 0;
    }

  bool exitrun=false;
  if (!parse(parser, argc, argv, ::std::cerr)) exitrun=true;
  if (isSetLong(parser, "help") || isSetLong(parser, "version")) exitrun=false;	// print help or version and exit
    

  char * mugsyinstallstr = std::getenv("MUGSY_INSTALL");
  if(mugsyinstallstr==NULL || strlen(mugsyinstallstr)==0){
    std::cerr << "ERROR: Environment variable MUGSY_INSTALL must be set to the installation directory for mugsy" << std::endl;
    exit(1);
  }
  assert(mugsyinstallstr != NULL);
  std::string mugsyinstall = std::string(mugsyinstallstr);
  assert(mugsyinstall.length()>0);
#ifdef DEBUGGING
  std::cerr << "Using MUGSY_INSTALL=" << mugsyinstall << std::endl;
#endif
  //Check for chaining executable
  struct stat st;
  if(stat(std::string(mugsyinstall+"/synchain-mugsy").c_str(),&st) == 0){
    //present
  }
  else{
    std::cerr << "ERROR: MUGSY_INSTALL/synchain-mugsy not found. check installation at MUGSY_INSTALL=" << mugsyinstall << std::endl;
    exitrun=true;
  }
  
  if(exitrun){
    return 1;
  }
  // Basic command line options
  String<char> alphabet = "dna";
  // Initialize scoring matrices
  _initScoreMatrix(parser, Dna5());
  return 0;
}


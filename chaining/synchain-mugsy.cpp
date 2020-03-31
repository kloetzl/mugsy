//USAGE:mugsy-chaining max-distance min-lcblen [min-lcblenstats] < anchors.projection
//

//Mugsy chaining algorithm to partition a graph of mult-genome anchors
//into collinear "syntenic" segments

//
//Input 
//----- 
//Projection format is 
//anchor1 anchor2 seqindex dist genomeindex orient1 orient2 beg1 end1 beg2 end2
//eg
//0 1 0 0 0 + + 0 196 196 15348
//1 3 0 1 0 + + 196 15348 15349 20373


//The anchor graph is an directed graph where each vertex is a
//multi-genome anchor and each edge connects adjacent anchors on one
//or more genomes.  This input projection already should list anchors
//that are adjacent on a genome within distance $max-distance.  The
//anchor graph will be built such that edges are stored for adjacent
//anchors in at least one genome. 

//A series of heuristics is applied to identify paths in the graph
//that correspond to collinear regions ignoring micro-rearrangments <
//max-distance.

//The regions may be overlapping with the degree of overlap determined
//by max-distance

//General outline

//(1)Build anchor graph

//(2)Initial clustering

//(2.1) Identify vertices with more than 2 edges and mask all incident
//edges in the graph. These are syntenic breakpoints. Some will be
//micro- events that we ignore later

//(2.2) Calculate connected components. The remaining edges correspond
//to vertices with exactly two vertices and comprise runs of synteny.

//(2.3) Run mincut to break paths that traverse breakpoints. Edges
//indicate on synteny on some genomes but do not ensure all incident
//anchors are syntenic. We use a maxflow-mincut procedure to determine
//which edges to break such that the LCBs respect max-distance and do
//not include inversions.

//(2.4) Merge adjacent LCBs. The procedures of (2.2) and (2.3) will
//over-parition the graph. Merge adjacent LCBs that have compatible
//anchor orientations and respect max-distance

//(2.5) Mask short LCBs after merge. LCBs < minlen after merging are
//masked from the graph. Next, the vertices are projected along each
//of the member sequences and additional edges are added to the anchor
//graph. The clustering of (2.1) and (2.2) is repeated to identify a
//new set of LCBs.  This step allows for ignoring short LCBs that may
//be breaking synteny.

//(2.6) Run mincut to restore invariants.

//(2.7) Merge

//At this step the LCBs are . Two additional iterations of masking short LCBs and merging are run to try to cluster additional bps.


//
//S. Angiuoli - UMD CS, 2009

#define NDEBUG 
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string>
#include <fstream>
#include <vector>
#include <queue>
#include <list>
#include <bitset>
#include <algorithm>
#include <ext/hash_set> //__gnu_cxx namespace
//#include <tr1/unordered_set>

#include <boost/graph/graphviz.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/adjacency_matrix.hpp>
#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>

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
#include <boost/property_map/property_map.hpp>
#include <boost/graph/breadth_first_search.hpp>

#include <boost/graph/edmonds_karp_max_flow.hpp>
#include <boost/graph/boykov_kolmogorov_max_flow.hpp>
#include <boost/graph/push_relabel_max_flow.hpp>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/read_dimacs.hpp>
#include <boost/graph/graph_utility.hpp>

// Archivers
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp> 

#include <boost/config.hpp>
#include <boost/pending/queue.hpp>

using namespace boost;
using namespace std;

//Maximum number of input genomes 
//There is no limit on the number of sequences per genome
//Used to set size of std::bitset<> only
//TODO, replace with boost::dynamic_bitset to avoid setting a limit
#define MAXGENOMES 256

//Print LCB stats
#define LCBSTATS 1
//Print timings for subsets
#define TIMING 1
//Use max-flow,min-cut 
#define CALCFLOW
#define CUTLCBEDGESONLY



//Debug creates 
//#define DEBUG 1
// print sequences and coords in graphviz output
//#define PRINTSEQS 
// draw flow network, cannot be combined with printseqs
//#define PRINTFLOW


//Defining this option removes all edges labelled in a single sequence
//only. Such edges represent can link to unaligned or non-syntenic
//regions

#define TRIMEDGES

//This is useful for simplying the graph and improves performance for draft genomes but changes the
//the algorithm 
//When defined, sequence specific indels > distance parameter will break blocks automatically
//When undefined, such indels are broken only during mincut


//Mugsy codes 
#include "graph.h"
#include "filters.h"
#include "file.h"
#include "lcbchecks.h" //isLabelCollinear,isLabelMaxGap,checkLCBGaps,checkLCBOrient,sameLabel,sameOrient,setLCBOrient
#include "mincut.h" //breakLCBmincutconnect


//variables for testing

//Define to store and print edge labels w/ distances
//in dot output
//Undefine for release to save space
//#define STORE_EDGE_LABELS 0
//#define V_DEBUG 0

 
//Remove misoriented vertices from an LCB
template<typename TGraph>
void fixMisOrientedLCBs(TGraph &g, 
			LCB & lcb,
			VertexSet &maskedLCBs,
			EdgeSet &maskedEdges,
			SequenceGenomeMap & sequence2genome){
#ifdef DEBUG
    std::cerr << "Trimming misoriented vertices in  lcb with " << lcb.size() << " vertices" << std::endl;
#endif
  std::vector<Vertex> badV;
  setLCBOrient(g,lcb,badV,sequence2genome);
  std::vector<Vertex>::iterator it,it_end;
  it_end = badV.end();
  for(it=badV.begin();it!=it_end;++it){
    maskedLCBs.insert(*it);
    typename graph_traits<TGraph>::out_edge_iterator ei, edge_end;
    typename graph_traits<TGraph>::in_edge_iterator ei2, edge_end2;
    tie(ei,edge_end) = out_edges(*it,g);
    for(;ei!=edge_end;++ei){
      maskedEdges.insert(std::make_pair(source(*ei,g),target(*ei,g)));
      put(edge_category,g,*ei,CYAN);
    }
    tie(ei2,edge_end2) = in_edges(*it,g);
    for(;ei2!=edge_end2;++ei2){
      maskedEdges.insert(std::make_pair(source(*ei2,g),target(*ei2,g)));
      put(edge_category,g,*ei2,CYAN);
    } 
  }
}



template<typename TGraph, typename TGraphB, typename TLCBMap>
void updateAdjacency(TGraph &g,
		     TGraphB &baseg,
		     std::set<Label> &seqidxSet,
		     VertexLabelIntervalMap &coordinates, 
		     TLCBMap & lcborientmap,
		     unsigned int distance,
		     EdgeSet&maskedEdges,
		     std::vector<int> &  ccvmap,		
		     std::vector<LCB> & componentMap,
		     SequenceGenomeMap & sequence2genome){
  //Graph of LCBs filtered by sequences 
  typedef typename property_map<TGraph, vertex_vlabelmask_t>::type VertexLabelMask; 
  typedef typename property_map<TGraph, vertex_orientmask_t>::type VertexOrientMask;
  typedef typename property_map<TGraph, edge_labelmask_t>::type EdgeLabelMask;
  
  typename property_map < TGraph, vertex_orientmask_t >::type orientmaskmap = get(vertex_orientmask,g);
  typename property_map < TGraph, vertex_vlabelmask_t >::type labelmaskmap = get(vertex_vlabelmask,g);
  typename property_map < TGraph, vertex_orient_t>::type orientmap = get(vertex_orient,g);
  typename property_map < TGraph, edge_labelmask_t >::type elabelmaskmap = get(edge_labelmask,g);
  typename property_map < TGraph, vertex_label_t >::type labelmap = get(vertex_label,g);
  typename property_map < TGraph, vertex_len_t >::type lenmap = get(vertex_len,g);
  typename property_map < TGraph, vertex_genome_t >::type genomemap = get(vertex_genome,g);

  //Variables
  //

  Edge e1;
  bool found;
  unsigned int numnewedges=0;

  //sort by coordinates by position on each sequence in the graph

#ifdef DEBUG
  std::cerr << "Updating adjacency edges in alignment graph" << std::endl;
#endif
  //TODO looping over all seqs is a bottleneck for draft genomes
  //Refactor by looping over graph and saving map [seqidx]->[vertex set]
  std::map<Label,std::vector<typename TGraph::vertex_descriptor> > seqVertexMap;
  for(typename boost::graph_traits<TGraph>::vertex_iterator 
	vit = vertices(g).first;vit!=vertices(g).second;++vit){
    //assert(labelmap.find(*vit)!=labelmap.end());
    assert(labelmap[*vit].size() > 0);
    for(LabelSet::iterator sit = labelmap[*vit].begin();sit!=labelmap[*vit].end();++sit){
      seqVertexMap[*sit].push_back(*vit);
    }
  }
  
  std::set<Label> skipseqs;
  for(typename std::map<Label,std::vector<typename TGraph::vertex_descriptor> >::iterator mit = seqVertexMap.begin();mit!=seqVertexMap.end();++mit){
    unsigned int spanlen=0;
    Label seqidx = mit->first;
    for(typename std::vector<typename TGraph::vertex_descriptor>::iterator vit=mit->second.begin();vit!=mit->second.end();++vit){
      typename TGraph::vertex_descriptor v = *vit;
      assert(coordinates.find(std::make_pair(v,seqidx))!=coordinates.end());
      if(coordinates.find(std::make_pair(v,seqidx))!=coordinates.end()){
	spanlen = spanlen + get(vertex_len,g,v);
      }
      else{
	assert(false);
      }
    }
    if(spanlen==0){
      skipseqs.insert(mit->first);
    }
  }

  for(typename std::map<Label,std::vector<typename TGraph::vertex_descriptor> >::iterator mit = seqVertexMap.begin();mit!=seqVertexMap.end();++mit){
    Label seqidx = mit->first;
    assert(sequence2genome.find(seqidx)!=sequence2genome.end());
    Label genomeidx = sequence2genome[seqidx];
    if(skipseqs.find(seqidx)==skipseqs.end()){
      //sort(sortedV.begin(),sortedV.end(),coordsorder(&coordinates,seqidx));
      sort(mit->second.begin(),mit->second.end(),coordsorder(&coordinates,seqidx));
      
      //
      //(8.1)Check and add any new edges between adjacent alignment blocks in an LCB
      for(std::vector<Vertex>::iterator it2 = mit->second.begin();it2!=mit->second.end();++it2){
	if(it2+1!=mit->second.end() && ccvmap[*it2]!=ccvmap[*(it2+1)]){//only consider new edges that bridge clusters
	  //check still on same sequence,genome
	  assert(labelmap[*it2].find(seqidx) != labelmap[*it2].end());
	  assert(genomemap[*it2].find(genomeidx) != genomemap[*it2].end());
	  assert(labelmap[*(it2+1)].find(seqidx) != labelmap[*(it2+1)].end());
	  assert(genomemap[*(it2+1)].find(genomeidx) != genomemap[*(it2+1)].end());
	  assert(coordinates.find(std::make_pair(*(it2+1),seqidx)) != coordinates.end());
	  assert(coordinates.find(std::make_pair(*it2,seqidx)) != coordinates.end());
	  //new edge exists only if dist < distance threshold
	  //int dist = abs(coordinates[std::make_pair(*it2,seqidx)].second - coordinates[std::make_pair(*(it2+1),seqidx)].first);
	  int dist = coordinates[std::make_pair(*(it2+1),seqidx)].first - coordinates[std::make_pair(*it2,seqidx)].second;
#ifdef NDEBUG
	  
#else
	  BitMask sharedlabels = (labelmaskmap[*it2]&labelmaskmap[*(it2+1)]);
#endif
	  assert(isLabelCollinearMask(sharedlabels,
				      orientmaskmap[*it2],
				      orientmaskmap[*(it2+1)]) 
		 == 
		 isLabelCollinear(orientmap[*it2],
				  orientmap[*(it2+1)],
				  sequence2genome));
	  //Additional checks to ensure that we only add "good" edges, between vertices on the same genomes
	  if(dist <= (int)distance 
	     && isLabelCollinear(orientmap[*it2],
				 orientmap[*(it2+1)],
				 sequence2genome) 
	     && isLabelMaxGap(*it2,*(it2+1),orientmap[*it2],orientmap[*(it2+1)],coordinates,distance,sequence2genome)){
	    //make sure that we do not introduce a rearrangment
	    //make sure that we do not introduce a long gap
	    
	    LCB newlcb;	
	    newlcb.insert(newlcb.end(),componentMap[ccvmap[*it2]].begin(),componentMap[ccvmap[*it2]].end());
	    newlcb.insert(newlcb.end(),componentMap[ccvmap[*(it2+1)]].begin(),componentMap[ccvmap[*(it2+1)]].end());
	    
	    //
	    BitMask longlabelmask=setSpanMask(newlcb,lenmap,labelmap,sequence2genome);

	    //Two check required
	    //Check if orientation of lcb1 and lcb2 are congruent
	    //Check if orientation of it2 and it2+1 are congruent	  
	    //TODO
	    //The checkLCBOrient(masks) will not consider the case where a single vertex 
	    //can be flipped to match the orientation
	    //TODO
	    //Only checking overall lcb orientation currently
	    //checkPairOrient(vlabelmap,vorientmap,*it1,*it2])
	    if(checkLCBOrient(g,newlcb,longlabelmask,sequence2genome)
	       && checkLCBOrient(lcborientmap,ccvmap[*it2],ccvmap[*(it2+1)],longlabelmask)
	       && checkLCBGaps(g,newlcb,ccvmap,coordinates,distance,sequence2genome)){
	      tie(e1,found) = edge(*it2,*(it2+1), g);
#ifdef DEBUG
		std::cerr << "Adding new edge between " << get(vertex_name,g,*it2) << "-" 
			  << get(vertex_name,g,*(it2+1)) << std::endl;
#endif
	      numnewedges++;
	      if(found){
		//TODO
		//addEdgeLabel(g,e1,genomeidx);
		if(!elabelmaskmap[e1].test(genomeidx)){
#if defined(STORE_EDGE_LABELS)
		  labelmap[e1].insert(std::make_pair(genomeidx,dist));
#endif
		  assert(!elabelmaskmap[e1].test(genomeidx));
		}
		else{
		  numnewedges--;
		}
		elabelmaskmap[e1].set(genomeidx,1);
	      }
	      else{
		//TODO
		//addEdgeLabel(g,e1,genomeidx)
		tie(e1,found) = edge(*(it2+1),*it2, g);
		if(found){
		  if(!elabelmaskmap[e1].test(genomeidx)){
#if defined(STORE_EDGE_LABELS)
		    labelmap[e1].insert(std::make_pair(genomeidx,dist));
#endif
		    assert(!elabelmaskmap[e1].test(genomeidx));
		  }
		  else{
		    numnewedges--;
		  }
		  elabelmaskmap[e1].set(genomeidx,1);
		}
		else{
#if defined(STORE_EDGE_LABELS)
		  LabelMap plabels;
		  plabels[genomeidx] = dist;
		  tie(e1,found) = add_edge(*it2,*(it2+1),EdgeProperties(plabels),baseg);
#else
		  tie(e1,found) = add_edge(*it2,*(it2+1),EdgeProperties(),baseg);
#endif
		  //TODO
		  //addEdgeLabel(g,e1,genomeidx)
		  elabelmaskmap[e1].set(genomeidx,1);
		  /*
		  //Remove any mask on this edges, only necessary if g is a filtered graph
		  if(maskedEdges.find(std::make_pair(*it2,*(it2+1)))!=maskedEdges.end()){
		  maskedEdges.erase(maskedEdges.find(std::make_pair(*it2,*(it2+1))));
		  }
		  if(maskedEdges.find(std::make_pair(*(it2+1),*it2))!=maskedEdges.end()){
		  maskedEdges.erase(maskedEdges.find(std::make_pair(*(it2+1),*it2)));
		  }
		  */
		}
	      }
	    }
	  }
	}
      }
    }
    else{
      //std::cerr << "Skipping merge on seq:" << seqidx 
      //<< " spanlen:" << spanlen << std::endl;
    }
  }
  //std::cerr << "Added " << numnewedges << " edges" << std::endl;
  setedgemasks(g,distance,coordinates,sequence2genome);
  setvertexmasks(g,sequence2genome);
  //std::cerr << "Finished setting edge and vertex masks" << std::endl;

}
//
//Merge adjacent lcbs
//This sub will only merge entire LCBs that are congruent
//TODO, consider edges from longest LCBs first
// Populate edgelcbmap with max(lcblenmap[source(e)],lcblenmap[target(e)])
// lcblenmap[edgelcbmap[e1]] < lcblenmap[edgelcbmap[e2]] 
template<typename TGraph, typename TLCBMap> 
int mergeLCBsGreedy(TGraph & g,
		    std::vector<int> & ccvmap,
		    std::vector<LCB> & componentMap,
		    TLCBMap & lcborientmap,
		    VertexLabelIntervalMap & coordinates,
		    EdgeSet & maskedEdges,
		    unsigned int maxgap,
		    SequenceGenomeMap & sequence2genome){
#ifdef DEBUG
  std::cerr << "Merging LCBs. Total count " << componentMap.size() << std::endl;
#endif
  typename property_map < TGraph, vertex_label_t >::type labelmap = get(vertex_label,g);
  typename property_map < TGraph, vertex_len_t >::type lenmap = get(vertex_len,g);
  typename property_map < TGraph, vertex_orient_t >::type omap = get(vertex_orient, g);
  typename property_map < TGraph, vertex_len_t >::type lmap = get(vertex_len, g);
  typename graph_traits<TGraph>::out_edge_iterator ei, edge_end;
  typename graph_traits<TGraph>::in_edge_iterator ei2, edge_end2;

  LCBLabelIntervalMap lcbcoords;  
  
  std::vector<int> ccremap=ccvmap;
  std::set<std::pair<int,int> > searches;
  int lcbcount=componentMap.size();
  int nummerges=0;

  //Capture LCB length
  std::map<int,int> lcblenMap; //lcbid->max_seq_span
  std::vector<int> lcbidx;
  for(unsigned int k=0;k<componentMap.size();++k){
#ifdef DEBUG
    std::cerr << "Component " << k << std::endl;
#endif
    if(componentMap[k].size()>0){
      int bplen=0;
      unsigned int len = get_LCB_length(componentMap[k],omap,lmap,coordinates,lcbcoords,k,bplen,sequence2genome,0); 
#ifdef DEBUG
      std::cerr << " len:" << len << std::endl;
#endif
      lcblenMap[k] = len;
    }
    else{
      lcblenMap[k]=0;
    }
    lcbidx.push_back(k);
  }
  
  //Sort LCBs on length
  sort(lcbidx.begin(),lcbidx.end(),lencmp(lcblenMap));
  //Greedy merge adjacent LCBs from largest to smallest
  for(std::vector<int>::reverse_iterator cit = lcbidx.rbegin();cit != lcbidx.rend();++cit){
#ifdef DEBUG
    std::cerr << "Greedy merge LCB:" << *cit << " len:" << lcblenMap[*cit] << std::endl;
#endif
    if(componentMap[*cit].size()>0){
      std::vector<Vertex> lcbv = componentMap[*cit];
      for(LCB::iterator vit = lcbv.begin();vit!=lcbv.end();++vit){
	std::vector<Edge> lcbedges;
	tie(ei,edge_end) = out_edges(*vit,g);
	for(;ei!=edge_end;++ei){
	  assert(source(*ei,g)==*vit);
	  lcbedges.push_back(*ei);
	}
	
	tie(ei2,edge_end2) = in_edges(*vit,g);
	for(;ei2!=edge_end2;++ei2){
	  assert(target(*ei2,g)==*vit);
	  lcbedges.push_back(*ei2);
	}
#ifdef DEBUG
	std::cerr << "Edges " << lcbedges.size() << std::endl;
#endif
	for(vector<Edge>::iterator eit=lcbedges.begin();eit!=lcbedges.end();++eit){
	  Vertex sv = source(*eit,g);
	  Vertex tv = target(*eit,g);
#ifdef DEBUG
	  std::cerr << "Vertex " << sv << "-" << tv << std::endl;
#endif
	  int sidx = ccvmap[sv];
	  int tidx = ccvmap[tv];
	  //
	  //If edge connects two components consider merging
	  //if compatible
	  if(sidx!=tidx
	     && searches.find(std::make_pair(sidx,tidx)) == searches.end()
	     && searches.find(std::make_pair(tidx,sidx)) == searches.end()
	     && componentMap[tidx].size()> 0
	     && componentMap[sidx].size()> 0
	     ){
	    //Make sure that there is no edge already connecting these LCBs
	    //This merge must be run after computing connected components
	    //
	    assert(maskedEdges.find(std::make_pair(sv,tv)) != maskedEdges.end()
		   ||maskedEdges.find(std::make_pair(tv,sv)) != maskedEdges.end());
	    //
	    //Mark that CC pair has been searched
	    searches.insert(std::make_pair(sidx,tidx));
	    
	    LCB newlcb;
	    newlcb.insert(newlcb.end(),componentMap[sidx].begin(),componentMap[sidx].end());
	    newlcb.insert(newlcb.end(),componentMap[tidx].begin(),componentMap[tidx].end());
	    
	    BitMask longlabelmask=setSpanMask(newlcb,lenmap,labelmap,sequence2genome);
	    
	    //TODO
	    //The funcs called in this loop are a bottleneck according to gprof 
	    //Most time spent copying OrientedLabelSet
	    //First this is a large loop, all edges. 
	    //checkLCBgaps/checkLCBOrient makes copies of vertex properties like orientedlabelset
	    //checkLCBGaps creates and sorts vectors
	    if(checkLCBOrient(g,newlcb,longlabelmask,sequence2genome) //pairwise check for all vertices
	       && checkLCBOrient(lcborientmap,sidx,tidx,longlabelmask) //check consistency with lcb orient
	       && checkLCBGaps(g,newlcb,ccvmap,coordinates,maxgap,sequence2genome)){
	      //Save LCB
	      //std::cerr << "New LCB: "; 
	      
	      for(LCB::iterator vit2=newlcb.begin();vit2!=newlcb.end();++vit2){
		assert(*vit2<ccvmap.size());
		ccvmap[*vit2]=lcbcount;
	      }
	      //std::cerr << std::endl;
	      componentMap.push_back(newlcb);
	      //Clear out old LCB
	      componentMap[sidx] = LCB();
	      componentMap[tidx] = LCB();
	      
	      std::vector<Vertex> badV;
	      //If joined LCBs have same labels, orietation, no need to recalc
	      assert(lcborientmap.find(sidx)!=lcborientmap.end());
	      assert(lcborientmap.find(tidx)!=lcborientmap.end());
	      if(lcborientmap[sidx]==lcborientmap[tidx]){
		lcborientmap[lcbcount]=lcborientmap[tidx];
	      }
	      else{
		//Set the label and orientation for the new lcb
#ifdef DEBUG
		std::cerr << lcborientmap[sidx].first << std::endl;
		std::cerr << lcborientmap[sidx].second << std::endl << std::endl;
		std::cerr << lcborientmap[tidx].first << std::endl;
		std::cerr << lcborientmap[tidx].second << std::endl << std::endl;
#endif
		lcborientmap[lcbcount]=setLCBOrient(g,newlcb,badV,sequence2genome);
	      }
	      //Remove mask on edge linking two lcbs
	      EdgeSet::iterator mit = maskedEdges.find(std::make_pair(sv,tv));
	      if(mit != maskedEdges.end()){
		maskedEdges.erase(mit);
		Edge e1;
		bool found;
		tie(e1,found) = edge(mit->first,mit->second,g);
		assert(found);
		put(edge_category,g,e1,ORANGERED);
	      }
	      else{
		mit = maskedEdges.find(std::make_pair(tv,sv));
		assert(mit != maskedEdges.end());
		maskedEdges.erase(mit);
		Edge e1;
		bool found;
		tie(e1,found) = edge(mit->first,mit->second,g);
		assert(found);
		put(edge_category,g,e1,ORANGERED);
	      }
#ifdef DEBUG
	      std::cerr << "Merging LCB:"<<sidx<< " with LCB:"<<tidx<< " into LCB:"<<lcbcount << std::endl;
#endif
	      nummerges++;
	      lcbcount++;
	      assert(lcbcount==(int)componentMap.size());
	    }
	    else{
	      //skip LCB
#ifdef DEBUG
	      std::cerr << "Skipping merge of LCB:"<<sidx
			<< " with LCB:"<<tidx
			<< " from edge " 
			<< get(vertex_name,g,sv) << "-" << get(vertex_name,g,tv) 
			<< std::endl;
#endif
	    }
	  }
	}
      }
    }
  }
  return nummerges;
}
//
//Merge adjacent lcbs
//This sub will only merge entire LCBs that are congruent
//TODO, consider edges from longest LCBs first
// Populate edgelcbmap with max(lcblenmap[source(e)],lcblenmap[target(e)])
// lcblenmap[edgelcbmap[e1]] < lcblenmap[edgelcbmap[e2]] 
template<typename TGraph, typename TLCBMap> 
int mergeLCBs(TGraph & g,
	      std::vector<int> & ccvmap,
	      std::vector<LCB> & componentMap,
	      TLCBMap & lcborientmap,
	      VertexLabelIntervalMap & coordinates,
	      EdgeSet & maskedEdges,
	      unsigned int maxgap,
	      SequenceGenomeMap & sequence2genome){
#ifdef DEBUG
  std::cerr << "Merging LCBs. Total count " << componentMap.size() << std::endl;
#endif
  typename property_map < TGraph, vertex_label_t >::type labelmap = get(vertex_label,g);
  typename property_map < TGraph, vertex_len_t >::type lenmap = get(vertex_len,g);
  
  std::vector<int> ccremap=ccvmap;
  std::set<std::pair<int,int> > searches;
  int lcbcount=componentMap.size();
  int nummerges=0;
  typename boost::graph_traits<TGraph>::edge_iterator eit, edge_end;
  edge_end=edges(g).second;
  for(eit=edges(g).first;eit!=edge_end;++eit){//all edges in g
    Vertex sv = source(*eit,g);
    Vertex tv = target(*eit,g);
    int sidx = ccvmap[sv];
    int tidx = ccvmap[tv];
    //
    //If edge connects two components consider merging
    //if compatible
    if(sidx!=tidx
       && searches.find(std::make_pair(sidx,tidx)) == searches.end()
       && searches.find(std::make_pair(tidx,sidx)) == searches.end()
       ){
      //Make sure that there is no edge already connecting these LCBs
      //This merge must be run after computing connected components
      //
      assert(maskedEdges.find(std::make_pair(sv,tv)) != maskedEdges.end()
	     ||maskedEdges.find(std::make_pair(tv,sv)) != maskedEdges.end());
      //
      //Mark that CC pair has been searched
      searches.insert(std::make_pair(sidx,tidx));

      LCB newlcb;
      newlcb.insert(newlcb.end(),componentMap[sidx].begin(),componentMap[sidx].end());
      newlcb.insert(newlcb.end(),componentMap[tidx].begin(),componentMap[tidx].end());
      
      BitMask longlabelmask=setSpanMask(newlcb,lenmap,labelmap,sequence2genome);
      
      //TODO
      //The funcs called in this loop are a bottleneck according to gprof 
      //Most time spent copying OrientedLabelSet
      //First this is a large loop, all edges. 
      //checkLCBgaps/checkLCBOrient makes copies of vertex properties like orientedlabelset
      //checkLCBGaps creates and sorts vectors
      if(checkLCBOrient(g,newlcb,longlabelmask,sequence2genome) //pairwise check for all vertices
	 && checkLCBOrient(lcborientmap,sidx,tidx,longlabelmask) //check consistency with lcb orient
	 && checkLCBGaps(g,newlcb,ccvmap,coordinates,maxgap,sequence2genome)){
	//Save LCB
	//std::cerr << "New LCB: "; 
	for(LCB::iterator vit=newlcb.begin();vit!=newlcb.end();++vit){
	  assert(*vit<ccvmap.size());
	  ccvmap[*vit]=lcbcount;
	  //std::cerr << get(vertex_name,g,*vit) << " ";
	}
	//std::cerr << std::endl;
	componentMap.push_back(newlcb);
	std::vector<Vertex> badV;
	//If joined LCBs have same labels, orietation, no need to recalc
	assert(lcborientmap.find(sidx)!=lcborientmap.end());
	assert(lcborientmap.find(tidx)!=lcborientmap.end());
	if(lcborientmap[sidx]==lcborientmap[tidx]){
	  lcborientmap[lcbcount]=lcborientmap[tidx];
	}
	else{
	  //Set the label and orientation for the new lcb
#ifdef DEBUG
	  std::cerr << lcborientmap[sidx].first << std::endl;
	  std::cerr << lcborientmap[sidx].second << std::endl << std::endl;
	  std::cerr << lcborientmap[tidx].first << std::endl;
	  std::cerr << lcborientmap[tidx].second << std::endl << std::endl;
#endif
	  lcborientmap[lcbcount]=setLCBOrient(g,newlcb,badV,sequence2genome);
	}
	//Remove mask on edge linking two lcbs
	EdgeSet::iterator mit = maskedEdges.find(std::make_pair(sv,tv));
	if(mit != maskedEdges.end()){
	  maskedEdges.erase(mit);
	  Edge e1;
	  bool found;
	  tie(e1,found) = edge(mit->first,mit->second,g);
	  assert(found);
	  put(edge_category,g,e1,ORANGERED);
	}
	else{
	  mit = maskedEdges.find(std::make_pair(tv,sv));
	  assert(mit != maskedEdges.end());
	  maskedEdges.erase(mit);
	  Edge e1;
	  bool found;
	  tie(e1,found) = edge(mit->first,mit->second,g);
	  assert(found);
	  put(edge_category,g,e1,ORANGERED);
	}
#ifdef DEBUG
	    std::cerr << "Merging LCB:"<<sidx<< " with LCB:"<<tidx<< " into LCB:"<<lcbcount << std::endl;
#endif
	    nummerges++;
	    lcbcount++;
	    assert(lcbcount==(int)componentMap.size());
      }
      else{
	//skip LCB
#ifdef DEBUG
	std::cerr << "Skipping merge of LCB:"<<sidx
		  << " with LCB:"<<tidx
		  << " from edge " 
		  << get(vertex_name,g,sv) << "-" << get(vertex_name,g,tv) 
		  << std::endl;
#endif
      }
    }
  }
  return nummerges;
}

//
//Completely remove the LCB from the graph (by adding to maskedLCBs)
void removeLCB(LCB & lcb, 
	       std::set<std::pair<Vertex,bool> > &breakpoints, 
	       VertexSet &maskedLCBs){
  std::set<std::pair<Vertex,bool> >::iterator it2;
  for(LCB::iterator vit = lcb.begin();vit!=lcb.end();++vit){
#ifdef DEBUG
      std::cerr << "Removing vertex " << *vit << std::endl;
#endif
    maskedLCBs.insert(*vit);
  }
}


//
//Mark possible syntenic breakpoints in graph g, storing in maskedEdges
//Breakpoints can arise from 
//(1)Change in label
//(2)Change in orientation
//(3)Flux, whereever indegree!=1 or outdegree!=1
template<typename TGraph, typename BPMap1, typename BPMap2, typename VMap1>
void markBreakpoints(TGraph &g, 
		     BPMap1 &breakpoints, 
		     BPMap2 &maskedEdges, 
		     VMap1 &vertexList, 
		     SequenceGenomeMap &sequence2genome){
  typename graph_traits<TGraph>::vertex_iterator i, end;
  typename graph_traits<TGraph>::out_edge_iterator ei, edge_end;
  typename graph_traits<TGraph>::in_edge_iterator ei2, edge_end2;
  typename property_map < TGraph, vertex_orientmask_t >::type vorientmap = get(vertex_orientmask, g);
  typename property_map < TGraph, vertex_vlabelmask_t >::type vlabelmap = get(vertex_vlabelmask, g);
  typename property_map < TGraph, vertex_orient_t >::type vmap = get(vertex_orient, g);
  typename property_map < TGraph, edge_labelmask_t >::type elabelmap = get(edge_labelmask, g);
  int bptype1=0;
  int bptype2=0;
  int bptype3=0;
  int keepmerge=false;

  for(typename boost::graph_traits<TGraph>::vertex_iterator 
	vit = vertices(g).first;vit!=vertices(g).second;++vit){
      Vertex v = *vit;
      if(vertexList.size()>0 && vertexList.find(v)==vertexList.end())
	continue;
#ifdef DEBUG
	std::cerr << "Checking for breakpoints on vertex v:" << get(vertex_name,g,v) << std::endl;
#endif
	//ei = out_edges(v, g).first;
	// bptype==0 no breakpoint
	// bptype==1 incoming bp, end a region
	// bptype==2 outgoing bp, start a region
	bool inlinebp=false;
	bool fluxbp=false;
	bool inbp=false;
	bool outbp=false;
	bool ismerge=false;

      if(in_degree(v,g)==1){
	tie(ei2,edge_end2) = in_edges(v,g);
	assert(target(*ei2,g)==v);
	//Check same labels
	if(sameLabel(vlabelmap[v],vlabelmap[source(*ei2,g)],elabelmap[*ei2])){	
	}
	else{
	  //Some type of case (3) flux
	  if(isLabelCollinear(vmap[v],vmap[source(*ei2,g)],sequence2genome)){
	  }
	  else{
	    fluxbp=true;
	    inbp=true;
	    maskedEdges.insert(std::make_pair(source(*ei2,g),v));
	    put(edge_category,g,*ei2,GREEN);
	    bptype1++;
	  }
	}
	if(isLabelCollinear(vmap[v],vmap[source(*ei2,g)],sequence2genome)){
	  assert(sameOrient(vorientmap[v]&elabelmap[*ei2],vorientmap[source(*ei2,g)]&elabelmap[*ei2],vlabelmap[v]&elabelmap[*ei2]));
	}
	else{
	  //Some type of case (2) orientation change
	  fluxbp=true;
	  inlinebp=false;
	  inbp=true;
	  maskedEdges.insert(std::make_pair(source(*ei2,g),v));
	  put(edge_category,g,*ei2,PURPLE);
	  bptype2++;
	}
      }
      if(out_degree(v,g)==1){
	tie(ei,edge_end) = out_edges(v,g);
	assert(source(*ei,g)==v);
	if(sameLabel(vlabelmap[v],vlabelmap[target(*ei,g)],elabelmap[*ei])){
	}
	else{
	  //Some type of case (3) flux
	  if(isLabelCollinear(vmap[v],vmap[target(*ei,g)],sequence2genome)){
	  }
	  else{
	    fluxbp=true;
	    //inlinebp=true;
	    outbp=true;
	    maskedEdges.insert(std::make_pair(v,target(*ei,g)));
	    put(edge_category,g,*ei,GREEN);
	    bptype1++;
	  }
	}
	if(isLabelCollinear(vmap[v],vmap[target(*ei,g)],sequence2genome)){
	  assert(sameOrient(vorientmap[v]&elabelmap[*ei],
			    vorientmap[target(*ei,g)]&elabelmap[*ei],
			    vlabelmap[v]&elabelmap[*ei]));
	}
	else{
	  //Some type of case (2) orientation change
	  maskedEdges.insert(std::make_pair(v,target(*ei,g)));
	  put(edge_category,g,*ei,PURPLE);
	  fluxbp=true;
	  inlinebp=false;
	  outbp=true;
	  bptype2++;
	}
      }
      if(in_degree(v,g)>1){
	//Some type of case (3) flux
	fluxbp=true;
	inbp=true;
	tie(ei2,edge_end2) = in_edges(v,g);

	for(;ei2!=edge_end2;++ei2){
	  assert(target(*ei2,g)==v);
	  //maskedEdges.insert(std::make_pair(source(*ei2,g),target(*ei2,g)));
#ifdef DEBUG
	  std::cerr << "Adding bp " << get(vertex_name,g,source(*ei2,g)) << "-" << get(vertex_name,g,target(*ei2,g)) << std::endl;
#endif
	  if(isLabelCollinear(vmap[source(*ei2,g)],vmap[target(*ei2,g)],sequence2genome)){
	    //Previously merged, keep
	    if(keepmerge && get(edge_category,g,*ei2)==ORANGERED){
	      ismerge=true;
	    }
	    else{
	      put(edge_category,g,*ei2,RED);
	      maskedEdges.insert(std::make_pair(source(*ei2,g),target(*ei2,g)));
	    }
	  }
	  else{
	    //Some type of case (2) orientation change
	    put(edge_category,g,*ei2,PURPLE);
	    maskedEdges.insert(std::make_pair(source(*ei2,g),target(*ei2,g)));
	  }
	}
	bptype3++;
      }
      if(out_degree(v,g)>1){
	//Some type of case (3) flux
	fluxbp=true;
	outbp=true;
	tie(ei,edge_end) = out_edges(v,g);
	for(;ei!=edge_end;++ei){
	  assert(source(*ei,g)==v);
	  //maskedEdges.insert(std::make_pair(source(*ei,g),target(*ei,g)));
#ifdef DEBUG
	  std::cerr << "Adding bp " << get(vertex_name,g,source(*ei,g)) << "-" << get(vertex_name,g,target(*ei,g)) << std::endl;
#endif
	  if(isLabelCollinear(vmap[source(*ei,g)],vmap[target(*ei,g)],sequence2genome)){
	    if(keepmerge && get(edge_category,g,*ei)==ORANGERED){
	      ismerge=true;
	    }
	    else{
	      put(edge_category,g,*ei,RED);
	      maskedEdges.insert(std::make_pair(source(*ei,g),target(*ei,g)));
	    }
	  }
	  else{
	    //Some type of case (2) orientation change
	    put(edge_category,g,*ei,PURPLE);
	    maskedEdges.insert(std::make_pair(source(*ei,g),target(*ei,g)));
	  }
	}
	bptype3++;
      }
      if(out_degree(v,g)==0){
	outbp=true;
      }
      if(in_degree(v,g)==0){
	inbp=true;
      }
      if(fluxbp){
	breakpoints.insert(std::make_pair(v,false));
      }      
    } 
#ifdef DEBUG
    std::cerr << "Marked breakpoints. type1 " << bptype1 << " type2 " << bptype2 << " type3 " << bptype3 << std::endl;
#endif
}

//
//Remove all breakpoints except
//PURPLE orientation changing
//BLUE mincuts
template<typename TGraph>
void clearInlineBreakpoints(TGraph & g,
			    EdgeSet  &maskedEdges){
  vector<EdgeSet::iterator > eraseMask;
  EdgeSet::iterator mit;
  for(mit = maskedEdges.begin();mit!=maskedEdges.end();++mit){
    Edge e;
    bool found;
    tie(e,found) = edge(mit->first,mit->second,g);
    assert(found);
    if(get(edge_category,g,e)!=PURPLE && get(edge_category,g,e)!=BLUE){
      eraseMask.push_back(mit);
    }
    tie(e,found) = edge(mit->second,mit->first,g);
    if(found){
      if(get(edge_category,g,e)!=PURPLE && get(edge_category,g,e)!=BLUE){
	eraseMask.push_back(mit);
      }
    }
  }
  vector<EdgeSet::iterator >::iterator eit;
  for(eit=eraseMask.begin();eit!=eraseMask.end();++eit){
    maskedEdges.erase(*eit);
  }
}


//############
//Connected components
template<typename TGraph, typename TGraphBase, typename TComponentMap, typename TVertexMap, typename TLCBMap>
inline
int calc_components_undirected(TGraph & fg, 
			       TGraphBase & g,
			       TComponentMap & componentMap, 
			       TVertexMap & c, 
			       TLCBMap & lcborientmap,
			       SequenceGenomeMap & sequence2genome){

  typedef adjacency_list<vecS,vecS,undirectedS,VertexProperties,EdgeProperties> TLGraph;
  typedef typename TLGraph::vertex_descriptor TLVertex;
  typedef typename TLGraph::edge_descriptor TLEdge;

  typedef typename TGraph::vertex_descriptor TVertex;
  typedef typename TGraphBase::edge_descriptor TEdgeBase;

  typedef typename boost::graph_traits<TGraph>::edge_iterator TEdgeIterator;
  typedef typename boost::graph_traits<TGraph>::vertex_iterator TVertexIterator;
  typedef typename boost::graph_traits<TLGraph>::vertex_iterator TLVertexIterator;

  typedef typename boost::graph_traits<TGraphBase>::edge_iterator TEdgeBaseIterator;
  typedef typename boost::graph_traits<TGraphBase>::vertex_iterator TVertexBaseIterator;

  bool inserted;
  //Undirected graph(currlcbg) is required here for the CC algorithm
  //TODO  Performance enhancement refactor
  //      Avoid building a second graph and run CC on directed graph(fg)
  adjacency_list<vecS,vecS,undirectedS,VertexProperties,EdgeProperties> currlcbg;
  std::map<VertexName, TLVertex> currlcbv;
  std::map<VertexName, TLVertex>::iterator pos;
  //Map between undirected graph(currlcbg) vertices and directed graph(fg) vertices
  std::map<TLVertex,TVertex> vmap;
  TLVertex news,newt;
  TEdgeIterator starte,ende;
  TVertexIterator startv,endv;
  TLEdge ne;

  tie(startv,endv)=vertices(fg);
  for(TVertexIterator vit = startv;vit!=endv;++vit){
    VertexName sname = get(vertex_name,fg,*vit);
    assert(currlcbv.find(sname)==currlcbv.end());
    tie(pos, inserted) = currlcbv.insert(std::make_pair(sname, TLVertex()));
    assert(inserted);
    news = add_vertex(sname,currlcbg);
    currlcbv[sname]=news;
    assert(vmap.find(news)==vmap.end());
    vmap[news]=*vit;
  }

  tie(starte,ende)=edges(fg);
  for(TEdgeIterator eit = starte;eit!=ende;++eit){
    TEdgeBase e = *eit;
    VertexName tname = get(vertex_name,fg,target(e,fg));
    VertexName sname = get(vertex_name,fg,source(e,fg));
    assert(currlcbv.find(sname)!=currlcbv.end());
    assert(currlcbv.find(tname)!=currlcbv.end());
    news=currlcbv[sname];
    newt=currlcbv[tname];
    tie(ne, inserted) = add_edge(news,newt,currlcbg);
    if(inserted){
      
    }
  }

  assert(currlcbv.size()==num_vertices(currlcbg));
  c.clear();
  c.resize(num_vertices(currlcbg));  
  int numComponents = connected_components(currlcbg,&c[0]);
  
  //Save mapping of componentNum->vector<Vertex>
  componentMap.clear();
  assert(componentMap.size()==0);
  componentMap.resize(numComponents);
  for(TLVertexIterator vit = vertices(currlcbg).first;vit!=vertices(currlcbg).second;++vit){ 
    assert(vmap.find(*vit)!=vmap.end());
    //This ensures lcbidx=c[vertex]
    componentMap[c[*vit]].push_back(vmap[*vit]);
  }  
  //
  //Save mask for the LCB
  setLCBOrient(g,lcborientmap,componentMap,sequence2genome);
  return numComponents;
}




//
//Calculate some summary statistics
//numLCBs,minlcblen(bp),totallcblen(bp),avglen(bp)
template<typename TGraph>
void summaryStats(TGraph & fglcbsyn,
		  std::vector<LCB> &componentMap, 
		  VertexLabelIntervalMap &coordinates, 
		  unsigned int minlength,
		  int & numc,
		  unsigned int & minlen,
		  int & totallen,
		  unsigned int & avglen,
		  int & maxv,
		  SequenceGenomeMap & sequence2genome){
  
  LCBLabelIntervalMap lcbcoords;  
  avglen=0;
  minlen=std::numeric_limits<unsigned int>::max();
  totallen=0;
  numc=0;
  maxv=0;
  for(unsigned int k=0;k<componentMap.size();++k){
    if(componentMap[k].size()>0){
      maxv = (maxv > (int)componentMap[k].size()) ? maxv : componentMap[k].size();
      assert(componentMap[k].size()>0);
      OrientedLabelSet label = get(vertex_orient, fglcbsyn, componentMap[k][0]);
      property_map < LCBSynFilterGraph, vertex_orient_t >::type omap = get(vertex_orient, fglcbsyn);
      property_map < LCBSynFilterGraph, vertex_len_t >::type lmap = get(vertex_len, fglcbsyn);
      int bplen=0;
      unsigned int len = get_LCB_length(componentMap[k],omap,lmap,coordinates,lcbcoords,k,bplen,sequence2genome,minlength); 
      if(len>0){
	if(len>=minlength){
	  minlen = (len < minlen) ? len : minlen;
	  avglen+=len;
	  totallen+=bplen;
	  numc++;
	}
	else{
	}
      }
    }
  }
}



int main(int argc, char* argv[])
{

  //Number of iterations to run
  unsigned int MAXITERS=5;
  unsigned int MAXSTABLE=1;

  //Input graph
  Graph g;

  //Key parameters
  unsigned int distance=0; //maximum gap length between anchors
  unsigned int shortlcblen=0; //maximum length of LCBs that are masked during chaining

  unsigned int minlength=0; //for reporting stats only
  unsigned int minanchor=0; //minimum anchor length
  unsigned int minprintlength=0;

  //Ensure chains do not overlap by removing overlapping regions
  bool removeoverlaps=false;

  //Lookups
  NameVertexMap name2vertex,name2vertexcomp;  
  NameLabelMap sequence2index,genome2index;  
  LabelNameMap index2sequence;
  SequenceGenomeMap sequence2genome;
  
  //Map of coordinates for each anchor
  VertexLabelIntervalMap coordinates;
  VertexLabelIntervalMap::iterator cpos;

#ifdef TIMING
  time_t now;
  time(&now);
  time_t lasttime=now;
#endif

  if(argc<=3){
    cerr << "USAGE:mugsy-chaining max-distance min-lcbspan min-statslen < anchors.projection" << std::endl;
    exit(1);
  }

  if(argc>1){
    assert(atoi(argv[1])>=0);
    distance = atoi(argv[1]);
  }
  if(argc>2){
    assert(atoi(argv[2])>=0);
    shortlcblen = atoi(argv[2]);
  }
  if(argc>3){
    minlength = atoi(argv[3]);
  }
  assert(distance>0);
  assert(minlength>=0);
  cerr << "#Using custom distance " << distance << endl;
  cerr << "#Using custom minlength " << minlength << endl;
  
  cerr << "#Parsing graph from stdin" << endl;
  if(0){
    //TODO
    //Check file format
    //Allow for unprojected,projected list of blocks
    //Read blocks and build alignment graph perform projection over
    //each sequence and connect blocks that are adjacent on any given
    //sequence at distance < d
    read_blocks(std::cin,
		g,
		name2vertex,
		genome2index,
		sequence2index,
		coordinates,
		distance);
  }
  else{
    //Read a projection of anchors and build anchor graph
    //Only consider anchors that are adjacent < d
    read_pairwiseprojection(std::cin,
			    g,
			    name2vertex,
			    genome2index,
			    sequence2index,
			    coordinates,
			    sequence2genome,
			    distance,
			    minanchor);
    //Save coordinates for each anchor in coordinates map
    updateCoordinates(coordinates,sequence2genome); 
  }


  //Reverse sequence2index map
  for(NameLabelMap::iterator i = sequence2index.begin();i!=sequence2index.end();++i){
#ifdef DEBUG
    std::cerr << "Seq idx:" << i->second << " " << i->first << std::endl;
#endif
    index2sequence[i->second] = i->first;
  }

  //Restrict to a set of labels
  LabelSet labels;
  if(argc>4){
    for(int i=4;i<argc;++i){
      NameLabelMap::iterator it = sequence2index.find(argv[i]);
      if(it != sequence2index.end()){
	cerr << "#Restricting outputs to sequence label " << argv[i] << endl;
	labels.insert(it->second);
      }
      else{
	cerr << "#Invalid sequence label " << argv[i] << endl;
	assert(false);
      }
    }
  }
  
  cerr << "#Num of vertices " << num_vertices(g) << endl;
  cerr << "#Num of edges " << num_edges(g) << endl;

  //Set edge and vertex masks for fast pattern matching
  cerr << "#Setting edge and vertex masks" << endl;
  setedgemasks(g,distance,coordinates,sequence2genome);
  setvertexmasks(g,sequence2genome);

#ifdef TRIMEDGES
  //Remove edges connnected in only one label. This simplifies the
  //graph by removing flux contributed by a single genome only.
  std::vector<boost::graph_traits<Graph>::edge_descriptor> eraseEdges;
  for(boost::graph_traits<Graph>::edge_iterator 
	eit = edges(g).first;eit!=edges(g).second;++eit){
    Edge e = *eit;

    BitMask emask = get(edge_labelmask,g,*eit);
#if defined(STORE_EDGE_LABELS)
    LabelMap inlabels = get(edge_label,g,*eit);
    assert(inlabels.size()==emask.count());
#endif    
    if(emask.count()<=1){
      eraseEdges.push_back(*eit);
    }
  }

  for(std::vector<boost::graph_traits<Graph>::edge_descriptor>::iterator eit=eraseEdges.begin();eit!=eraseEdges.end();++eit){
    //std::cerr << "Removing edge:" << get(vertex_name,g,source(*eit,g)) << "-" << get(vertex_name,g,target(*eit,g)) << std::endl;
    remove_edge(*eit,g);
  }
#endif

  property_map < Graph, vertex_name_t >::type vertex_name_map = get(vertex_name, g);
  
  //Variables
  Edge e1;
  int itercount=0;
  //Extract all sequences
  std::set<Label> seqidxSet;
  for(NameLabelMap::iterator it = sequence2index.begin();it!=sequence2index.end();++it){
    //it->first is the sequence name
    //it->second is the index
    seqidxSet.insert(it->second);
  }
#ifdef TIMING
  time(&now);
  std::cerr << "TIME_INIT:" << now-lasttime << std::endl;
  lasttime=now;
#endif
  //Initial clustering: build anchor graph, cut, merge, maskshort, recluster
  //Building the anchor graph
  //(1) Create filtered graph that supports breakpoints and maskedLCBs
  std::cerr << "Building alignment graph and initial clutering" << std::endl;
  std::set<std::pair<Vertex,bool> > breakpoints;
  VertexSet maskedLCBs;

  //Store vertex pair rather than edge_descriptor to avoid problems with 
  //stale edge descriptors and lack of < operator needed for std::set
  EdgeSet maskedEdges; 

  LCBLabelIntervalMap lcbcoords;  

  //Filter graph has predicates for
  //-Masked edges
  //-Masked LCBs
  //
  //Edge filters
  synbp_edge_filter<Graph> synefilter(&maskedEdges,&g);
  LCB_edge_filter<Graph> lcbefilter(&maskedLCBs,&g);
  compound_edge_filter<LCB_edge_filter<Graph>, synbp_edge_filter<Graph> > 
    cmpefilter(lcbefilter,synefilter);
  //Vertex filters
  LCB_vertex_filter<Graph> lcbvfilter(&maskedLCBs);
  //The graph
  LCBSynFilterGraph fglcbsyn(g,cmpefilter,lcbvfilter);

  //
  //(2.1) Find and mark all breakpoints in the graph
  //Breakpoint types (stored in edge_category)
  //RED - potential syntenic brkpt due to multiple incoming/outgoing edges
  //PURPLE - change in orientation between sequences in adjacent blocks/vertices
  //GREEN - other flux such loss of homology in a single genome
  set<Vertex> dummySet;
  markBreakpoints(g,breakpoints,maskedEdges,dummySet,sequence2genome);
#ifdef DEBUG
  std::cerr << "Marked " << breakpoints.size() << " breakpoints" << std::endl;
#endif
  //
  //(2.2) Calculate LCBs using connected components
  //This initial clustering is expected to produce a over-segmented
  //set of LCBs. Later steps in the clustering will collapse LCBs
  std::vector<LCB> componentMap;
  std::vector<int> ccvmap(num_vertices(fglcbsyn));

#ifdef DEBUG
    do_write_graphviz(g, std::string("gout.input.dot"),ccvmap,coordinates,maskedEdges,maskedLCBs);
#endif

  //TODO
  //Replace with non BitMask version or use boost::dynamic_bitset
  //Lookup for lcbid->(labelmask,orientmask)
  std::map<int,std::pair<BitMask,BitMask> > lcborientmap;

  int numComponents = calc_components_undirected(fglcbsyn,g,componentMap,ccvmap,lcborientmap,sequence2genome);
  itercount=numComponents;
#ifdef TIMING
  time(&now);
  std::cerr << "TIME_CLUST1:" << now-lasttime << std::endl;
  lasttime=now;
#endif

  unsigned int avglen,minlen;
  int totallen,numc,maxv;
  int allbps;

#ifdef LCBSTATS
  //Calculate stats
  summaryStats(fglcbsyn,componentMap,coordinates,minprintlength, 
	       numc,minlen,totallen,avglen,maxv,sequence2genome); //using minprintlength
  allbps=totallen;
  if(numc>0){
    std::cerr << "LCB summary orig " << numc << " min:" << minlen << " coverage:" << totallen << "(" << (float)totallen/allbps << ")" << " avg_bp:" << avglen/numc << " maxv:" << maxv<< std::endl;
  }
  summaryStats(fglcbsyn,componentMap,coordinates,minlength, 
	       numc,minlen,totallen,avglen,maxv,sequence2genome); //using minlength
  if(numc>0){
    std::cerr << "LCB summary orig " << numc << " min:" << minlen << " coverage:" << totallen << "(" << (float)totallen/allbps << ")" << " avg_bp:" << avglen/numc << " maxv:" << maxv<< std::endl;
  }
#endif
#ifdef DEBUG
  do_write_graphviz(g, std::string("gout.orig.dot"),ccvmap,coordinates,maskedEdges,maskedLCBs);
  do_write_graphviz(fglcbsyn, std::string("gout.orig.dot.filtered"),ccvmap,coordinates,maskedEdges,maskedLCBs);
#endif
  
  std::cerr << "Partitioning graph to maintain contraints" << std::endl;
  //(2.3) Breaks LCBs based on gap lengths, mismatched orient, and mult seqs same genome
  int cutattempts=0;
  int origbreaks = breakLCBmincutconnect(componentMap,ccvmap,maskedEdges,g,fglcbsyn,distance,coordinates,seqidxSet,name2vertex,sequence2genome);
#ifdef DEBUG
  std::cerr << "Num orig breaks " << origbreaks << std::endl;
#endif
  numComponents = calc_components_undirected(fglcbsyn,g,componentMap,ccvmap,lcborientmap,sequence2genome);
#ifdef TIMING
  time(&now);
  std::cerr << "TIME_MINCUT1:" << now-lasttime << std::endl;
  lasttime=now;
#endif
#ifdef DEBUG
  do_write_graphviz(g, std::string("gout.mincut1.dot"),ccvmap,coordinates,maskedEdges,maskedLCBs);
  do_write_graphviz(fglcbsyn, std::string("gout.mincut1.dot.filtered"),ccvmap,coordinates,maskedEdges,maskedLCBs);
#endif
  
#ifdef LCBSTATS
  //Calculate stats
  summaryStats(fglcbsyn,componentMap,coordinates,minprintlength, 
	       numc,minlen,totallen,avglen,maxv,sequence2genome); //using minprintlength
  if(numc>0){
    std::cerr << "LCB summary post-cuts (" << origbreaks << " cuts) " << numc << " min:" << minlen << " coverage:" << totallen << "(" << (float)totallen/allbps << ")" << " avg_bp:" << avglen/numc << " maxv:" << maxv<< std::endl;
  }
  summaryStats(fglcbsyn,componentMap,coordinates,minlength, 
	       numc,minlen,totallen,avglen,maxv,sequence2genome); //using minlength
  if(numc>0){
    std::cerr << "LCB summary post-cuts (" << origbreaks << " cuts) " << numc << " min:" << minlen << " coverage:" << totallen << "(" << (float)totallen/allbps << ")" << " avg_bp:" << avglen/numc << " maxv:" << maxv<< std::endl;
  }
#endif
  
  
  int lcbidx=0;
  
#ifdef DEBUG
  //Preceeding step breakLCBmincut and CC should not introduce
  //bad edges so check predicates
  cutattempts+=1000;
  int morebreaks = breakLCBmincutconnect(componentMap,ccvmap,maskedEdges,g,fglcbsyn,distance,coordinates,seqidxSet,name2vertex,sequence2genome,cutattempts);
  std::cerr << "Num orig breaks " << morebreaks << std::endl;
  assert(morebreaks==0);

  for(std::vector<LCB >::iterator it = componentMap.begin();it!=componentMap.end();++it){
    
    //checkSeqsPerLCB(g,*it)
    std::map<Label,std::set<Label> > seqspergenomeMap; //tracks the number of seqs per genome in an LCB
    std::map<Label,std::set<Label> >::iterator gpos;
    bool inserted;
    property_map < Graph, vertex_label_t >::type vlabelmap = get(vertex_label,g);
    std::cerr << " LCB " << lcbidx << std::endl;
    for(LCB::iterator vit = it->begin();vit!=it->end();++vit){
      std::cerr << " V:" << *vit << std::endl;
      printlabel(get(vertex_orient,g,*vit));
      std::cerr << std::endl;
      for(LabelSet::iterator sit = vlabelmap[*vit].begin();sit!=vlabelmap[*vit].end();++sit){
	//std::cerr << " seqidx:" << *sit << " genomeidx:" << sequence2genome[*sit] << std::endl;
	tie(gpos, inserted) = seqspergenomeMap.insert(std::make_pair(sequence2genome[*sit],std::set<Label>()));
	gpos->second.insert(*sit);
	assert(gpos->second.size()==1);
      }
    }
    
    if(checkLCBGaps(g,*it,ccvmap,coordinates,distance,sequence2genome)){}
    else{
      std::cerr << "Bad gap" << std::endl;
      assert(false);
    }
    if(checkLCBOrient(g,*it,sequence2genome)){}
    else{
      std::cerr << "Misoriented LCB" << std::endl;
      //TODO, add orientation condition to mincut
      //assert(false);
    }
    lcbidx++;
  }
#endif
  
  //
  //(2.4)Attempt to merge lcbs that are adjacent on two or more genomes
  //   and do not introduce rearrangements, gaps
  //
  std::cerr << "Merging adjacent LCBs" << std::endl;
  //Update lcbcoords
  lcbcoords.clear();
#ifdef DEBUG
  for(unsigned int k=0;k<componentMap.size();++k){
    if(componentMap[k].size()>0){
      property_map < LCBSynFilterGraph, vertex_orient_t >::type omap = get(vertex_orient, fglcbsyn);
      property_map < LCBSynFilterGraph, vertex_len_t >::type lmap = get(vertex_len, fglcbsyn);
      int bplen=0;
      unsigned int len = get_LCB_length(componentMap[k],omap,lmap,coordinates,lcbcoords,k,bplen,sequence2genome); 
      assert(len>=0);
    }
  }
#endif

  //Breakpoints are stored in maskedEdges. mergeLCBs clears breakpoints between connected,adjacent and congruent LCBs
  //TODO, put this in loop
  int nummerges=-1;
  int totalnummerges=0;
  while(nummerges!=0){
    nummerges = mergeLCBsGreedy(g,ccvmap,componentMap,lcborientmap,coordinates,maskedEdges,distance,sequence2genome);
    totalnummerges = totalnummerges+nummerges;
  }
  //int origmerges = mergeLCBs(g,ccvmap,componentMap,lcborientmap,coordinates,maskedEdges,distance,sequence2genome);
#ifdef DEBUG
  std::cerr << "Num orig merges " << totalnummerges << std::endl;
#endif
  numComponents = calc_components_undirected(fglcbsyn,g,componentMap,ccvmap,lcborientmap,sequence2genome);  
#ifdef TIMING
  time(&now);
  std::cerr << "TIME_MERGE1:" << now-lasttime << std::endl;
  lasttime=now;
#endif
#ifdef DEBUG
  do_write_graphviz(g, std::string("gout.merge1.dot"),ccvmap,coordinates,maskedEdges,maskedLCBs);
  do_write_graphviz(fglcbsyn, std::string("gout.merge1.dot.filtered"),ccvmap,coordinates,maskedEdges,maskedLCBs);
#endif
#ifdef DEBUG
  //TESTING
  //Ensure merge didn't introduce large gaps
  //origbreaks = breakLCBmincutconnect(componentMap,maskedEdges,g,fglcbsyn,distance,coordinates,seqidxSet,name2vertex,sequence2genome);
  //std::cerr << "Num orig breaks " << origbreaks << std::endl;
  //assert(origbreaks==0);
  //Preceeding step breakLCBmincut and CC should not introduce
  //bad edges so check predicates
  lcbidx=0;
  for(std::vector<LCB >::iterator it = componentMap.begin();it!=componentMap.end();++it){
    if(checkLCBGaps(g,*it,ccvmap,coordinates,distance,sequence2genome)
       && checkLCBOrient(g,*it,sequence2genome)){}
    else{
      //TODO, fix orient check in mincut
      //property_map < Graph, vertex_label_t >::type labelmap = get(vertex_label,g);
      //property_map < Graph, vertex_len_t >::type lenmap = get(vertex_len,g);
      //BitMask longlabelmask=setSpanMask(*it,lenmap,labelmap,sequence2genome);
      //assert(!checkLCBOrient(g,*it,longlabelmask,sequence2genome));
      //assert(false);
    }
    lcbidx++;
  }
#endif
  
#ifdef LCBSTATS
  //Calculate stats
  summaryStats(fglcbsyn,componentMap,coordinates,minprintlength, 
	       numc,minlen,totallen,avglen,maxv,sequence2genome); //using minlength==0
  if(numc>0){
    std::cerr << "LCB summary post-merge (" << totalnummerges << " merges) " << numc << " min:" << minlen << " coverage:" << totallen << "(" << (float)totallen/allbps << ")" << " avg_bp:" << avglen/numc << " maxv:" << maxv<< std::endl;
  }
  summaryStats(fglcbsyn,componentMap,coordinates,minlength, 
	       numc,minlen,totallen,avglen,maxv,sequence2genome); //using minlength
  if(numc>0){
    std::cerr << "LCB summary post-merge (" << totalnummerges << " merges) " << numc << " min:" << minlen << " coverage:" << totallen << "(" << (float)totallen/allbps << ")" << " avg_bp:" << avglen/numc << " maxv:" << maxv<< std::endl;
  }
#endif
  
#ifdef DEBUG
  cutattempts+=1000;
  int newbreaks = breakLCBmincutconnect(componentMap,ccvmap,maskedEdges,g,fglcbsyn,distance,coordinates,seqidxSet,name2vertex,sequence2genome,cutattempts);
  std::cerr << "Breaks after merge " << newbreaks << std::endl;
  assert(newbreaks==0);
#endif
  //(7) Remove breakpoints caused by short LCBs
  unsigned int threshold=shortlcblen;//bp
  std::cerr << "Masking short lcbs <= length " << threshold << std::endl;
  
  
  unsigned int numremoved=0;
  std::vector<LCB> currRemovedLCB = componentMap;
  for(unsigned int k=0;k<componentMap.size();++k){
    if(componentMap[k].size()>0){
      assert(componentMap[k].size()>0);
      //checkLCB(componentMap[k],fglcbsyn);
      property_map < LCBSynFilterGraph, vertex_orient_t >::type omap = get(vertex_orient, fglcbsyn);
      property_map < LCBSynFilterGraph, vertex_len_t >::type lmap = get(vertex_len, fglcbsyn);
      unsigned int len = get_LCB_length(componentMap[k],omap,lmap,coordinates,lcbcoords,k,totallen,sequence2genome); 
#ifdef DEBUG
      std::cerr << "LCB " << k << " len:" << len << std::endl;
#endif
      if(len >=0 && len < threshold){
	//Remove LCB
	removeLCB(componentMap[k],breakpoints,maskedLCBs);
	currRemovedLCB[k].clear();
	numremoved++;
	for(LCB::iterator vit = componentMap[k].begin();vit!=componentMap[k].end();++vit){
	  put(vertex_relorder,g,*vit,len);
	}
      }
    }
  }
#ifdef DEBUG
  std::cerr << "Removed " << numremoved << " LCBs (len<" << threshold << ") containing " 
	    << maskedLCBs.size() << " vertices" << std::endl;
  std::cerr << "Remaining LCBs: " << numComponents-numremoved << std::endl;
#endif
  
  //(8) Update synteny graph
  //to connect vertices that are adjacent when ignoring short/masked LCBs
  //be sure to only add good edges to avoid over-merging clusters
  updateAdjacency(fglcbsyn,
		  g,
		  seqidxSet,
		  coordinates,
		  lcborientmap,
		  distance,
		  maskedEdges,
		  ccvmap,
		  componentMap,
		  sequence2genome);

#ifdef DEBUG
  //std::cerr << "Iteration 1 of CC. Num LCBs: " << numComponents << std::endl;
  do_write_graphviz(g, std::string("gout.dot"),ccvmap,coordinates,maskedEdges,maskedLCBs);
  do_write_graphviz(fglcbsyn, std::string("gout.dot.filtered"),ccvmap,coordinates,maskedEdges,maskedLCBs);
#endif
  
  //
  //(9) Recalculate breakpoints on updated graph
  //breakpoints.clear();
  EdgeSet keepmaskedEdges;
  for(EdgeSet::iterator it = maskedEdges.begin();it!=maskedEdges.end();++it){
    Edge e;
    bool found;
    tie(e,found) = edge(it->first,it->second,g);
    assert(found);
    //BLUE edges are previous cuts
    if(get(edge_category,g,e)==BLUE){
      keepmaskedEdges.insert(*it);
    }
  }
#ifdef DEBUG
  std::cerr << "Keeping " << keepmaskedEdges.size() << " breakpoints" << std::endl;
#endif
  maskedEdges.clear();
  maskedLCBs.clear();     
  //
  //Mark breakpoints with short LCBs "masked"
  markBreakpoints(fglcbsyn,breakpoints,maskedEdges,dummySet,sequence2genome);
  //markBreakpoints(fglcbsyn,breakpoints,maskedEdges,maskedLCBs,sequence2genome);
#ifdef DEBUG
  std::cerr << "Recalc breakpoints. Num:" << breakpoints.size() << std::endl;
#endif
  for(EdgeSet::iterator it = keepmaskedEdges.begin();it!=keepmaskedEdges.end();++it){
    maskedEdges.insert(*it);
  }
  
  maskedLCBs.clear();         
  numComponents = calc_components_undirected(fglcbsyn,g,componentMap,ccvmap,lcborientmap,sequence2genome);
#ifdef TIMING
  time(&now);
  std::cerr << "TIME_MASKSHORT1:" << now-lasttime << std::endl;
  lasttime=now;    
#endif
#ifdef DEBUG
  do_write_graphviz(g, std::string("gout.maskshort.dot"),ccvmap,coordinates,maskedEdges,maskedLCBs);
  do_write_graphviz(fglcbsyn, std::string("gout.maskshort.dot.filtered"),ccvmap,coordinates,maskedEdges,maskedLCBs);
#endif
  if(numComponents==0){
    //No components
    return 0;
  }
#ifdef LCBSTATS
  //Calculate stats
  summaryStats(fglcbsyn,componentMap,coordinates,minprintlength, 
	       numc,minlen,totallen,avglen,maxv,sequence2genome);
  assert(avglen>0);
  assert(numc>0);
  std::cerr << std::endl;
  
  //std::cerr << "Iteration 2 of CC. Num LCBs: " << numComponents << std::endl;
  std::cerr << "LCB summary post-maskshort+merge ("<< numremoved << " LCBs < " << threshold << ") " << numc << " min:" << minlen << " coverage:" << totallen << "(" << (float)totallen/allbps << ")" << " avg_bp:" << avglen/numc << " maxv:" << maxv<< std::endl;
  summaryStats(fglcbsyn,componentMap,coordinates,minlength, 
	       numc,minlen,totallen,avglen,maxv,sequence2genome);
  if(numc>0){
    std::cerr << "LCB summary post-maskshort+merge ("<< numremoved << " LCBs < " << threshold << ") " << numc << " min:" << minlen << " coverage:" << totallen << "(" << (float)totallen/allbps << ")" << " avg_bp:" << avglen/numc << " maxv:" << maxv<< std::endl;
  }
#endif

  itercount=MAXITERS; 
  int nobreaks=0; //number of iterations with no breaks
  while(itercount>0){
    /*
      (10) Break apart components that violate invariants
      For each component/LCB
      -Order the component by projecting blocks onto each member
      sequence in increasing order along the genome
      -Iterator over the projection checking the distance invariant at each iteration
      -If the "gap" between the current block previously seen block > distance (prev.max-curr.min>distance)
      Break/mask all edges that connect the current block with the previous blocks in the ordering.
      Save the broken edges in maskedEdges
    */
#ifdef DEBUG
    std::cerr << "Breaking LCBs that violate contraints" << std::endl;
#endif
  std::cerr << "Partitioning graph to maintain contraints" << std::endl;
  cutattempts+=1000;
  int breaks = breakLCBmincutconnect(componentMap,ccvmap,maskedEdges,g,fglcbsyn,distance,coordinates,seqidxSet,name2vertex,sequence2genome,cutattempts);
#ifdef DEBUG
  std::cerr << "Number breaks " << breaks << std::endl;
#endif
  if(breaks==0){
#ifdef DEBUG
    std::cerr << "Summary mincut shows no breaks necessary, ending iteration at " << MAXITERS-itercount << std::endl;
#endif
    nobreaks++;
    if(nobreaks>MAXSTABLE){
      break;
    }
  }
  //
  //(11) Recalc CC
  //
  numComponents = calc_components_undirected(fglcbsyn,g,componentMap,ccvmap,lcborientmap,sequence2genome);
#ifdef TIMING
  time(&now);
  std::cerr << "TIME_MINCUT2:" << now-lasttime << std::endl;
  lasttime=now;
#endif
#ifdef DEBUG
  do_write_graphviz(g, std::string("gout.mincut2.dot"),ccvmap,coordinates,maskedEdges,maskedLCBs);
  do_write_graphviz(fglcbsyn, std::string("gout.mincut2.dot.filtered"),ccvmap,coordinates,maskedEdges,maskedLCBs);
#endif
#ifdef DEBUG
  int lcbidx=0;
  //Preceeding step breakLCBmincut and CC should not introduce
  //bad edges so check predicates
  for(std::vector<LCB >::iterator it = componentMap.begin();it!=componentMap.end();++it){
    if(checkLCBGaps(g,*it,ccvmap,coordinates,distance,sequence2genome)
       && checkLCBOrient(g,*it,sequence2genome)){}
    else{
      std::cerr << "Bad gap or orient in LCB " << lcbidx << std::endl;
      //Assert(false);
    }
    lcbidx++;
  }
#endif
  VertexSet shortLCBs;
  
#ifdef DEBUG
  //
  //Update cc map to match output
  std::vector<int> xxxc(num_vertices(fglcbsyn));
  int gvlcbnum=0;
  int maskedlcbnum=0;
  for(std::vector<LCB >::iterator it = componentMap.begin();it!=componentMap.end();++it){
    if(it->size()>0){ 
      LCB clcb = *it;
      property_map < LCBSynFilterGraph, vertex_orient_t >::type omap = get(vertex_orient, fglcbsyn);
      property_map < LCBSynFilterGraph, vertex_len_t >::type lmap = get(vertex_len, fglcbsyn);
      unsigned int len = get_LCB_length(clcb,omap,lmap,coordinates,lcbcoords,gvlcbnum,totallen,sequence2genome); 
      if(len>=minlength){
	for(LCB::iterator vit = it->begin();vit!=it->end();++vit){
	  xxxc[*vit]=gvlcbnum;
	}
	gvlcbnum++;
      }
      else{
	maskedlcbnum--;
	for(LCB::iterator vit = it->begin();vit!=it->end();++vit){
	  xxxc[*vit]=maskedlcbnum;
	  shortLCBs.insert(*vit);
	}
      }
    }
  }
  
  do_write_graphviz(fglcbsyn, std::string("gout.dot.postmincut2"),xxxc,coordinates,maskedEdges,shortLCBs);
  do_write_graphviz(g, std::string("gout.dot.postmincut2.all"),xxxc,coordinates,maskedEdges,shortLCBs);
#endif
  
#ifdef LCBSTATS
  
  //Calculate stats      
  unsigned int avglen,minlen;
  int totallen,numc,maxv;
  summaryStats(fglcbsyn,componentMap,coordinates,minprintlength,
	       numc,minlen,totallen,avglen,maxv,sequence2genome);
  assert(avglen>0);
  assert(numc>0);
  std::cerr << "LCB summary post-cuts (" << breaks << " cuts) " << numc << " min:" << minlen << " coverage:" << totallen << "(" << (float)totallen/allbps << ")" << " avg_bp:" << avglen/numc << " maxv:" << maxv<< std::endl;
  summaryStats(fglcbsyn,componentMap,coordinates,minlength,
	       numc,minlen,totallen,avglen,maxv,sequence2genome);
  std::cerr << "LCB summary post-cuts (" << breaks << " cuts) " << numc << " min:" << minlen << " coverage:" << totallen << "(" << (float)totallen/allbps << ")" << " avg_bp:" << avglen/numc << " maxv:" << maxv<< std::endl;
#endif
  
  //(12) Mask short LCBs
  //
  //Mask short LCBs and fix misorients
  unsigned int numremoved=0;
  for(unsigned int k=0;k<componentMap.size();++k){
    if(componentMap[k].size()>0){
      assert(componentMap[k].size()>0);
      property_map < LCBSynFilterGraph, vertex_orient_t >::type omap = get(vertex_orient, fglcbsyn);
      property_map < LCBSynFilterGraph, vertex_len_t >::type lmap = get(vertex_len, fglcbsyn);
      int bplen=0;
      unsigned int len = get_LCB_length(componentMap[k],omap,lmap,coordinates,lcbcoords,k,bplen,sequence2genome); 
      if(len>=0 && len>=shortlcblen){
	//TODO, consider if misoriented vertices should be trimmed or kept till end
	//fixMisOrientedLCBs(g,componentMap[k],maskedLCBs,maskedEdges);
      }
      else{
	//
	//Mask the LCB
	removeLCB(componentMap[k],breakpoints,maskedLCBs);
	numremoved++;
      }
    }
  }
  
  //(13) Connect adjacent LCBs after masking short LCBs
  //
  //Update synteny graph with short LCBs masked 
  //be sure to only add good edges to avoid over-merging clusters
  updateAdjacency(fglcbsyn,
		  g,
		  seqidxSet,
		  coordinates,
		  lcborientmap,
		  distance,
		  maskedEdges,
		  ccvmap,
		  componentMap,
		  sequence2genome);
  
  //(14) Recalculate breakpoints on updated graph
  //breakpoints.clear();
  EdgeSet keepmaskedEdges;
  for(EdgeSet::iterator it = maskedEdges.begin();it!=maskedEdges.end();++it){
    Edge e;
    bool found;
    tie(e,found) = edge(it->first,it->second,g);
    assert(found);
    if(get(edge_category,g,e)==BLUE){
      keepmaskedEdges.insert(*it);
    }
  }
#ifdef DEBUG
  std::cerr << "Keeping " << keepmaskedEdges.size() << " breakpoints" << std::endl;
#endif
  maskedEdges.clear();
  maskedLCBs.clear();     
  //Mark breakpoints on original graph
  markBreakpoints(fglcbsyn,breakpoints,maskedEdges,dummySet,sequence2genome);
  //markBreakpoints(fglcbsyn,breakpoints,maskedEdges,maskedLCBs,sequence2genome);
#ifdef DEBUG
  std::cerr << "Recalc breakpoints. Num:" << breakpoints.size() << std::endl;
#endif
  for(EdgeSet::iterator it = keepmaskedEdges.begin();it!=keepmaskedEdges.end();++it){
    maskedEdges.insert(*it);
  }
  
  maskedLCBs.clear();  
  numComponents = calc_components_undirected(fglcbsyn,g,componentMap,ccvmap,lcborientmap,sequence2genome);
#ifdef TIMING
  time(&now);
  std::cerr << "TIME_MASK2:" << now-lasttime << std::endl;
  lasttime=now;
#endif
#ifdef DEBUG
  do_write_graphviz(g, std::string("gout.mask2.dot"),ccvmap,coordinates,maskedEdges,maskedLCBs);
  do_write_graphviz(fglcbsyn, std::string("gout.mask2.dot.filtered"),ccvmap,coordinates,maskedEdges,maskedLCBs);
#endif 
#ifdef LCBSTATS
  //Calculate stats
  summaryStats(fglcbsyn,componentMap,coordinates,minprintlength,
	       numc,minlen,totallen,avglen,maxv,sequence2genome);
  
  if(numc>0){
    std::cerr << "LCB summary post-maskshort ("<< numremoved << " LCBs < " << threshold << ") " << numc << " min:" << minlen << " coverage:" << totallen << "(" << (float)totallen/allbps << ")" << " avg_bp:" << avglen/numc << " maxv:" << maxv<< std::endl;
  } 
  summaryStats(fglcbsyn,componentMap,coordinates,minlength,
	       numc,minlen,totallen,avglen,maxv,sequence2genome);
  
  if(numc>0){
    std::cerr << "LCB summary post-maskshort ("<< numremoved << " LCBs < " << threshold << ") " << numc << " min:" << minlen << " coverage:" << totallen << "(" << (float)totallen/allbps << ")" << " avg_bp:" << avglen/numc << " maxv:" << maxv<< std::endl;
  }
#endif
  std::cerr << "Merging adjacent LCBs " << std::endl;

  //
  //(XX)Attempt to merge LCBs
  int nummerges=-1;
  int totalnummerges=0;
  while(nummerges!=0){
    nummerges = mergeLCBsGreedy(g,ccvmap,componentMap,lcborientmap,coordinates,maskedEdges,distance,sequence2genome);
    totalnummerges = totalnummerges+nummerges;
    //TODO, consider fixing misoriented vertices here
    //for(unsigned int k=0;k<componentMap.size();++k){
    //fixMisOrientedLCBs(g,componentMap[k],maskedLCBs,maskedEdges);
    //}
    
    //
    //(XX)Recalculate connected components
#ifdef DEBUG
    std::cerr << "Num merges:" << nummerges << std::endl;
#endif
    numComponents = calc_components_undirected(fglcbsyn,g,componentMap,ccvmap,lcborientmap,sequence2genome);      
#ifdef DEBUG
    std::cerr << "Recalc components: " << numComponents << std::endl;
#endif
  }
#ifdef TIMING
  time(&now);
  std::cerr << "TIME_MERGE2:" << now-lasttime << std::endl;
  lasttime=now;
#endif
#ifdef DEBUG
  do_write_graphviz(g, std::string("gout.merge2.dot"),ccvmap,coordinates,maskedEdges,maskedLCBs);
  do_write_graphviz(fglcbsyn, std::string("gout.merge2.dot.filtered"),ccvmap,coordinates,maskedEdges,maskedLCBs);
#endif      
#ifdef DEBUG
  for(std::vector<LCB >::iterator it = componentMap.begin();it!=componentMap.end();++it){
    if(it->size()>0){ 
      LCB clcb = *it;
      property_map < LCBSynFilterGraph, vertex_orient_t >::type omap = get(vertex_orient, fglcbsyn);
      property_map < LCBSynFilterGraph, vertex_len_t >::type lmap = get(vertex_len, fglcbsyn);
      unsigned int len = get_LCB_length(clcb,omap,lmap,coordinates,lcbcoords,0,totallen,sequence2genome); 
      if(len>=minlength){
      }
      else{
	for(LCB::iterator vit = it->begin();vit!=it->end();++vit){
	  //shortLCBs.insert(*vit);
	}
      }
    }
  }
  do_write_graphviz(g, std::string("gout.dot.postimerge."+lexical_cast<std::string>(MAXITERS-itercount)+".all"),ccvmap,coordinates,maskedEdges,shortLCBs);
#endif
#ifdef LCBSTATS
  //Calculate stats
  summaryStats(fglcbsyn,componentMap,coordinates,minprintlength,
	       numc,minlen,totallen,avglen,maxv,sequence2genome);
  if(numc>0){
    std::cerr << "LCB summary post-maskshort+merge ("<< totalnummerges << " merges) " << numc << " min:" << minlen << " coverage:" << totallen << "(" << (float)totallen/allbps << ")" << " avg_bp:" << avglen/numc << " maxv:" << maxv<< std::endl;
  }
  summaryStats(fglcbsyn,componentMap,coordinates,minlength,
	       numc,minlen,totallen,avglen,maxv,sequence2genome);
  if(numc>0){
    std::cerr << "LCB summary post-maskshort+merge ("<< totalnummerges << " merges) " << numc << " min:" << minlen << " coverage:" << totallen << "(" << (float)totallen/allbps << ")" << " avg_bp:" << avglen/numc << " maxv:" << maxv<< std::endl;
  }
#endif
  itercount--;
  }// end numiters
  
  //(9) Remove misoriented vertices
  VertexSet shortLCBs;
  shortLCBs.clear();
  for(std::vector<LCB >::iterator it = componentMap.begin();it!=componentMap.end();++it){
    if(it->size()>0){ 
      LCB clcb = *it;
      property_map < LCBSynFilterGraph, vertex_orient_t >::type omap = get(vertex_orient, fglcbsyn);
      property_map < LCBSynFilterGraph, vertex_len_t >::type lmap = get(vertex_len, fglcbsyn);
      unsigned int len = get_LCB_length(clcb,omap,lmap,coordinates,lcbcoords,0,totallen,sequence2genome); 
      if(len>=shortlcblen){
	  fixMisOrientedLCBs(g,clcb,maskedLCBs,maskedEdges,sequence2genome);
      }
      else{
#ifdef DEBUG
	for(LCB::iterator vit = it->begin();vit!=it->end();++vit){
	  shortLCBs.insert(*vit);
	}
#endif
      }
    }
  }
  
  numComponents = calc_components_undirected(fglcbsyn,g,componentMap,ccvmap,lcborientmap,sequence2genome);
  
  //Sanity check to ensure no large gaps
  //TODO
  //Try and avoid this condition
  //Currently remvong LCBs in fixMisOriented may introduce gaps > threshold
  //so we need to recut here. 
  int breaks=-1;
  while(breaks!=0){
    cutattempts += 1000;
    breaks=breakLCBmincutconnect(componentMap,ccvmap,maskedEdges,g,fglcbsyn,distance,coordinates,seqidxSet,name2vertex,sequence2genome,cutattempts);
#ifdef DEBUG
    std::cerr << "Num final breaks " << breaks << std::endl;
#endif
    numComponents = calc_components_undirected(fglcbsyn,g,componentMap,ccvmap,lcborientmap,sequence2genome);
  }
  
#ifdef DEBUG
  cutattempts+=1000;
  breaks=breakLCBmincutconnect(componentMap,ccvmap,maskedEdges,g,fglcbsyn,distance,coordinates,seqidxSet,name2vertex,sequence2genome,cutattempts);
  assert(breaks==0);
#endif
#ifdef LCBSTATS
  //Calculate stats
  summaryStats(fglcbsyn,componentMap,coordinates,minprintlength,
	       numc,minlen,totallen,avglen,maxv,sequence2genome);
  if(numc>0){
    std::cerr << "LCB summary final " << numc << " min:" << minlen 
	      << " coverage:" << totallen << "(" << (float)totallen/allbps << ")" << " avg_bp:" << avglen/numc << " maxv:" << maxv<< std::endl;
  }
  summaryStats(fglcbsyn,componentMap,coordinates,minlength,
	       numc,minlen,totallen,avglen,maxv,sequence2genome);
  if(numc>0){
    std::cerr << "LCB summary final " << numc << " min:" << minlen 
	      << " coverage:" << totallen << "(" << (float)totallen/allbps << ")" << " avg_bp:" << avglen/numc << " maxv:" << maxv<< std::endl;
  }
#endif
#ifdef DEBUG
  do_write_graphviz(g, std::string("gout.dot.final.all"),ccvmap,coordinates,maskedEdges,shortLCBs);
#endif
#ifdef TIMING
  time(&now);
  std::cerr << "TIME_MINCUT3:" << now-lasttime << std::endl;
  lasttime=now;
#endif
#ifdef DEBUG
	do_write_graphviz(g, std::string("gout.final.dot"),ccvmap,coordinates,maskedEdges,maskedLCBs);
	do_write_graphviz(fglcbsyn, std::string("gout.final.dot.filtered"),ccvmap,coordinates,maskedEdges,maskedLCBs);
#endif
  //
  //
  //Optionally remove overlapping LCBs, keep longest LCBs spanning each region
  int idx=0;
  //Index of lcb in componentMap
  lcbidx=0;
  std::map<int,int> lcboverlapMap; //lcbid->longer_overlapping_lcbid
  std::set<int> bestLCBs; //set of LCBs that are longest over at least one genomic segment
  std::map<int,int> lcblenMap; //lcbid->max_seq_span
  if(removeoverlaps){
    std::cerr << "Sorting LCBs along each seq. Num seqs " << seqidxSet.size() << std::endl;
    property_map < Graph, vertex_orient_t >::type orientmap = get(vertex_orient,g);
    property_map < Graph, vertex_len_t >::type lenmap = get(vertex_len,g);
    typedef iloc TLoc;
    LCB::iterator it;
    std::vector<std::vector<TLoc> > olaplcbs;
    olaplcbs.resize(seqidxSet.size()+1);
    for(std::vector<LCB >::iterator it = componentMap.begin();it!=componentMap.end();++it){
      std::cerr << "Looking at LCB " << lcbidx << " of size " << it->size() << std::endl;
      std::set<Label> currseqs;
      property_map < LCBSynFilterGraph, vertex_orient_t >::type omap = get(vertex_orient, fglcbsyn);
      property_map < LCBSynFilterGraph, vertex_len_t >::type lmap = get(vertex_len, fglcbsyn);
      for(LCB::iterator lit = it->begin();lit!=it->end();++lit){
	OrientedLabelSet::iterator it2_end = omap[*lit].end();
	for(OrientedLabelSet::iterator it2 = omap[*lit].begin();it2!=it2_end;++it2){
	  Label seqidx = it2->first;
	  currseqs.insert(seqidx);
	  }
      }
      unsigned int len = get_LCB_length(*it,omap,lmap,coordinates,lcbcoords,lcbidx,totallen,sequence2genome); 
      lcblenMap[lcbidx] = len;
      for(std::set<Label>::iterator it2 = currseqs.begin();it2!=currseqs.end();++it2){
	Label seqidx = *it2;
	assert(lcbcoords.find(std::make_pair(lcbidx,seqidx))!=lcbcoords.end());
	//Label genomeidx = sequence2genome[seqidx];
	TLoc t1,t2;
	t1.first = lcbcoords[std::make_pair(lcbidx,seqidx)].second;
	t1.second = 0;
	t1.blocknum=lcbidx;
	t2.first = lcbcoords[std::make_pair(lcbidx,seqidx)].first;
	t2.second = 1;
	t2.blocknum=lcbidx;
	if(t1.first-t2.first > 0){//{(int)shortlcblen){
	  std::cerr << "lcbidx: " << lcbidx << " sidx: " << seqidx << " " << olaplcbs.size() << std::endl;
	  assert(seqidx<olaplcbs.size());
	  olaplcbs[seqidx].push_back(t1);
	  olaplcbs[seqidx].push_back(t2);
	}
      }
      lcbidx++;
      idx++;  
    }
    
    std::cerr << "Iterating over LCBs to look for overlaps" << std::endl;
    for(int i=0;i<(int)seqidxSet.size()+1;i++){
      int open=0;
      assert(i<(int)olaplcbs.size());
      std::vector<TLoc> &ait=olaplcbs[i];
      std::cerr << "Sorting on seq " << i << std::endl;
      sort(ait.begin(),ait.end(),poscmp<TLoc>());
      std::cerr << "sorted" << std::endl;
      std::set<int> currlcbs;
      for(std::vector<TLoc>::iterator pit = ait.begin();pit!=ait.end();pit++){
	int currlen = lcblenMap[pit->blocknum];
	if(pit->second>0){
	  int longestlen = lcblenMap[pit->blocknum];
	  int bestlcb = pit->blocknum;
	  std::cerr << "LCB open " << pit->blocknum << " len " << currlen << std::endl;
	  if(open>0){
	    //in overlap
	    //only remove short LCBs that are overlapped by longer lcbs
	    //if(currlen<shortlcblen){
	    std::cerr << "Number overlaps " << currlcbs.size() << std::endl;
	    for(std::set<int>::iterator cit=currlcbs.begin();cit!=currlcbs.end();++cit){
	      std::cerr << "Overlapping " << *cit << " len " << lcblenMap[*cit] << std::endl;
	      assert(*cit!=pit->blocknum);
	      if(lcblenMap[*cit]>longestlen){
		longestlen = lcblenMap[*cit];
		bestlcb = *cit;
	      }
	      //overlapping lcb > current lcb
	      if(lcblenMap[*cit]>currlen){
		if(lcboverlapMap.find(pit->blocknum)!=lcboverlapMap.end()){
		  //
		  if(lcblenMap[*cit]>lcblenMap[lcboverlapMap[pit->blocknum]]){
		    lcboverlapMap[pit->blocknum] = *cit;
		  }
		}
		else{
		  lcboverlapMap[pit->blocknum] = *cit;
		}
	      }
	    }
	  }
	  open++;
	  assert(currlcbs.find(pit->blocknum)==currlcbs.end());
	  currlcbs.insert(pit->blocknum);
	  bestLCBs.insert(bestlcb);
	  std::cerr << "opened " << pit->blocknum << std::endl;
	}
	else{
	  open--;
	  assert(currlcbs.find(pit->blocknum)!=currlcbs.end());
	  assert(currlcbs.size()>0);
	  currlcbs.erase(pit->blocknum);
	  std::cerr << "closed " << pit->blocknum << std::endl;
	  int longestlen = lcblenMap[pit->blocknum];
	  int bestlcb = pit->blocknum;
	  if(open){
	    assert(currlcbs.size()>0);
	  }
	  for(std::set<int>::iterator cit=currlcbs.begin();cit!=currlcbs.end();++cit){
	    if(lcblenMap[*cit]>longestlen){
	      longestlen = lcblenMap[*cit];
	      bestlcb = *cit;
	    }
	    if(lcblenMap[*cit]>currlen){
	      if(lcboverlapMap.find(pit->blocknum)!=lcboverlapMap.end()){
		if(lcblenMap[*cit]>lcblenMap[lcboverlapMap[pit->blocknum]]){
		  lcboverlapMap[pit->blocknum] = *cit;
		}
	      }
	      else{
		lcboverlapMap[pit->blocknum] = *cit;
	      }
	    }
	  }
	  bestLCBs.insert(bestlcb);
	}
      }
    }
    
  }
  lcbidx=0;
  std::vector<LCB > validLCBs;
    if(removeoverlaps){
      for(std::vector<LCB >::iterator it = componentMap.begin();it!=componentMap.end();++it){
	if(bestLCBs.find(lcbidx)==bestLCBs.end()){ 
	  std::cerr << "LCB idx: " << lcbidx;
	  if(lcboverlapMap.find(lcbidx)!=lcboverlapMap.end()){
	    std::cerr << " overlaps " << lcboverlapMap[lcbidx];
	  }
	  LCB newlcb;
	  newlcb.insert(newlcb.end(),it->begin(),it->end());
	  newlcb.insert(newlcb.end(),componentMap[lcboverlapMap[lcbidx]].begin(),componentMap[lcboverlapMap[lcbidx]].end());
	  validLCBs.push_back(newlcb);
	}
	else{
	  validLCBs.push_back(*it);
	}
	lcbidx++;
      }
    }
    else{
      validLCBs = componentMap;
    }
#ifdef LCBSTATS
    //Calculate stats
    summaryStats(fglcbsyn,validLCBs,coordinates,minprintlength,
		 numc,minlen,totallen,avglen,maxv,sequence2genome);
    if(numc>0){
      std::cerr << "LCB summary final post-processing " << numc << " min:" << minlen 
		<< " coverage:" << totallen << "(" << (float)totallen/allbps << ")" << " avg_bp:" << avglen/numc << " maxv:" << maxv<< std::endl;
    }
    summaryStats(fglcbsyn,validLCBs,coordinates,minlength,
		 numc,minlen,totallen,avglen,maxv,sequence2genome);
    if(numc>0){
      std::cerr << "LCB summary final post-processing " << numc << " min:" << minlen 
		<< " coverage:" << totallen << "(" << (float)totallen/allbps << ")" << " avg_bp:" << avglen/numc << " maxv:" << maxv<< std::endl;
    }
#endif

  //Write out LCBs
  //Format is 2 lines per LCB
  //I seq1 orient1 coords1 ... seqN orientN coordsN
  //V feat1 feat2 .... featN
    assert((unsigned int)numComponents==componentMap.size());
    lcbidx=0;
    idx=0;
    unsigned int maxlcblen=0;
    property_map < Graph, vertex_orient_t >::type vmap = get(vertex_orient,g);
    //for(std::vector<LCB >::iterator it = componentMap.begin();it!=componentMap.end();++it){
    for(std::vector<LCB >::iterator it = validLCBs.begin();it!=validLCBs.end();++it){
      if(it->size()>0){ 
	idx++;
	LCB clcb = *it;
	unsigned int len;
	property_map < LCBSynFilterGraph, vertex_orient_t >::type omap = get(vertex_orient, fglcbsyn);
	property_map < LCBSynFilterGraph, vertex_len_t >::type lmap = get(vertex_len, fglcbsyn);
	unsigned int nlen = get_LCB_length(clcb,omap,lmap,coordinates,lcbcoords,idx,totallen,sequence2genome);
	len=nlen;
	if(len>=minprintlength){
	  maxlcblen = (len>maxlcblen ? len : maxlcblen);
	  unsigned int numcomps=0; 
	  if(checkLCBOrient(fglcbsyn,*it,sequence2genome)){
	    if(checkLCBGaps(g,*it,ccvmap,coordinates,distance,sequence2genome)){
	      std::cout << "I ";
	      //
	      //Save mask for the LCB
	      BitMask labelsmask;
	      BitMask orientmask;
	      std::vector<Vertex> badV;
	      tie(labelsmask,orientmask) = setLCBOrient(g,*it,badV,sequence2genome);
	      SeqSet currlabelset;
	      for(LCB::iterator vit = it->begin();vit!=it->end();++vit){
		currlabelset.insert(vmap[*vit].begin(),vmap[*vit].end());
	      }
	      for(OrientedLabelSet::iterator oit = currlabelset.begin();oit != currlabelset.end();++oit){
		Label seqidx = oit->first;
		Label genomeidx = sequence2genome[seqidx];
		assert(labelsmask.test(genomeidx));
		std::cout << index2sequence[seqidx] <<" " << (orientmask.test(genomeidx) ? '+' : '-') << " ";
		std::cout << lcbcoords[std::make_pair(idx,seqidx)].first << "-" << lcbcoords[std::make_pair(idx,seqidx)].second << " ";
	      }
	      std::cout << " ;" << std::endl;
	      std::cout << "V ";
	      for(LCB::iterator vit = it->begin();vit!=it->end();++vit){
		std::cout << get(vertex_name,fglcbsyn,*vit) << " ";
		numcomps++;
	      }
	      std::cout << " ;" << std::endl; 
	    }
	    else{
	      std::cerr << "SKIPPING LCB:" << idx << " Bad gap" << std::endl;
	      assert(false);
	    }
	  }
	  else{
	    std::cerr << "BAD LCB:" << idx << " Mis-matched lable orientation" << std::endl;
	    assert(false);
	  }
	}
	else{
	  std::cerr << "SKIPPING LCB:" << lcbidx << " len:" << nlen << " < " << minprintlength << std::endl;
	}
      }
      lcbidx++;
    }
    std::cerr << "Max LCB length " << maxlcblen << std::endl;
#ifdef TIMING
    time(&now);
    std::cerr << "TIME_POSTPROC:" << now-lasttime << std::endl;
    lasttime=now;    
#endif
    return 0;
}


//################
//General utilities
//
//

unsigned int getIntervalDist(int s1, int e1, int s2, int e2){
  //Contained
  if(s1>s2 && s1<e2){
    return 0;
  }
  else{
    if(s2>s1 && s2<e1){
      return 0;
    }
    else{
      if(s1<s2){
	assert(s2-s1>=0);
	return (unsigned int)s2-s1;
      }
      else{
	assert(s1-s2>=0);
	return (unsigned int)s1-s2;
      }
    }
  }
}


// Graph filters
//
//Filter edges and exclude all edges that don't match sequences in matchlabel
//Fast edge filter using bit masks predefined on mindist in setedgesmask
template <typename EdgeLabelMap, typename LabelContainer>
struct distance_label_filter_bv {
  distance_label_filter_bv() { }
  distance_label_filter_bv(EdgeLabelMap label, LabelContainer &lc) 
    : m_label(label),matchlabels(lc) { 
    
  }
  template <typename Edgetype>
  bool operator()(const Edgetype& e) const {
    if((m_label[e]&matchlabels)==matchlabels){
      return true;
    }
    else
      return false;
  }
  EdgeLabelMap m_label;
  LabelContainer matchlabels;
};

template <typename TEdgeLabelMap, typename TLabelSet>
struct edge_label_filter {
  edge_label_filter() { }
  edge_label_filter(TEdgeLabelMap l, 
		    TLabelSet s)
    : labelmap(l),
      matchlabels(s){}
  template <typename Edgetype>
  bool operator()(const Edgetype& e) const {
    assert(matchlabels.size()!=0);
    if(std::includes(matchlabels.begin(),matchlabels.end(),labelmap[e].begin(),labelmap[e].end())){
      return true;
    }
  }
  TEdgeLabelMap  labelmap;
  TLabelSet  matchlabels;
  
};

template <typename TGraph>
struct snode_efilter {
  snode_efilter()
    :G(NULL),m_snodes(NULL)
  {}
  snode_efilter(std::set<typename TGraph::vertex_descriptor> *m,TGraph *gin)
    :m_snodes(m),G(gin)
  {}
  template <typename Edgetype>
  bool operator()(const Edgetype& e) const {
    if(m_snodes->find(source(e,*G))==m_snodes->end()
       && m_snodes->find(target(e,*G))==m_snodes->end()){
      return false;
    }
    else{
      return true;
    }
  }
  std::set<typename TGraph::vertex_descriptor> *m_snodes;
  TGraph *G;
};


/*
//
//Filter vertices and exclude all vertices that do not contain a subset of sequences defined by matchlabels
//and (optional) orientation matchorient
//Passing a empty matchorient container ignores filtering by orientation
//Implemented using bitmasks preset based on distance within setvertexmasks
//TODO
//Refactor to another container besides BitMask and ensure proper filtering
template <typename VertexLabelMaskMap, typename VertexOrientMaskMap, typename OrientContainer>
struct orient_filter_bv {
  orient_filter_bv() { }
  orient_filter_bv(VertexLabelMaskMap label, 
		   OrientContainer &lc, 
		   VertexOrientMaskMap omask, 
		   OrientContainer &oc, 
		   OrientContainer &roc) : m_label(label),
					   matchlabels(lc),
					   m_orientmask(omask),
					   matchorient(oc),
					   matchrevorient(roc) { }
  template <typename Vertextype>
  bool operator()(const Vertextype& v) const {
    assert(matchlabels!=0);
    if((m_label[v]&matchlabels)==matchlabels){ 
      OrientContainer orientmask = (m_orientmask[v]&matchlabels);
      //!matchorient.any() means all - orients, which is disallowed
      //Using this state to specify shortcircuit "ignore orient mask"
      if(!matchorient.any() || orientmask==matchorient){
	return true;
      }
      else{
	if(!matchrevorient.any() || orientmask==matchrevorient){
	  return true;
	}
	else{
	  return false;
	}
      }
    }
    else{
      return false;
    }
  }

  VertexLabelMaskMap m_label;
  OrientContainer matchlabels;
  VertexOrientMaskMap m_orientmask;
  OrientContainer matchorient;
  OrientContainer matchrevorient;
};
*/

template <typename TVertexLabelMap, typename TLabelSet>
struct vertex_label_filter {
  vertex_label_filter() { }
  vertex_label_filter(TVertexLabelMap l, 
		      TLabelSet s)
    : labelmap(l),
      matchlabels(s){}
  template <typename Vertextype>
  bool operator()(const Vertextype& v) const {
    assert(matchlabels.size()!=0);
    //assert((labelmap[v].find(*(matchlabels.begin())) != labelmap[v].end())
    //==
    //(std::includes(labelmap[v].begin(),labelmap[v].end(),matchlabels.begin(),matchlabels.end())));
    if(labelmap[v].find(*(matchlabels.begin())) != labelmap[v].end()){
      //std::includes(labelmap[v].begin(),labelmap[v].end(),matchlabels.begin(),matchlabels.end())){
      return true;
    }
    else{
      return false;
    }
  }
  TVertexLabelMap  labelmap;
  TLabelSet  matchlabels;
  
};

//
//Filter edges and exclude all edges in EdgeSet
//where EdgeSet is a pair of vertices
template <typename TGraph>
struct synbp_edge_filter {
  typedef typename boost::graph_traits<TGraph>::edge_descriptor Edge;
  typedef typename boost::graph_traits<TGraph>::vertex_descriptor Vertex;
  synbp_edge_filter()
    :G(NULL)
  {}
  synbp_edge_filter(EdgeSet *m, TGraph *gin)
    :maskededges(m),G(gin)
  {}
  template <typename Edgetype>
  bool operator()(const Edgetype& e) const {
    if(maskededges->find(std::make_pair(source(e,*G),target(e,*G)))!=maskededges->end()
       || maskededges->find(std::make_pair(target(e,*G),source(e,*G)))!=maskededges->end()){      
      //TODO
      //This additional check target,source should not be necessary
      return false;
    }
    else{
      return true;
    }
  }
  EdgeSet *maskededges;
  TGraph *G;
};

//
//Filter edges and exclude all incident edges to vertices in VertexSet
template <typename TGraph>
struct LCB_edge_filter {
  typedef typename boost::graph_traits<TGraph>::edge_descriptor Edge;
  typedef typename boost::graph_traits<TGraph>::vertex_descriptor Vertex;
  LCB_edge_filter()
    :G(NULL)
  {}
  LCB_edge_filter(VertexSet *m, TGraph *gin)
    :maskedvertices(m),G(gin)
  {}
  template <typename Edgetype>
  bool operator()(const Edgetype& e) const {
    if(maskedvertices->find(source(e,*G))!=maskedvertices->end()
       ||
       maskedvertices->find(target(e,*G))!=maskedvertices->end()){
      return false;
    }
    else{
      return true;
    }
  }
  VertexSet *maskedvertices;
  TGraph *G;
};

//
//Filter vertices and exclude all vertices in VertexSet
template <typename TGraph>
struct LCB_vertex_filter {
  LCB_vertex_filter() { }
  LCB_vertex_filter(VertexSet *m)
    :maskedvertices(m)
  {}
  template <typename Vertextype>
  bool operator()(const Vertextype& v) const {
    if(maskedvertices->find(v)!=maskedvertices->end()){
      return false;
    }
    else{
      return true;
    }
  }
  VertexSet *maskedvertices;
};

//
//Create compound edge filter by chaining two edge filters
template <typename TFilter1, typename TFilter2>
struct compound_edge_filter{
  compound_edge_filter(){}
  compound_edge_filter(TFilter1 &tf1, TFilter2 &tf2)
    :filter1(tf1),filter2(tf2)
  {}
  template <typename Edgetype>
  bool operator()(const Edgetype& e) const {
    if(filter1(e) && filter2(e)){
      assert(filter1(e));
      assert(filter2(e));
    }
    else{
      assert(!filter1(e)||!filter2(e));
    }
    return filter1(e) && filter2(e);
  }
  TFilter1 filter1;
  TFilter2 filter2;
};

//
//Create compound vertex filter by chaining two vertex filters
template <typename TFilter1, typename TFilter2>
struct compound_vertex_filter{
  compound_vertex_filter(){}
  compound_vertex_filter(TFilter1 &tf1, TFilter2 &tf2)
    :filter1(tf1),filter2(tf2)
  {} 
  template <typename Vertextype>
  bool operator()(const Vertextype& v) const {
    if(filter1(v) && filter2(v)){
      assert(filter1(v));
      assert(filter2(v));
    }
    else{
      assert(!filter1(v)||!filter2(v));
    }
    return filter1(v) && filter2(v);
  }
  TFilter1 filter1;
  TFilter2 filter2;
};

//
//Define LCB as a set of connected vertices
typedef std::vector<Vertex> LCB;
typedef std::map<std::pair<unsigned int,Label>,Interval > LCBLabelIntervalMap; 

//
//Graph types
//Synteny graph contains connected subgraphs, each an LCB
typedef filtered_graph<Graph, 
		       synbp_edge_filter<Graph> >  FilterSynGraph; 
typedef filtered_graph<Graph, 
		       LCB_edge_filter<Graph>, 
		       LCB_vertex_filter<Graph> > FilterLCBGraph;
//Synteny graph that supports masking/filtering of LCBs
typedef filtered_graph<Graph, 
		       compound_edge_filter<LCB_edge_filter<Graph>, synbp_edge_filter<Graph> >, 
		       LCB_vertex_filter<Graph> > LCBSynFilterGraph; 


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

class lencmp
{
public:
  lencmp(std::map<int,int> & m)
    :lenmap(&m)
  {}
  bool operator()( const int i1, const int i2) const {
    assert(lenmap->find(i1)!=lenmap->end());
    assert(lenmap->find(i2)!=lenmap->end());

    if (lenmap->find(i1)->second < lenmap->find(i2)->second){
      return true;
    }
    else{
      return false;
    }
  }
  std::map<int,int> *lenmap;
};




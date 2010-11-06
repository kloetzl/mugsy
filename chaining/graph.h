//Types
//############
//BGL requires crazy amount of code to define custom graph properties
//
//Edge properties
struct label_t {
  typedef edge_property_tag kind;
};

struct genome_t { 
    typedef edge_property_tag kind;
};

//TODO f
//Make edge_labelmask set based
//on genomeidx rather than seqidx
struct labelmask_t{
  typedef edge_property_tag kind;
};

struct visted_t{
  typedef edge_property_tag kind;
};

struct stringname_t{
  typedef edge_property_tag kind;
};

//Vertex properties
struct orient_t{
  typedef vertex_property_tag kind;
};

struct orientmask_t{
  typedef vertex_property_tag kind;
};

struct vlabelmask_t{
  typedef vertex_property_tag kind;
};

struct chains_t{
  typedef vertex_property_tag kind;
};

struct relorder_t{
  typedef vertex_property_tag kind;
};

//vertex_name_t and edge_weight_t are already defined by default
enum edge_label_t { edge_label = 10001 };
enum edge_labelmask_t { edge_labelmask = 10004 };
enum edge_visited_t { edge_visited = 10005 };
enum edge_stringname_t { edge_stringname = 10011 };
enum edge_category_t { edge_category = 10012 };
enum vertex_orient_t { vertex_orient = 10006 };
enum vertex_label_t { vertex_label = 10016 };
enum vertex_genome_t { vertex_genome = 10015 };
enum vertex_vlabelmask_t { vertex_vlabelmask = 10007 };
enum vertex_orientmask_t { vertex_orientmask = 10008 };
enum vertex_len_t { vertex_len = 10009 };
enum vertex_relorder_t { vertex_relorder = 10010 };

namespace boost {
  BOOST_INSTALL_PROPERTY(edge, label);
  BOOST_INSTALL_PROPERTY(edge, labelmask);
  BOOST_INSTALL_PROPERTY(edge, visited);
  BOOST_INSTALL_PROPERTY(edge, stringname);
  BOOST_INSTALL_PROPERTY(edge, category);
  BOOST_INSTALL_PROPERTY(vertex, relorder);
  BOOST_INSTALL_PROPERTY(vertex, label);
  BOOST_INSTALL_PROPERTY(vertex, vlabelmask);
  BOOST_INSTALL_PROPERTY(vertex, orientmask);
  BOOST_INSTALL_PROPERTY(vertex, orient);
  BOOST_INSTALL_PROPERTY(vertex, genome);
  BOOST_INSTALL_PROPERTY(vertex, len);
}
//End custom properties code
//###########

//Label is an index that corresponds to a genome sequence for complete
//genomes or a species for incomplete genomes.  using a short limits
//to 65,535 labels in an attempt to save some space
//typedef unsigned short int Label;
typedef unsigned int Label;

//DNA sequence orientation -,+ == false,true
typedef bool Orientation;

//Label,distance map specifies the proximity between two
//anchors/blocks along a sequence whose index is Label
typedef std::map<Label,int> LabelMap;
//typedef __gnu_cxx::hash_map<Label,int> LabelMap;

//Set of labels
//typedef std::set<Label> LabelSet;
typedef boost::unordered_set<Label> LabelSet;

//MAXGENOMES,BITMAX: Critical parameters that limit the number of genome labels
//bitset is used for fast pattern matching of subsets
//Increasing the size of this parameter will degrade performance
//even for small numbers of sequences
//TODO: Refactor. Compare or replace with use of STL includes,set_intersection,...
//or use boost::dynamic_bitset
#define BITMAX MAXGENOMES
typedef std::bitset<BITMAX> BitMask;

typedef pair<Label,Orientation> OrientedLabel;

struct orientedlabelhasher {
  size_t operator()(const OrientedLabel& v) const { return hash<Label>()(v.first); }
};

//faster for bitsets that can fit in ulong
//otherwise will throw an overflow exception
struct bitsethasher_ulong {
  size_t operator()(const BitMask& v) const { return hash<long unsigned int>()(v.to_ulong()); }
};

//general for bitsets of any size
struct bitsethasher_string {
  size_t operator()(const BitMask& v) const { return hash<std::string>()(v.to_string()); }
};

struct orientedlabelcmp {
  bool operator()( const OrientedLabel& s1, const OrientedLabel& s2 ) const {
    if(s1.first == s2.first){
      return s1.second < s2.second;
    }
    else{
      return s1.first < s2.first;
    }
  }
};

struct hasheq
{
  bool operator()(const OrientedLabel& s1, const OrientedLabel& s2) const
  {
    return s1==s2;
  }
};

//Edges in the anchor graph are marked with a set of labels called an
//OrientedLabelSet.  The cardinality of this set is the number of
//member sequences. Each sequence is labeled by an integer index (type
//Label) and a boolean orientation (type Orientation) which are paired
//to defined an OrientedLabel.
//
//TODO objects of this type are copied all over the place in the current impl
//need to refactor to improve performance
typedef std::set<OrientedLabel,orientedlabelcmp > OrientedLabelSet;
//typedef boost::unordered_set<OrientedLabel> OrientedLabelSet;
//typedef __gnu_cxx::hash_set<OrientedLabel,orientedlabelhasher> OrientedLabelSet;

//Names of blocks/anchors and graph vertices: VertexName, VertexID
//The program input includes a set of anchors across two or more
//genomes.  These anchors are also referred to as blocks.  Each block
//is stored as a vertex in a graph.  The identifier provided for each
//block in the user provided input is stored as the VertexName. It is
//assumed that the VertexName is a unique identifier for a block.  An
//additional internal identifier for each block, VertexID, is used by
//the graph library.  For a given block, the VertexName and VertexID
//may not be equivalent.  It is also possible to change the typedef
//for VertexName to std::string to support string names for blocks.
typedef unsigned int VertexName;
typedef unsigned int VertexID; //TODO, reconcile,replace with Graph::vertex_descriptor

//Coordinates and Intervals
typedef int Coordinate;
typedef std::pair<Coordinate,Coordinate> Interval;


//typedef BitMask SeqSet;
typedef OrientedLabelSet SeqSet;
//bit set per genome
typedef LabelSet GenomeSet;

//BGL requires properties in this nested format
//VertexProperties
typedef property<vertex_genome_t, GenomeSet> VertexGenome;
typedef property<vertex_relorder_t, int, VertexGenome> VertexRelOrder;
typedef property<vertex_len_t, int, VertexRelOrder> VertexLen;
typedef property<vertex_orientmask_t, BitMask, VertexLen> VertexOrientMask;
typedef property<vertex_vlabelmask_t, BitMask, VertexOrientMask> VertexLabelMask;
typedef property<vertex_orient_t, SeqSet, VertexLabelMask> VertexOrientedLabel;
typedef property<vertex_label_t, LabelSet, VertexOrientedLabel> VertexLabel;
typedef property<vertex_name_t, VertexName, VertexLabel> VertexProperties;

//Define graph properties
//TODO
//Replace edge properties, such as EdgeMask and LabelMap, with index to save space
//typedef property<edge_category_t, std::string, property<edge_weight_t,int> > EdgeCategory;
//typedef property<edge_stringname_t, std::string> EdgeStringName;
//typedef property<edge_visited_t, bool, EdgeCategory > EdgeVisited;
//typedef property<edge_category_t, std::string > EdgeCategory;


//BLACK - collinear and syntenic edge between segments, indegree==outdegree==1
//RED - collinear edge that traverses a syntenic breakpoint, degree!=1
//PURPLE - non-collinear edge indicating possible reversal, change in orientation
//GREEN - non-collinear edge indicative of some other flux
//BLUE - edge removed during a mincut to split an LCB
//ORANGERED - previously masked edge that was reintroduced during a merge
//YELLOW - new edge introduced during a mask short, merge adjacent iteration

enum EdgeCats { RED,PURPLE,CYAN,BLUE,GREEN,ORANGERED,YELLOW };
typedef property<edge_category_t, EdgeCats > EdgeCategory;
#if defined(STORE_EDGE_LABELS)
typedef property<edge_labelmask_t, GenomeSet, EdgeCategory > EdgeMask;
typedef property<edge_label_t, LabelMap, EdgeMask > EdgeProperties;
#else
typedef property<edge_labelmask_t, BitMask,EdgeCategory > EdgeProperties;
#endif

//TODO
//Determine if edge,vertex storage is better as vecS,listS,multisetS
typedef boost::adjacency_list < 
  vecS,               // Store out-edges of each vertex in a std::set
  vecS,               // Store vertex set in a std::vector 
  bidirectionalS,     // The file dependency graph is directed, support for in_edges as well as out_edges
  VertexProperties,   // vertex properties 
  EdgeProperties      // edge properties
  > Graph;

typedef Graph::vertex_descriptor Vertex;
typedef Graph::edge_descriptor Edge;

//Lookups
typedef std::map<std::string,Label> NameLabelMap; 
typedef std::map<Label, std::string> LabelNameMap; 
typedef std::map<VertexName, Vertex> NameVertexMap;
typedef std::map<Label,Label> SequenceGenomeMap;
//typedef std::map<pair<Vertex,Label>,Interval > VertexLabelIntervalMap; 
typedef boost::unordered_map<std::pair<Vertex,Label>,Interval > VertexLabelIntervalMap; 
typedef boost::unordered_map<Vertex,Interval > VertexIntervalMap;

typedef adjacency_list_traits < setS, vecS, directedS > LTraits;

#ifdef DEBUG
//Add edgecategory for printing
typedef adjacency_list < 
  listS, //to allow for removal of edges 
  vecS, 
  directedS,
  property < vertex_name_t, VertexName,
	     property < vertex_index_t, long,
			property < vertex_color_t, boost::default_color_type,
				   property < vertex_distance_t, long,
					      property < vertex_predecessor_t, LTraits::edge_descriptor > > > > >,    
  property < edge_capacity_t, long,
  property < edge_residual_capacity_t, long,
  property < edge_reverse_t, LTraits::edge_descriptor, 
  EdgeCategory  > > > > LGraph;
#else
typedef adjacency_list < 
  listS, //to allow for removal of edges 
  vecS, 
  directedS,
  property < vertex_name_t, VertexName,
	     property < vertex_index_t, long,
			property < vertex_color_t, boost::default_color_type,
				   property < vertex_distance_t, long,
					      property < vertex_predecessor_t, LTraits::edge_descriptor > > > > >,    
  property < edge_capacity_t, long,
	     property < edge_residual_capacity_t, long,
			property < edge_reverse_t, LTraits::edge_descriptor > > > > LGraph;
#endif

typedef LGraph::vertex_descriptor LVertex;
typedef LGraph::edge_descriptor EVertex;

typedef std::set<std::pair<Vertex,Vertex> > EdgeSet;
//typedef boost::unordered_set<std::pair<Vertex,Vertex> > EdgeSet;

typedef std::set<Vertex> VertexSet;
//typedef boost::unordered_set<Vertex> VertexSet;


struct iloc{
  int first;
  int second;
  int blocknum;
};

//Subs

void printtime(){
  time_t now;
  time(&now);
  struct tm *current = localtime(&now);
  current = localtime(&now);
  std::cerr << "TIME " << current->tm_hour << ":" << current->tm_min << ":" << current->tm_sec << std::endl;
}

class cutsdist
{
public:
  cutsdist(std::map<std::pair<LVertex,LVertex>,unsigned int > *m)
    :distmap(m)
  {}
  bool operator()( const std::pair<LVertex,LVertex>& v1, const std::pair<LVertex,LVertex>& v2 ) const {
    assert(distmap->find(v1)!=distmap->end());
    assert(distmap->find(v2)!=distmap->end());
    //.first is fmin coordinate
    if(distmap->find(v1)->second > distmap->find(v2)->second){
      return true;
    }
    else{
      return false;
    }
  }
  std::map<std::pair<LVertex,LVertex>,unsigned int > *distmap;
};

class coordsorder
{
public:
  coordsorder(VertexLabelIntervalMap *c, Label s)
    :coords(c),currentSeq(s)
  {}
  bool operator()( const LVertex& v1, const LVertex& v2 ) const {
    assert(coords->find(std::make_pair(v1,currentSeq))!=coords->end());
    assert(coords->find(std::make_pair(v2,currentSeq))!=coords->end());
    //.first is fmin coordinate
    if (coords->find(std::make_pair(v1,currentSeq))->second.first < coords->find(std::make_pair(v2,currentSeq))->second.first){
      return true;
    }
    else{
      return false;
    }
  }
  VertexLabelIntervalMap *coords;
  Label currentSeq;
};

class coordsorder_vertex
{
public:
  coordsorder_vertex(VertexIntervalMap *c)
    :coords(c)
  {}
  bool operator()( const LVertex& v1, const LVertex& v2 ) const {
    assert(coords->find(v1)!=coords->end());
    assert(coords->find(v2)!=coords->end());
    //.first is fmin coordinate
    if (coords->find(v1)->second.first < coords->find(v2)->second.first){
      return true;
    }
    else{
      return false;
    }
  }
  VertexIntervalMap *coords;
};

class matchmaporder
{
public:
  matchmaporder(std::map<Vertex,int> *c)
    :matchmap(c)
  {}
  bool operator()( const Vertex& v1, const Vertex& v2 ) const {
    assert(matchmap->find(v1)!=matchmap->end());
    assert(matchmap->find(v2)!=matchmap->end());
    if (matchmap->find(v1)->second > matchmap->find(v2)->second){
      return true;
    }
    else{
      return false;
    }
  }
  std::map<Vertex,int> *matchmap;
};

//################
//DEBUGGING UTILS
//
//
void printlabel(OrientedLabelSet & i){
  for(OrientedLabelSet::iterator j=i.begin();j!=i.end();++j){
    cerr << j->first << ":" << j->second << " ";
  }
}
void printlabel(OrientedLabelSet & i, LabelNameMap & index2sequence){
  for(OrientedLabelSet::iterator j=i.begin();j!=i.end();++j){
    cerr << index2sequence[j->first] << ":" << j->second << " ";
  }
}
//
//################################

template<class TGraph, class CoordMap, class SequenceGenomeMap>
void setedgemasks(TGraph & g, int distance, CoordMap & coordinates, SequenceGenomeMap & sequence2genome){
  typename property_map < TGraph, vertex_orient_t >::type orientmap = get(vertex_orient,g);
#if defined(STORE_EDGE_LABELS)
  typename property_map < TGraph, edge_label_t >::type elabelmap = get(edge_label,g);
#endif
  typename property_map < TGraph, edge_labelmask_t >::type elabelmaskmap = get(edge_labelmask,g);
  typename boost::graph_traits<TGraph>::edge_iterator eit,eit_end;
  eit_end = edges(g).second;
  for(eit = edges(g).first;eit!=eit_end;++eit){
      Edge e = *eit;
      Vertex sv = source(e,g);
      Vertex tv = target(e,g);

      //add extra edge labels
      OrientedLabelSet::iterator it_end = orientmap[sv].end();
      for(OrientedLabelSet::iterator it = orientmap[sv].begin();it != it_end;++it){
	Label seqidx = it->first;
	assert(sequence2genome.find(seqidx)!=sequence2genome.end());
	Label genomeidx = sequence2genome[seqidx];
	Orientation orient = it->second;
	Orientation rorient = (orient) ? false : true;
	//If sv--tv connected in label seqidx
	//Need to check original and reverse orientations of tv to ensure a match to sv
	if(orientmap[tv].find(*it) != orientmap[tv].end()
	   || orientmap[tv].find(std::make_pair(genomeidx,rorient)) != orientmap[tv].end()){
	  //vertices share label
	  //If edge label does not include seqidx, then we need to update
	  //if(elabelmap[*eit].find(genomeidx)==elabelmap[*eit].end()){
	  if(!elabelmaskmap[*eit].test(genomeidx)){
	    assert(!elabelmaskmap[*eit].test(genomeidx));
	    int dist=-1;
	    //Need to check if coordinates for seqidx
	    if(coordinates.find(std::make_pair(source(*eit,g),seqidx))!=coordinates.end()
	       && coordinates.find(std::make_pair(target(*eit,g),seqidx))!=coordinates.end()){
	      if(it->second==true){
		if(coordinates[std::make_pair(target(*eit,g),seqidx)].first >= coordinates[std::make_pair(source(*eit,g),seqidx)].second){
		  dist = coordinates[std::make_pair(target(*eit,g),seqidx)].first - coordinates[std::make_pair(source(*eit,g),seqidx)].second;
		}
		else{
		  dist = coordinates[std::make_pair(source(*eit,g),seqidx)].first - coordinates[std::make_pair(target(*eit,g),seqidx)].second;
		}
		//assert(dist>=0);
		if(dist<0){
		  //std::cout << source(e,g) << "-" << target(e,g) << " " 
		  //<< coordinates[std::make_pair(source(*eit,g),seqidx)].first << "-" << coordinates[std::make_pair(source(*eit,g),seqidx)].second
		  //<< " " 
		  //<< coordinates[std::make_pair(target(*eit,g),seqidx)].first << "-" << coordinates[std::make_pair(target(*eit,g),seqidx)].second
		  //<< " " << dist << ":COORDS " << std::endl;
		  dist=0;
		}
	      }
	      else{
		if(coordinates[std::make_pair(target(*eit,g),seqidx)].first >= coordinates[std::make_pair(source(*eit,g),seqidx)].second){
		  dist = coordinates[std::make_pair(target(*eit,g),seqidx)].first - coordinates[std::make_pair(source(*eit,g),seqidx)].second;
		}
		else{
		  dist = coordinates[std::make_pair(source(*eit,g),seqidx)].first - coordinates[std::make_pair(target(*eit,g),seqidx)].second;
		}
		if(dist<0){
		  //std::cout << source(e,g) << "-" << target(e,g) << " " 
		  //<< coordinates[std::make_pair(source(*eit,g),seqidx)].first << "-" << coordinates[std::make_pair(source(*eit,g),seqidx)].second
		  //<< " " 
		  //<< coordinates[std::make_pair(target(*eit,g),seqidx)].first << "-" << coordinates[std::make_pair(target(*eit,g),seqidx)].second
		  //<< " " << dist << ":COORDS" << std::endl;
		  dist=0;
		}
		//assert(dist>=0);
	      }
	      assert(dist>=0);
	      if(dist<=distance){
#if defined(STORE_EDGE_LABELS)
		elabelmap[*eit].insert(std::make_pair(genomeidx,dist));
#endif
		elabelmaskmap[*eit].set(genomeidx,1);
	      }
	      else{
		//elabelmaskmap[*eit].set(genomeidx,0);
	      }
	    }
	    else{

	    }
	  }
	}
      }
#if defined(STORE_EDGE_LABELS)
      BitMask mask;
      unsigned int numset=0;
      for(LabelMap::iterator i1 = labelmap[*eit].begin();i1 != labelmap[*eit].end();++i1){
	if(i1->second<=distance){
	  assert(i1->first>=0&&i1->first<BITMAX);
	  mask.set(i1->first,1);
	  numset++;
	}
      }
      assert(numset==labelmap[*eit].size());
      assert(mask.any());
      if(mask!=elabelmaskmap[*eit]){
	std::cerr << "Mask     " << mask << std::endl;
	std::cerr << "EdgeMask " << elabelmaskmap[*eit] << std::endl;
      }
      assert(mask==elabelmaskmap[*eit]);
      //put(edge_labelmask,g,*eit,mask); 
#endif
  }
}

//Sets the following graph vertex properties: vertex_vlabelmask, vertex_orientmask
template<class TGraph>
void setvertexmasks(TGraph & g, SequenceGenomeMap & sequence2genome){
  typename property_map < TGraph, vertex_orient_t >::type vmap = get(vertex_orient,g);
  typename property_map < TGraph, vertex_vlabelmask_t >::type lmaskmap = get(vertex_vlabelmask,g);
  typename property_map < TGraph, vertex_orientmask_t >::type omaskmap = get(vertex_orientmask,g);

  typename boost::graph_traits<TGraph>::vertex_iterator vit_end = vertices(g).second;
  for(typename boost::graph_traits<TGraph>::vertex_iterator 
	vit = vertices(g).first;vit!=vit_end;++vit){
    Vertex v = *vit;
    OrientedLabelSet::iterator o_end = vmap[v].end();
    for(OrientedLabelSet::iterator o=vmap[v].begin();o!=o_end;++o){
      Label seqidx = o->first;
      Label genomeidx = sequence2genome[seqidx];
      assert(genomeidx>=0&&genomeidx<BITMAX);
      //set lmask for each genome
      lmaskmap[v].set(genomeidx,1);
      //set omask for all + orientation
      if(o->second == true){//+ orientation
	omaskmap[v].set(genomeidx,1);
      }
    }
#ifdef DEBUG
    cerr << "VERTEX: " << get(vertex_name,g,v) << endl;
    cerr << "VERTEXIDX: " << v << endl;
    cerr << "LMASK :" << lmaskmap[v]  << endl;
    cerr << "OMASK :" << omaskmap[v] << endl;
#endif
  }
}

void updateCoordinates(VertexLabelIntervalMap & coordinates,
		       SequenceGenomeMap & sequence2genome){
  
  std::map<Label,Coordinate> maxcoord;
  for(VertexLabelIntervalMap::iterator it=coordinates.begin();it!=coordinates.end();it++){
    Label seqidx = it->first.second;
    assert(it->second.first<=it->second.second);
    maxcoord[seqidx] = (it->second.second > maxcoord[seqidx]) ? it->second.second : maxcoord[seqidx];
  }
  std::map<Label,std::vector<Label> > genome2sequence;
  for(SequenceGenomeMap::iterator sit=sequence2genome.begin();sit!=sequence2genome.end();sit++){
    genome2sequence[sit->second].push_back(sit->first);
  }
  std::map<Label,Coordinate> seqoffset;
  for(std::map<Label,std::vector<Label> >::iterator git=genome2sequence.begin();git!=genome2sequence.end();git++){
    //Label genomeidx = git->first;
    Coordinate curroffset=0;
    for(std::vector<Label>::iterator sit=git->second.begin();sit!=git->second.end();sit++){
      Label seqidx = *sit;
      seqoffset[seqidx] = curroffset;
      curroffset = curroffset+maxcoord[seqidx];
    }
  }
  VertexLabelIntervalMap newcoordinates;
  for(VertexLabelIntervalMap::iterator it=coordinates.begin();it!=coordinates.end();it++){
    Label seqidx = it->first.second;
    Vertex v = it->first.first;
    Label genomeidx = sequence2genome[seqidx];
    Coordinate newbeg=seqoffset[genomeidx]+it->second.first;
    Coordinate newend=seqoffset[genomeidx]+it->second.second;
    //update coordinate map
    pair<Vertex,Label> key = std::make_pair(v,seqidx);
    pair<Vertex,Label> value = std::make_pair(newbeg,newend);
    //assert(newcoordinates.find(key)==newcoordinates.end());
    //if(newcoordinates.find(key)!=newcoordinates.end()){
    //std::cerr << "Duplicate V:"<<v<<" seqidx:" <<seqidx << " genomeidx:"<<genomeidx << " " <<newbeg << "-" << newend << std::endl;
    //}
    //else{
    //std::cerr << "Storing V:"<<v<<" seqidx:" <<seqidx << " genomeidx:"<<genomeidx << " " <<newbeg << "-" << newend << std::endl;
    //}
    newcoordinates[key]=value;
    assert(newcoordinates[key].first==newbeg);
    assert(newcoordinates[key].second==newend);
  }
  coordinates=newcoordinates;
}

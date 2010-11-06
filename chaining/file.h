

//################################
//File IO
/*Block format is 6 column
anchor seqindex genomeindex orient beg end
*/ 
void read_blocks(std::istream &in, 
		 Graph & g, 
		 NameVertexMap & name2vertex, 
		 NameLabelMap & genome2index, 
		 NameLabelMap & sequence2index, 
		 VertexLabelIntervalMap & coordinates,
		 int distance){

  NameLabelMap::iterator pos1; 
  NameVertexMap::iterator pos;
  VertexLabelIntervalMap::iterator pos2;  
  bool inserted; 
  Label seqindex=0;
  Vertex news;
  Edge e1,newe;
  int edges=0;
  bool found;

  std::string line;
  typedef tokenizer<char_separator<char> > Tok;
  int field=0;
  VertexName sname=0;
  std::string sequence,genome;
  Orientation sorient=false;
  Coordinate sbeg=0,send=0;
  int dist=0;
  std::string sorientstr;
  OrientedLabelSet newso;
  
  property_map < Graph, vertex_orient_t >::type orientmap = get(vertex_orient,g);
  property_map < Graph, vertex_len_t >::type lenmap = get(vertex_len,g);
#if defined(STORE_EDGE_LABELS)
  property_map < Graph, edge_label_t >::type labelmap = get(edge_label,g);
#endif
  property_map < Graph, edge_labelmask_t >::type elabelmaskmap = get(edge_labelmask,g);

  
  vector<int> ordercounts(BITMAX);

  while (getline(in, line)) {
    Tok tok(line, char_separator<char>(" "));
    field=0;
    for (Tok::iterator id = tok.begin(); id != tok.end(); ++id) {
      switch(field){
      case 0:
	try{
	  sname  = lexical_cast<VertexName>(*id);
	  field++;
	}
	catch (std::exception e){
	}
	break;
      case 1:
	try{
	  sequence = lexical_cast<std::string>(*id);
	  field++;
	}
	catch (std::exception e){
	}
	break;
      case 2:
	try{
	  genome = lexical_cast<std::string>(*id);
	  field++;
	}
	catch (std::exception e){
	}
	break;
      case 3:
	try {
	  sorientstr = lexical_cast<std::string>(*id);
	  if(sorientstr == "+"){
	    sorient = true;
	  }
	  else{
	    assert(sorientstr=="-");
	    sorient = false;
	  }
	  field++;
	}
	catch (std::exception e){
	}
	break;
      case 4:
	try {
	  sbeg = lexical_cast<Coordinate>(*id);
	  field++;
	}
	catch (std::exception e){
	}
	break;
      case 5:
	try {
	  send = lexical_cast<Coordinate>(*id);
	  field++;
	}
	catch (std::exception e){
	}
	break;
      }
    //6 column table
    //7 column table includes orientations
    //11 column table includes coordinates
      if(field==6){
	//Set either returns existing index or inserts
	tie(pos1, inserted) = sequence2index.insert(std::make_pair(sequence, 0));
	if (inserted) {
	  pos1->second = sequence2index.size();
	  seqindex = pos1->second;
	  assert(seqindex>=0&&seqindex<BITMAX);
	}
	else{
	  seqindex = pos1->second;
	  assert(seqindex>=0&&seqindex<BITMAX);
	}
	
	//Add source vertex
	//Set either returns existing vertex or inserts
	tie(pos, inserted) = name2vertex.insert(std::make_pair(sname, Vertex()));
	if (inserted) {
	  news = add_vertex(VertexProperties(sname),g);
	  pos->second = news;
	} else{
	  news = pos->second;
	}
	//Add oriented label
	orientmap[news].insert(make_pair(seqindex,sorient));
      
	//Save coordinates
	assert(sbeg!=send);
	tie(pos2, inserted) = coordinates.insert(std::make_pair(make_pair(news,seqindex), std::make_pair(sbeg,send)));
	if(lenmap[news]>0){
	  assert(lenmap[news]==pos2->second.second-pos2->second.first);
	}
	lenmap[news] = pos2->second.second-pos2->second.first;
      } 
      else{
	//Ignoring line
      }
    }
    //iterate over all seqs
    NameLabelMap::iterator sit,sit_end;
    sit_end = sequence2index.end();
    for(sit = sequence2index.begin();sit!=sit_end;++sit){
      Label seqidx = sit->second;
      //store all vertices with this seqlabel
      list<Vertex> sortedV;
      boost::graph_traits<Graph>::vertex_iterator vit,vit_end;
      vit_end = vertices(g).second;
      VertexIntervalMap currcoords;
      for(vit=vertices(g).first;vit!=vit_end;++vit){
	VertexLabelIntervalMap::iterator cit = coordinates.find(std::make_pair(*vit,seqidx));
	if(cit!=coordinates.end()){
	  sortedV.push_back(*vit);
	  currcoords.insert(std::make_pair(*vit,
					   cit->second));
	}
      }
      //sort
      sortedV.sort(coordsorder_vertex(&currcoords));
      Vertex currvertex,prevvertex;
      list<Vertex>::iterator it,it_end;
      for(it=sortedV.begin();it!=sortedV.end();++it){
	currvertex = *it;
	if(it==sortedV.begin()){
	  prevvertex=currvertex;
	}
	else{
	  dist = abs(coordinates[std::make_pair(prevvertex,seqidx)].second - coordinates[std::make_pair(currvertex,seqidx)].first);
	  //Add edge if
	  assert(dist>=0);
	  if(dist <= distance){
	    tie(e1,found) = edge(prevvertex,currvertex, g);
	    if(found){
	      //existing edge between prevvertex--currvertex
	      //add attributes from Graph g edge,e to Graph gcomp edge,e1
#if defined(STORE_EDGE_LABELS)
	      labelmap[e1].insert(std::make_pair(seqindex,dist));
#endif
	      elabelmaskmap[e1].set(seqindex,1);
	    }
	    else{
	      //Code to handle directed graph where
	      //reverse orientation
	      //TODO
	      //Need to consider case where
	      //this edge is mis-oriented introducing an artificial breakpoint
	      //in the chain
	      tie(e1,found) = edge(currvertex,prevvertex,g);
	      if(found){
#if defined(STORE_EDGE_LABELS)
		labelmap[e1].insert(std::make_pair(seqindex,dist));
#endif
		elabelmaskmap[e1].set(seqindex,1);
	      }
	      else{
		bool inserted;
		Edge e1;
#if defined(STORE_EDGE_LABELS)
		LabelMap labels;
		labels[seqindex] = dist;
		tie(e1, inserted) = add_edge(prevvertex,currvertex,EdgeProperties(labels),g);
#else
		tie(e1, inserted) = add_edge(prevvertex,currvertex,EdgeProperties(),g);
#endif
		assert(inserted);
		elabelmaskmap[e1].set(seqindex,1);
	      }
	    }
	    edges++;
	  }
	}
      }
    }
  }
}


/*
Projection input is 
anchor1 anchor2 seqindex dist genomeindex orient1 orient2 beg1 end1 beg2 end2
eg
0 1 0 0 0 + + 0 196 196 15348
1 3 0 1 0 + + 196 15348 15349 20373
*/

void read_pairwiseprojection(std::istream &in, 
			     Graph & g, 
			     NameVertexMap & name2vertex, 
			     NameLabelMap & genome2index, 
			     NameLabelMap & sequence2index, 
			     VertexLabelIntervalMap & coordinates,
			     SequenceGenomeMap & sequence2genome,
			     int distance,
			     int minanchor){

  NameVertexMap::iterator pos;
  NameLabelMap::iterator pos1; 
  VertexLabelIntervalMap::iterator pos2;  
  SequenceGenomeMap::iterator pos3;

  bool inserted; 
  Label seqindex,genomeindex;
  Vertex news, newt;
  Edge e1,newe;
  int edges=0;
  bool found;

  std::string line;
  typedef tokenizer<char_separator<char> > Tok;
  int field=0;
  VertexName sname=0,tname=0;
  std::string sequence,genome;
  Orientation sorient=false,torient=false;
  Coordinate sbeg=0,send=0,tbeg=0,tend=0;
  int dist=0;
  std::string sorientstr;
  std::string torientstr;

  OrientedLabelSet newso,newto;
  
  property_map < Graph, vertex_orient_t >::type orientmap = get(vertex_orient,g);
  property_map < Graph, vertex_label_t >::type vlabelmap = get(vertex_label,g);
  property_map < Graph, vertex_genome_t >::type genomemap = get(vertex_genome,g);
  property_map < Graph, vertex_len_t >::type lenmap = get(vertex_len,g);

  property_map < Graph, edge_labelmask_t >::type elabelmaskmap = get(edge_labelmask,g);
#if defined(STORE_EDGE_LABELS)
  property_map < Graph, edge_label_t >::type labelmap = get(edge_label,g);
#endif


  
  vector<int> ordercounts(BITMAX);

  while (getline(in, line)) {
    //std::cerr << line << std::endl;
    Tok tok(line, char_separator<char>(" "));
    field=0;
    for (Tok::iterator id = tok.begin(); id != tok.end(); ++id) {
      switch(field){
      case 0:
	try{
	  sname  = lexical_cast<VertexName>(*id);
	  field++;
	}
	catch (std::exception e){
	}
	break;
      case 1:
	try{
	  tname  = lexical_cast<VertexName>(*id);
	  field++;
	}
	catch (std::exception e){
	}
	break;
      case 2:
	try{
	  sequence = lexical_cast<std::string>(*id);
	  field++;
	}
	catch (std::exception e){
	}
	break;
      case 3:
	try{
	  dist = lexical_cast<long int>(*id);
	  field++;
	}
	catch (std::exception e){
	}
	break;
      case 4:
	try{
	  genome = lexical_cast<std::string>(*id);
	  field++;
	}
	catch (std::exception e){
	}
	break;
      case 5:
	try {
	  sorientstr = lexical_cast<std::string>(*id);
	  if(sorientstr == "+"){
	    sorient = true;
	  }
	  else{
	    assert(sorientstr=="-");
	    sorient = false;
	  }
	  field++;
	}
	catch (std::exception e){
	}
	break;
      case 6:
	try {
	  torientstr = lexical_cast<std::string>(*id);
	  if(torientstr == "+"){
	    torient = true;
	  }
	  else{
	    assert(torientstr=="-");
	    torient = false;
	  }
	  field++;
	}
	catch (std::exception e){
	}
	break;
      case 7:
	try {
	  sbeg = lexical_cast<Coordinate>(*id);
	  field++;
	}
	catch (std::exception e){
	}
	break;
      case 8:
	try {
	  send = lexical_cast<Coordinate>(*id);
	  field++;
	}
	catch (std::exception e){
	}
	break;
      case 9:
	try {
	  tbeg = lexical_cast<Coordinate>(*id);
	  field++;
	}
	catch (std::exception e){
	}
	break;
      case 10:
	try {
	  tend = lexical_cast<Coordinate>(*id);
	  field++;
	}
	catch (std::exception e){
	}
	break;
      }
    }
    //5 column table is minimum input
    //7 column table includes orientations
    //11 column table includes coordinates
    if(field==5 || field==7 || field==11){
      if(field==5){
	cerr << "Incomplete file "<< endl;
	sorient = true;
	torient = true;
      }
      if(abs(tend-tbeg)>=minanchor && abs(send-sbeg)>=minanchor){
	//Set either returns existing index or inserts
	tie(pos1, inserted) = sequence2index.insert(std::make_pair(sequence, 0));
	if (inserted) {
	  pos1->second = sequence2index.size();
	  seqindex = pos1->second;
	  //assert(seqindex>=0&&seqindex<BITMAX);
	}
	else{
	  seqindex = pos1->second;
	  //assert(seqindex>=0&&seqindex<BITMAX);
	}
	
	//Add source vertex
	//Set either returns existing vertex or inserts
	tie(pos, inserted) = name2vertex.insert(std::make_pair(sname, Vertex()));
	if (inserted) {
	  news = add_vertex(VertexProperties(sname),g);
	  pos->second = news;
	} else{
	  news = pos->second;
	}
	
	//Add genome
	tie(pos1, inserted) = genome2index.insert(std::make_pair(genome, 0));
	if (inserted) {
	  pos1->second = genome2index.size();
	  genomeindex = pos1->second;
	  assert(genomeindex>=0&&genomeindex<BITMAX);
	}
	else{
	  genomeindex = pos1->second;
	  assert(genomeindex>=0&&genomeindex<BITMAX);
	}
	
	//Save sequence2genome lookup
	//TODO
	//Use genomeindex
	//TEMP HACK to support 1 seq per genome
	//genomeindex=seqindex;
	//</HACK>
	
	genomemap[news].insert(genomeindex);
	
	tie(pos3, inserted) = sequence2genome.insert(std::make_pair(seqindex, genomeindex));
	if (inserted) {
	}
	else{
	  assert(seqindex==pos3->first);
	  assert(genomeindex=pos3->second);
	}
	
	//Add oriented label
	orientmap[news].insert(make_pair(seqindex,sorient));
	vlabelmap[news].insert(seqindex);
	
	//Add target vertex
	tie(pos, inserted) = name2vertex.insert(std::make_pair(tname, Vertex()));
	if (inserted) {
	  newt = add_vertex(VertexProperties(tname),g);
	  pos->second = newt;
	} else{
	  newt = pos->second;
	}
	
	genomemap[newt].insert(genomeindex);
	orientmap[newt].insert(make_pair(seqindex,torient));
	vlabelmap[newt].insert(seqindex);
	//Save coordinates
	if(field==11){
	  assert(sbeg!=send);
	  assert(tbeg!=tend);
	  tie(pos2, inserted) = coordinates.insert(std::make_pair(make_pair(news,seqindex), std::make_pair(sbeg,send)));
	  //std::cerr << "Coords source V:" << news << " " << sbeg << "-" << send << std::endl;
	  if(lenmap[news]>0){
	    //assert(lenmap[news]==pos2->second.second-pos2->second.first);
	    lenmap[news] = (lenmap[news]>pos2->second.second-pos2->second.first) ? lenmap[news] : pos2->second.second-pos2->second.first;
	  }
	  else{
	    lenmap[news] = pos2->second.second-pos2->second.first;
	  }
	  tie(pos2, inserted) = coordinates.insert(std::make_pair(make_pair(newt,seqindex), std::make_pair(tbeg,tend)));
	  //std::cerr << "Coords target V:" << newt << " " << tbeg << "-" << tend << std::endl;
	  if(lenmap[newt]>0){
	    //assert(lenmap[newt]==pos2->second.second-pos2->second.first);
	    lenmap[newt] = (lenmap[newt]>pos2->second.second-pos2->second.first) ? lenmap[newt] : pos2->second.second-pos2->second.first;
	  }
	  else{
	    lenmap[newt] = pos2->second.second-pos2->second.first;
	  }
	}
	//Add edge if
	assert(dist>=0);
	if(dist <= distance){
	  tie(e1,found) = edge(news,newt, g);
	  if(found){
	    //existing edge between news--newt
	    //add attributes from Graph g edge,e to Graph gcomp edge,e1
#if defined(STORE_EDGE_LABELS)
	    labelmap[e1].insert(std::make_pair(genomeindex,dist));
#endif
	    elabelmaskmap[e1].set(genomeindex,1);
	  }
	  else{
	    //Code to handle directed graph where
	    //reverse orientation
	    //TODO
	    //Need to consider case where
	    //this edge is mis-oriented introducing an artificial breakpoint
	    //in the chain
	    tie(e1,found) = edge(newt,news,g);
	    if(found){
#if defined(STORE_EDGE_LABELS)
	      labelmap[e1].insert(std::make_pair(genomeindex,dist));
#endif
	      elabelmaskmap[e1].set(genomeindex,1);
	    }
	    else{
	      bool inserted;
	      Edge e1;
#if defined(STORE_EDGE_LABELS)
	      LabelMap labels;
	      labels[genomeindex] = dist;
	      tie(e1, inserted) = add_edge(news,newt,EdgeProperties(labels),g);
#else
	      tie(e1, inserted) = add_edge(news,newt,EdgeProperties(),g);
#endif
	      assert(inserted);
	      elabelmaskmap[e1].set(genomeindex,1);
	    }
#ifdef DEBUG
	std::cerr << "Added edge " << sname << " " << tname << std::endl;
#endif
	  }
	  edges++;
	}
	else{
#ifdef DEBUG
	std::cerr << "Skipping edge dist>distance " << line << std::endl;
#endif
	}
      }
      else{
#ifdef DEBUG
	std::cerr << "Skipping short anchor " << line << std::endl;
#endif
	//Skipping short anchor
      }
    }
    else{
#ifdef DEBUG
	std::cerr << "Ignoring line " << line << std::endl;
#endif
      //Ignoring line
    }
  }
}
//Add this to the dot file to force drawing of labels in large graphs
//node [fontsize="9",margin="0.0,0.0",fixedsize=true];

template<typename TGraph, typename Tedgelabelmap> 
  void do_write_graphviz(TGraph &g, 
			 std::string fname, 
			 std::vector<int> & cc, 
			 VertexLabelIntervalMap & coordinates, 
			 EdgeSet & maskedEdges,
			 VertexSet & maskedLCBs,
			 Tedgelabelmap & edgelabelmap,
			 bool expandlabel){

  typedef typename TGraph::vertex_descriptor TVertex;
  typedef typename TGraph::edge_descriptor TEdge;
  //property_map < Graph,edge_stringname_t >::type edgelabelmap;// = get(edge_stringname, g);
  typename property_map < TGraph,edge_category_t >::type  ecatmap = get(edge_category,g);

  //Set up dynamic properties for graphviz
  boost::dynamic_properties dp;
  dp.property("id", get(vertex_name, g));


  std::map<TEdge,std::string> edgecatmap;
  //Need to set edge, vertex labels
  //Build label map
  std::map<TVertex, std::string> vertexlabelmap;
  std::map<TEdge, std::string> efmap;
  std::map<TVertex, std::string> vfmap;
  std::map<TVertex, std::string> shapemap;
  std::map<TEdge, std::string> linemap;
  //std::map< Graph::edge_descriptor, std::string> edgelabelmap;
  for(typename boost::graph_traits<TGraph>::vertex_iterator 
	vit = vertices(g).first;vit!=vertices(g).second;++vit){
    Vertex v = *vit;
    std::ostringstream labelstring;
#ifdef PRINTSEQS
    labelstring << get(vertex_name,g,v) << " " 
		<< v
		<< " "
		<< "CC" << cc[v];
    //if(get(vertex_relorder,g,v)){
    //labelstring << " TL" << get(vertex_relorder,g,v);
    //}
    labelstring << "\\n";
    if(expandlabel){
      OrientedLabelSet olabel = get(vertex_orient,g,v);
      for(OrientedLabelSet::iterator it = olabel.begin();it!=olabel.end();++it){
	//TODO support genomeidx
	labelstring << "S"  << it->first << ":" << (it->second ? '+' : '-') 
		    << ":" 
		    << coordinates[std::make_pair(v,it->first)].first << "-" << coordinates[std::make_pair(v,it->first)].second << "\\n";
      }
      
    }
#else
    labelstring << v << "\\n";
#endif
    vertexlabelmap[v]=labelstring.str();
    vfmap[v]="6";
    if(maskedLCBs.find(v)!=maskedLCBs.end()){
      shapemap[v]="diamond";
    }
    else{
      shapemap[v]="circle";
    }
  } 
  for(typename boost::graph_traits<TGraph>::edge_iterator 
	eit = edges(g).first;eit!=edges(g).second;++eit){
    TEdge e = *eit;
    std::ostringstream labelstring;
    unsigned int numset=0;
#if defined(STORE_EDGE_LABELS)
    LabelMap currlm = get(edge_label,g,e);
    //labelstring << "MASKS:" << get(edge_labelmask,g,e) << "\\n";
    for(LabelMap::iterator it = currlm.begin(); it!=currlm.end(); ++it){
      labelstring << it->first << ":" << it->second << " ";
      numset++;
    }
    assert(numset==currlm.size());
#endif
    if(maskedEdges.find(std::make_pair(source(e,g),target(e,g)))!=maskedEdges.end() 
       || maskedEdges.find(std::make_pair(target(e,g),source(e,g)))!=maskedEdges.end()){
      //only true if g is of type filteredgraph assert(false);
      linemap[e] = "dashed,bold";
    }
    else{
      linemap[e] = "solid";
    }
    //edgelabelmap[e] = labelstring.str();
    efmap[e]="6";
    switch (ecatmap[e]){
    case RED: //default
      edgecatmap[e] = "red";
      break;
    case GREEN:
      edgecatmap[e] = "green";
      break;
    case BLUE: //cut by mincut
      edgecatmap[e] = "blue";
      break;
    case ORANGERED: //introduced via a merge
      edgecatmap[e] = "yellow";
      break;
    case PURPLE: //change in relative orientation
      edgecatmap[e] = "purple";
      break;
    case CYAN:
      edgecatmap[e] = "cyan";
      break;
    default:
      assert(false);
      break;
    }
  }
  
  boost::associative_property_map< std::map<TVertex, std::string> >
    vlabel_map(vertexlabelmap);
  //boost::associative_property_map< std::map<Vertex, std::string> >
  //vfontmap(vfmap);
  //boost::associative_property_map< std::map<Edge, std::string> >
  //efontmap(efmap);
  boost::associative_property_map< std::map<TVertex, std::string> >
    bshapemap(shapemap);
  boost::associative_property_map< std::map<TEdge, std::string> >
    blinemap(linemap);
  //boost::associative_property_map< std::map<TEdge, Tedgelabel> >
  //elmap(edgelabelmap);
  boost::associative_property_map< std::map<TEdge, std::string> >
    ecmap(edgecatmap);
 
  dp.property("label",vlabel_map);
  //dp.property("label",edgelabelmap);
  dp.property("label",edgelabelmap);
  dp.property("shape",bshapemap);
  dp.property("style",blinemap);
  //dp.property("fontsize",vfontmap);
  //dp.property("fontsize",efontmap);
  dp.property("color",ecmap);
  dp.property("color",ecmap);

  //dp.property("rankdir","LR");
  //dp.property("rotate","90");
  //Open file
  std::ofstream gout;
  gout.open(fname.c_str());
  std::string node_id("id");

  std::map<std::string,std::string> graph_attr, vertex_attr, edge_attr;
  graph_attr["rankdir"] = "LR";
  graph_attr["rotate"] = 90;
  write_graphviz(gout, g, 
		 dynamic_vertex_properties_writer(dp,node_id),
		 dynamic_properties_writer(dp),
		 make_graph_attributes_writer(graph_attr,vertex_attr,edge_attr));
		 //graph::detail::node_id_property_map<Vertex>(dp,node_id));

  gout.close();
}


template<typename TGraph> 
void do_write_graphviz(TGraph &g, 
		       std::string fname, 
		       std::vector<int> & cc, 
		       VertexLabelIntervalMap & coordinates, 
		       EdgeSet & maskedEdges,
		       VertexSet & maskedLCBs,
		       bool expandlabel=true){
  std::map<Edge,std::string> nullmap;
  boost::associative_property_map< std::map<Edge,std::string > > edgelabelmap(nullmap);
  
  do_write_graphviz(g,fname,cc,coordinates,maskedEdges,maskedLCBs,edgelabelmap,expandlabel);

}

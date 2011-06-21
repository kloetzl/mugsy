//##########
//breakLCBmincutconnect()
//Interpret anchor graph as a flow network
//Use mincut,max-flow to partition the graph to fullfill criteria
//1)gaps <= distance
//2)no conflicting orientations
//3)at most one sequence per genome (important for draft data)
//Cut edges are tagged as BLUE in the input graph. 

template<typename TGraph, typename TGraph2>
inline
int breakLCBmincutconnect(std::vector<LCB > &componentMap,
			  std::vector<int> &ccvmap,
			  EdgeSet&maskedEdges,
			  TGraph g,
			  TGraph2 fglcbsyn,
			  unsigned int distance,
			  VertexLabelIntervalMap &coordinates,
			  std::set<Label> &seqidxSet,
			  NameVertexMap &name2vertex,
			  SequenceGenomeMap & sequence2genome,
			  int filenumoffset=0){
  bool found=false;
  int lcbcount=0;

  bool reusesupernodes=false;
  //int SEARCH_RADIUS=std::numeric_limits<unsigned int>::max();
  int DEFAULT_CAP=1;
  int numcuts=0;
  int MINSPANLEN=0;
  int cutcount=0;
  //Determine cuts over each LCB
  for(std::vector<LCB >::iterator it = componentMap.begin();it!=componentMap.end();++it){
#ifdef DEBUG
      std::cerr << "mincut LCB count:" << lcbcount << " num vertices " << it->size() << std::endl;
#endif
    lcbcount++;
    std::map<VertexName, LVertex> currlcbv;
    std::map<VertexName, LVertex>::iterator pos;
    std::map<LVertex,Vertex> vmap;
    typename property_map < TGraph, vertex_orientmask_t >::type orientmaskmap = get(vertex_orientmask,g);
    typename property_map < TGraph, vertex_vlabelmask_t >::type labelmaskmap = get(vertex_vlabelmask,g);
    typename property_map < TGraph, vertex_orient_t>::type orientmap = get(vertex_orient,g);

    std::map<std::pair<VertexName,VertexName>,std::pair<LVertex,LVertex> > cuts;
    std::map<std::pair<LVertex,LVertex>,std::pair<VertexName,VertexName> > revcuts;
    std::set<std::pair<LVertex,LVertex> > cutsnodes;
    std::map<std::pair<LVertex,LVertex>,unsigned int> cutsdistmap;
    std::map<LVertex,int> supernodes;
    //std::set<LVertex> supernodes2;
    int snodeedges=0;


    //The connectivity graph used for max flow, min cut
    LGraph currlcbg;
    int supercount=0;
    property_map < LGraph, edge_capacity_t >::type
      capacity = get(edge_capacity, currlcbg);
    property_map < LGraph, edge_reverse_t >::type 
      rev = get(edge_reverse, currlcbg);
    property_map < LGraph, edge_residual_capacity_t >::type
      residual_capacity = get(edge_residual_capacity, currlcbg);
    
    //The set of vertices for the LCB
    LCB blockV = *it;
    //Get sequence labels for this lcb
    //Slows us down for complete genomes, but helps us with perf for draft genomes
    std::set<Label> currseqidxSet;
    typename property_map < TGraph, vertex_label_t >::type vlabelmap = get(vertex_label,g);
    std::map<Label,std::set<Label> > seqspergenomeMap; //tracks the number of seqs per genome in an LCB
    std::map<Label,std::set<Label> >::iterator gpos;
    std::set<Label>::iterator spos,spos2;
    std::map<Label,std::set<LVertex> > seqsvertex;
    std::set<LVertex>::iterator vpos,vpos2;
    bool inserted;
    for(LCB::iterator vit = blockV.begin();vit!=blockV.end();++vit){
      for(LabelSet::iterator sit = vlabelmap[*vit].begin();sit!=vlabelmap[*vit].end();++sit){
	currseqidxSet.insert(*sit);
	tie(gpos, inserted) = seqspergenomeMap.insert(std::make_pair(sequence2genome[*sit],std::set<Label>()));
	gpos->second.insert(*sit);
      }
    }
    std::vector<LGraph::edge_descriptor> disconnecting_set;
    for(std::set<Label>::iterator it2 = currseqidxSet.begin(); it2 != currseqidxSet.end(); ++it2){
      assert(seqidxSet.find(*it2) != seqidxSet.end());
      Label seqidx = *it2;
      //Label genomeidx = sequence2genome[seqidx];
      std::list<LVertex> sortedV;
      unsigned int spanlen=0;
      //need to create custom coords map for LGraph
      VertexIntervalMap currcoords;

      //Create special graph currlcbv to represent the current lcb
      //createLCBGraph(g,currlcbv,blockV);
      
      for(LCB::iterator vit = blockV.begin();vit!=blockV.end();++vit){
	Vertex v=*vit;
	VertexName sname = get(vertex_name,g,v);
	LVertex news;
	//
	//Insert vertex into currlcbg if needed
	tie(pos, inserted) = currlcbv.insert(std::make_pair(sname, LVertex()));
	if(inserted){
	  news = add_vertex(sname,currlcbg);
	  currlcbv[sname]=news;
	  vmap[news]=v;
	}
	else{
	  news = pos->second;
	}

	//
	//Save coordinate information for news
	if(coordinates.find(std::make_pair(v,seqidx))!=coordinates.end()){
	  assert(coordinates.find(std::make_pair(v,seqidx))->second.first<coordinates.find(std::make_pair(v,seqidx))->second.second);
	  sortedV.push_back(news);
	  spanlen = spanlen + get(vertex_len,g,v);
#ifdef DEBUG
	  std::cerr << "seqidx: " << seqidx << " V:" << get(vertex_name,g,v) 
		    << " len:" << get(vertex_len,g,news) << " spanlen:" << spanlen << std::endl;
#endif
	  currcoords.insert(std::make_pair(news,
					   coordinates.find(std::make_pair(v,seqidx))->second));
	  seqsvertex[seqidx].insert(news);
	}

	//
	//Add all edges for news
	//First make sure target vertex is part of currlcbg
	//graph_traits<LCBSynFilterGraph>::out_edge_iterator out_i, out_end;
	typename graph_traits<TGraph2>::out_edge_iterator out_i, out_end;
	for(tie(out_i, out_end) = out_edges(v, fglcbsyn); out_i != out_end; ++out_i){
#ifdef CUTLCBEDGESONLY
	  if(ccvmap[v] != ccvmap[target(*out_i,g)]){
	    std::cerr << "Skipping edge, outside LCB" << std::endl;
	    continue;
	  }
#endif
	  VertexName tname = get(vertex_name,g,target(*out_i,g));
	  LVertex newt;
	  tie(pos, inserted) = currlcbv.insert(std::make_pair(tname, LVertex()));
	  if(inserted){
	    newt = add_vertex(tname,currlcbg);
	    currlcbv[tname]=newt;
	    vmap[newt]=target(*out_i,g);
	  }
	  else{
	    newt = pos->second;
	  }
	  //Now add the forward and reverse edges
	  //and flow properties
	  LGraph::edge_descriptor e1,e2;
	  tie(e1, inserted) = edge(news,newt,currlcbg);
	  if(!inserted){
	    tie(e1, inserted) = edge(newt,news,currlcbg);
	    assert(!inserted);
	    tie(e1, inserted) = add_edge(news,newt,currlcbg);
	    assert(inserted);
	    tie(e2, inserted) = add_edge(newt,news,currlcbg);
	    assert(inserted);
	    //put(edge_reverse,currlcbg,e1,e2);
	    //put(edge_reverse,currlcbg,e2,e1);
	    rev[e1] = e2;
	    assert(rev[e1]==e2);
	    rev[e2] = e1;
	    assert(rev[e2]==e1);
	    //Capacity is set as number of labels on the edge
	    BitMask emask = get(edge_labelmask,g,*out_i);
	    //int minlen = (get(vertex_len,g,news) < get(vertex_len,g,newt)) ? get(vertex_len,g,news) : get(vertex_len,g,newt);
	    //int ecapacity = emask.count() * minlen;
	    int ecapacity = emask.count();
	    assert(ecapacity>=1);
#ifdef DEBUG
	    std::cerr << "mincutlcbg " << sname << "-" << tname << " capacity:" << ecapacity << std::endl;
#endif
	    capacity[e1]=ecapacity;//DEFAULT_CAP;
	    capacity[e2]=ecapacity;//DEFAULT_CAP;
	    residual_capacity[e1]=0;
	    residual_capacity[e2]=0;
	  }
	}
      }
      

      //Condition (1) split gaps
      if(num_vertices(currlcbg)>0 
	 && num_edges(currlcbg)>0 
	 && sortedV.size()>0
	 && spanlen>=MINSPANLEN){ //Check that span of seqs > MINSPANLEN to avoid breaking LCBs based on inconsistent short fragments
	assert(num_vertices(currlcbg)>0);
	assert(num_edges(currlcbg)>0);
	assert(num_vertices(currlcbg)>=sortedV.size());
	assert(currcoords.size()==sortedV.size());

	//Project order onto seq
	sortedV.sort(coordsorder_vertex(&currcoords));
	
	int prevcoord=-1;
	int currpos=0;
	LVertex prevvertex=0,currvertex=0;
	VertexName prevname=0,currname=0;
	LVertex currvertexlcb;
	std::vector<int> lcbcc(num_vertices(currlcbg));

	prevname=0;
	prevcoord=-1;
#ifdef DEBUG
	std::cerr << "Order by seqidx:" << seqidx << std::endl;
#endif
	for(std::list<LVertex>::iterator vit = sortedV.begin();vit!=sortedV.end();++vit){
	  currvertex=*vit;
	  currname = get(vertex_name,currlcbg,currvertex);
	  assert(get(vertex_name,currlcbg,currvertex)==get(vertex_name,g,vmap[currvertex]));
	  currvertexlcb = currlcbv[currname];

	  assert(coordinates.find(std::make_pair(name2vertex[currname],seqidx))!=coordinates.end());
	  //assert(coordinates.find(std::make_pair(name2vertex[currname],seqidx))->second==currcoords.find(std::make_pair(currvertex,seqidx))->second);
	  int currstart,currend;
	  tie(currstart,currend) = coordinates.find(std::make_pair(name2vertex[currname],seqidx))->second;
	  if(prevcoord==-1){
	    assert(vit==sortedV.begin());
	  }
	  else{
	    //assert(*(vit-1)==prevvertex);
	    //assert(currstart>=prevcoord);
	    int dist = currstart-prevcoord;
	    
#ifdef DEBUG
	    std::cerr << "seqidx:"  << seqidx << " dist:" << dist << " " 
		      << prevname   << "-" << currname << " "   
		      << prevvertex   << "-" << currvertex  
		      << " coords " << prevcoord << "-" << currstart
		      << " spanlenonseq: " << spanlen  
		      << " numV: " << num_vertices(currlcbg) << std::endl;
#endif
	    if(dist>(int)distance){
	      //Since the vertices are sorted by genomic position.
	      //All verticies begin()->currVertex are also at a dist>distance
#ifdef DEBUG
	      std::cerr << "Found GAP " << dist << ">" << distance << std::endl;
#endif
#ifdef CALCFLOW
		;
#else
		boost::graph_traits<LGraph>::edge_descriptor de2;
		tie(de2,found) = edge(currvertex,prevvertex,currlcbg);
		if(!found){
		  tie(de2,found) = edge(prevvertex,currvertex,currlcbg);
		  if(found){
		    disconnecting_set.push_back(de2);
		  }
		}
		else{
		  disconnecting_set.push_back(de2);
		}

#endif	     
	      //Convert into multi-source multi-sink problem
	      //Add super-source and super sink nodes
	      LVertex ssource,ssink;
	      if(reusesupernodes && 
		 cuts.find(std::make_pair(prevname,currname)) != cuts.end()){
		ssource = cuts.find(std::make_pair(prevname,currname))->second.first;
		ssink = cuts.find(std::make_pair(prevname,currname))->second.second;
		//std::cerr << "Found prev source sink " << ssource << "-" << ssink << std::endl;
		assert(cutsdistmap.find(std::make_pair(ssource,ssink))!=cutsdistmap.end());
		if(dist<(int)cutsdistmap[std::make_pair(ssource,ssink)]){
		  cutsdistmap[std::make_pair(ssource,ssink)]=dist;
		}
	      }
	      else{
		ssource = add_vertex(std::numeric_limits<int>::max()-supercount,currlcbg);
		currlcbv[std::numeric_limits<int>::max()-supercount]=ssource;
		ssink = add_vertex(std::numeric_limits<int>::max()-supercount-1,currlcbg);
		currlcbv[std::numeric_limits<int>::max()-supercount-1]=ssink;
		supercount+=2;
		cuts[std::make_pair(prevname,currname)] = std::make_pair(ssource,ssink);
		revcuts[std::make_pair(ssource,ssink)] = std::make_pair(prevname,currname);
		cutsdistmap[std::make_pair(ssource,ssink)]=dist;
	      }
#ifdef DEBUG
		std::cerr << "Source,sink " << ssource << "-" << ssink << " for cut " << prevname << "-" << currname << endl;
#endif
	      
	      
	      std::list<LVertex>::iterator sinkend = sortedV.end();
	      std::list<LVertex>::iterator vit3 = vit;
	      //std::vector<LVertex>::iterator sinkend = (int(SEARCH_RADIUS+currpos)<(int)sortedV.size()) ? vit+SEARCH_RADIUS : sortedV.end();
	      for(;vit3!=sinkend;++vit3){
		//addFlowEdge(ssink,*vit3,currlcbg);
		graph_traits < LGraph >::edge_descriptor e1,e2;
		tie(e1, inserted) = add_edge(ssink,*vit3,currlcbg);
		if(inserted){
		  //std::cerr << "Adding edge " << get(vertex_name,currlcbg,*vit3) << " --> sink:" << ssink << std::endl;
		  snodeedges++;
		  tie(e2, inserted) = add_edge(*vit3,ssink,currlcbg);
		  assert(inserted);
		  snodeedges++;
		  //put(edge_reverse,currlcbg,e1,e2);
		  //put(edge_reverse,currlcbg,e2,e1);
		  rev[e1] = e2;
		  assert(rev[e1]==e2);
		  rev[e2] = e1;
		  assert(rev[e2]==e1);
		  capacity[e1]=0;
		  capacity[e2]=std::numeric_limits<int>::max();
		  residual_capacity[e1]=0;
		  residual_capacity[e2]=0;
		}
		else{
		  tie(e2, inserted) = add_edge(*vit3,ssink,currlcbg);
		  assert(!inserted);
		}
	      }

	      std::list<LVertex>::iterator sourceend = sortedV.begin();
	      //std::list<LVertex>::iterator sourceend = ((int)(currpos-SEARCH_RADIUS)>0) ? vit-SEARCH_RADIUS : sortedV.begin();
	      //std::cerr << "Curr pos " << currpos << std::endl;
	      std::list<LVertex>::iterator vit2=vit;
	      for(--vit2;vit2!=sourceend;--vit2){
		graph_traits < LGraph >::edge_descriptor e1,e2;
		tie(e1, inserted) = add_edge(ssource,*vit2,currlcbg);
		if(inserted){
		  //std::cerr << "Adding edge source:" << ssource << " --> " << get(vertex_name,currlcbg,*vit2) << std::endl;
		  snodeedges++;
		  tie(e2, inserted) = add_edge(*vit2,ssource,currlcbg);
		  snodeedges++;
		  assert(inserted);
		  //put(edge_reverse,currlcbg,e1,e2);
		  //put(edge_reverse,currlcbg,e2,e1);
		  rev[e1] = e2;
		  assert(rev[e1]==e2);
		  rev[e2] = e1;
		  assert(rev[e2]==e1);
		  capacity[e1]=std::numeric_limits<int>::max();
		  capacity[e2]=0;
		  residual_capacity[e1]=0;
		  residual_capacity[e2]=0;
		}
		else{
		  tie(e2, inserted) = add_edge(*vit2,ssource,currlcbg);
		  assert(!inserted);
		}
	      }
	      if(vit2==sortedV.begin()){
		graph_traits < LGraph >::edge_descriptor e1,e2;
		tie(e1, inserted) = add_edge(ssource,*vit2,currlcbg);
		if(inserted){
		  //std::cerr << "Adding edge source:" << ssource << " --> " << get(vertex_name,currlcbg,*vit2) << std::endl;
		  snodeedges++;
		  tie(e2, inserted) = add_edge(*vit2,ssource,currlcbg);
		  snodeedges++;
		  assert(inserted);
		  //put(edge_reverse,currlcbg,e1,e2);
		  //put(edge_reverse,currlcbg,e2,e1);
		  rev[e1] = e2;
		  assert(rev[e1]==e2);
		  rev[e2] = e1;
		  assert(rev[e2]==e1);
		  capacity[e1]=std::numeric_limits<int>::max();
		  capacity[e2]=0;
		  residual_capacity[e1]=0;
		  residual_capacity[e2]=0;
		}
		else{
		  tie(e2, inserted) = add_edge(*vit2,ssource,currlcbg);
		  assert(!inserted);
		}
	      }
	      if(cutsnodes.find(std::make_pair(ssource,ssink))==cutsnodes.end()){
		supernodes[ssource]++;
		supernodes[ssink]++;
		//supernodes2.insert(ssource);
		//supernodes2.insert(ssink);
		cutsnodes.insert(std::make_pair(ssource,ssink));
	      }
	    }
	  }
	  currpos++;
	  currvertex=*vit;
	  currname = get(vertex_name,currlcbg,currvertex);
	  assert(currvertex==(*vit));
	  prevvertex = currvertex;
	  //max coord of block
	  prevcoord = currend;
	  prevname = currname;
	}
      }
      else{
	//std::cerr << "Skipping merge on seq:" << seqidx 
	//<< " spanlen:" << spanlen << " < " << MINSPANLEN << std::endl;
      }
    }
#ifdef DEBUG
      std::cerr << "Graph built for lcbidx:" << lcbcount << " V:" << num_vertices(currlcbg) << " E:" << num_edges(currlcbg) << std::endl; 
#endif
    //
      //Condition (2) - conflicting orientation
      //Check orientation on this one seq      
      //Condition (3) - multiple seqs per genome
      //Need to break LCBs that have multiple seqs from the same genome
    LVertex prevvertex=0,currvertex=0;
    VertexName prevname=0,currname=0;
    LVertex currvertexlcb,prevvertexlcb;

    for(gpos = seqspergenomeMap.begin();gpos!=seqspergenomeMap.end();++gpos){
      if(gpos->second.size()>1){
#ifdef DEBUG
	std::cerr << "LCB with multiple seqs " << gpos->second.size() << " from same genome, splitting" << std::endl;
#endif
	std::vector<Label> seqs;
	for(spos = gpos->second.begin();spos!=gpos->second.end();++spos){//each seq1
	  seqs.push_back(*spos);
#ifdef DEBUG
	  std::cerr << "Seqs " << *spos << std::endl;
#endif
	  assert(sequence2genome[*spos]==gpos->first);
	}
	std::vector<LVertex> compv;
	//Split sequences from the same genome
	for(std::vector<Label>::iterator spos1 = seqs.begin();spos1!=seqs.end();++spos1){
	  //std::cerr << "S1" << *spos1 << std::endl;
	  assert(seqsvertex.find(*spos1)!=seqsvertex.end());
	  for(std::vector<Label>::iterator spos2 = spos1+1;spos2!=seqs.end();++spos2){
	    //std::cerr << "S2" << *spos2 << std::endl;
	    assert(seqsvertex.find(*spos2)!=seqsvertex.end());
	    assert(spos1!=spos2);
	    for(vpos = seqsvertex[*spos1].begin();vpos != seqsvertex[*spos1].end();++vpos){//each vertex seq1
	      compv.push_back(*vpos);
	      currvertex = *vpos;
	      currname = get(vertex_name,currlcbg,currvertex);
	      currvertexlcb = currlcbv[currname];
	      //std::cerr << *vpos << " name:" << currname << std::endl;
	      assert(get(vertex_name,currlcbg,currvertex)==get(vertex_name,g,vmap[currvertex]));
	      for(vpos2 = seqsvertex[*spos2].begin();vpos2 != seqsvertex[*spos2].end();++vpos2){//each vertex seq1
		prevvertex = *vpos2;
#ifdef CALCFLOW
		;
#else
		boost::graph_traits<LGraph>::edge_descriptor de2;
		tie(de2,found) = edge(currvertex,prevvertex,currlcbg);
		if(!found){
		  tie(de2,found) = edge(prevvertex,currvertex,currlcbg);
		  if(found){
		    disconnecting_set.push_back(de2);
		  }
		}
		else{
		  disconnecting_set.push_back(de2);
		}
#endif
		prevname = get(vertex_name,currlcbg,prevvertex);
		prevvertexlcb = get(vertex_name,currlcbg,prevvertex);
		//std::cerr << *vpos2 << " name:" << prevname << std::endl;
		LVertex ssource,ssink;
		ssource = add_vertex(std::numeric_limits<int>::max()-supercount,currlcbg);
		currlcbv[std::numeric_limits<int>::max()-supercount]=ssource;
		ssink = add_vertex(std::numeric_limits<int>::max()-supercount-1,currlcbg);
		currlcbv[std::numeric_limits<int>::max()-supercount-1]=ssink;
		supercount+=2;
		cuts[std::make_pair(prevname,currname)] = std::make_pair(ssource,ssink);
		revcuts[std::make_pair(ssource,ssink)] = std::make_pair(prevname,currname);
		cutsdistmap[std::make_pair(ssource,ssink)]=0;
			  
		graph_traits < LGraph >::edge_descriptor e1,e2;
		tie(e1, inserted) = add_edge(ssink,currvertex,currlcbg);

		if(inserted){
#ifdef DEBUG
		  std::cerr << "Added edge sink for multiple anchors same genome:" << currvertex << " <-- " << ssink << std::endl;
#endif
		  snodeedges++;
		  tie(e2, inserted) = add_edge(currvertex,ssink,currlcbg);
 		  assert(inserted);
		  snodeedges++;
		  rev[e1] = e2;
		  assert(rev[e1]==e2);
		  rev[e2] = e1;
		  assert(rev[e2]==e1);
		  capacity[e1]=0;
		  capacity[e2]=std::numeric_limits<int>::max();
		  residual_capacity[e1]=0;
		  residual_capacity[e2]=0;
		}
		else{
		  tie(e2, inserted) = add_edge(currvertex,ssink,currlcbg);
		  assert(!inserted);
		}
		tie(e1, inserted) = add_edge(ssource,prevvertex,currlcbg);
		if(inserted){
#ifdef DEBUG
		  std::cerr << "Adding edge source:" << ssource << " --> " << prevvertex << std::endl;
#endif
		  snodeedges++;
		  tie(e2, inserted) = add_edge(prevvertex,ssource,currlcbg);
		  snodeedges++;
		  assert(inserted);
		  rev[e1] = e2;
		  assert(rev[e1]==e2);
		  rev[e2] = e1;
		  assert(rev[e2]==e1);
		  capacity[e1]=std::numeric_limits<int>::max();
		  capacity[e2]=0;
		  residual_capacity[e1]=0;
		  residual_capacity[e2]=0;
		}
		else{
		  tie(e2, inserted) = add_edge(prevvertex,ssource,currlcbg);
		  assert(!inserted);
		}
		//Add ssource,ssink to cutset
		if(cutsnodes.find(std::make_pair(ssource,ssink))==cutsnodes.end()){
		  supernodes[ssource]++;
		  supernodes[ssink]++;
		  cutsnodes.insert(std::make_pair(ssource,ssink));
		}
	      }
	    }
	  }
	}
	/*
	//Check for misoriented vertices within an LCB and break
	for(std::vector<LVertex>::iterator vpos = compv.begin();vpos != compv.end();++vpos){//each vertex seq1
	  for(std::vector<LVertex>::iterator vpos2 = vpos+1;vpos2 != compv.end();++vpos2){//each vertex seq1
	    //Mismatched orient
	    BitMask sharedlabels = (labelmaskmap[vmap[*vpos]]&labelmaskmap[vmap[*vpos2]]);
	    assert(isLabelCollinearMask(sharedlabels,
					orientmaskmap[vmap[*vpos]],
					orientmaskmap[vmap[*vpos2]]) 
		   == 
		   isLabelCollinear(orientmap[vmap[*vpos]],
				    orientmap[vmap[*vpos2]],
				    sequence2genome));
	    if(! isLabelCollinearMask(sharedlabels,orientmaskmap[vmap[*vpos]],orientmaskmap[vmap[*vpos2]])){
	      std::cerr << "Breaking vertices with incompatible labeling " << vmap[*vpos] << "--" << vmap[*vpos2] << std::endl;
	      currvertex = *vpos;
	      currname = get(vertex_name,currlcbg,currvertex);
	      currvertexlcb = currlcbv[currname];
	      std::cerr << *vpos << " name:" << currname << std::endl;
	      assert(get(vertex_name,currlcbg,currvertex)==get(vertex_name,g,vmap[currvertex]));
	      prevvertex = *vpos2;
	      prevname = get(vertex_name,currlcbg,prevvertex);
	      prevvertexlcb = get(vertex_name,currlcbg,prevvertex);
	      std::cerr << *vpos2 << " name:" << prevname << std::endl;
	      LVertex ssource,ssink;
	      ssource = add_vertex(std::numeric_limits<int>::max()-supercount,currlcbg);
	      currlcbv[std::numeric_limits<int>::max()-supercount]=ssource;
	      ssink = add_vertex(std::numeric_limits<int>::max()-supercount-1,currlcbg);
	      currlcbv[std::numeric_limits<int>::max()-supercount-1]=ssink;
	      supercount+=2;
	      cuts[std::make_pair(prevname,currname)] = std::make_pair(ssource,ssink);
	      revcuts[std::make_pair(ssource,ssink)] = std::make_pair(prevname,currname);
	      cutsdistmap[std::make_pair(ssource,ssink)]=0;
	      
	      graph_traits < LGraph >::edge_descriptor e1,e2;
	      tie(e1, inserted) = add_edge(ssink,currvertex,currlcbg);
	      if(inserted){
		std::cerr << "Added edge sink:" << currvertex << " <-- " << ssink << std::endl;
		snodeedges++;
		tie(e2, inserted) = add_edge(currvertex,ssink,currlcbg);
		std::cerr << "Added edge sink:" << currvertex << " <-- " << ssink << std::endl;
		assert(inserted);
		snodeedges++;
		rev[e1] = e2;
		assert(rev[e1]==e2);
		rev[e2] = e1;
		assert(rev[e2]==e1);
		capacity[e1]=0;
		capacity[e2]=std::numeric_limits<int>::max();
		residual_capacity[e1]=0;
		residual_capacity[e2]=0;
	      }
	      else{
		tie(e2, inserted) = add_edge(currvertex,ssink,currlcbg);
		assert(!inserted);
	      }
	      std::cerr << "Added sink" << std::endl;
	      tie(e1, inserted) = add_edge(ssource,prevvertex,currlcbg);
	      if(inserted){
		std::cerr << "Adding edge source:" << ssource << " --> " << prevvertex << std::endl;
		snodeedges++;
		tie(e2, inserted) = add_edge(prevvertex,ssource,currlcbg);
		snodeedges++;
		assert(inserted);
		rev[e1] = e2;
		assert(rev[e1]==e2);
		rev[e2] = e1;
		assert(rev[e2]==e1);
		capacity[e1]=std::numeric_limits<int>::max();
		capacity[e2]=0;
		residual_capacity[e1]=0;
		residual_capacity[e2]=0;
	      }
	      else{
		tie(e2, inserted) = add_edge(prevvertex,ssource,currlcbg);
		assert(!inserted);
	      }
	      //Add ssource,ssink to cutset
	      if(cutsnodes.find(std::make_pair(ssource,ssink))==cutsnodes.end()){
		supernodes[ssource]++;
		supernodes[ssink]++;
		  cutsnodes.insert(std::make_pair(ssource,ssink));
	      }
	    }
	    else{
	      std::cerr << "Compatible labeling " << vmap[*vpos] << "--" << vmap[*vpos2] << std::endl;
	    }
	  }
	
	  }
	*/
      }
    }
  
    //
    //
    
#ifdef PRINTFLOW
    //Write graph
    std::vector<int> ccvmap; //empty
    VertexSet maskedLCBs; //empty
#ifdef PRINTSEQS
    ;
#else
    do_write_graphviz(currlcbg, std::string("gout.preflow"+lexical_cast<std::string>(cutcount+filenumoffset)+".dot"),ccvmap,coordinates,maskedEdges,maskedLCBs,capacity,false);
    std::cerr << "Writing " << std::string("gout.preflow"+lexical_cast<std::string>(cutcount+filenumoffset)+".dot") << std::endl;
#endif
#endif
    
    LGraph::edge_iterator ei,e_end;

    //property_map < LGraph, edge_reverse_t >::type revtest = get(edge_reverse,currlcbg);
    //for(tie(ei, e_end) = edges(currlcbg); ei != e_end; ++ei) {
    //assert(revtest[revtest[*ei]] == *ei); //check if the reverse edge map is build up properly
    //}
    //Evaluation order of cuts can matter
    //TODO, try smallest->largest and largest->smallest
    std::vector<std::pair<LVertex,LVertex> > cutsnodesuniq;
    for(std::set<pair<LVertex,LVertex> >::iterator cit = cutsnodes.begin(); cit!= cutsnodes.end();++cit){
      cutsnodesuniq.push_back(*cit);
    }
    sort(cutsnodesuniq.begin(),cutsnodesuniq.end(),cutsdist(&cutsdistmap));
    for(std::vector<std::pair<LVertex,LVertex> >::iterator cit = cutsnodesuniq.begin(); cit!= cutsnodesuniq.end();++cit){
      LVertex ssource = cit->first;
      LVertex ssink = cit->second;
#ifdef DEBUG
	std::cerr << "Attempting split " << get(vertex_name,currlcbg,ssource) << "(" << supernodes[ssource] << ")" 
		  << "-" << get(vertex_name,currlcbg,ssink)  << "(" << supernodes[ssink] << ")" 
		  << " due to edge " <<  revcuts[std::make_pair(ssource,ssink)].first
		  << "-" <<  revcuts[std::make_pair(ssource,ssink)].second 
		  << " dist:" << cutsdistmap[std::make_pair(ssource,ssink)] << std::endl;
#endif
      
      assert(supernodes[ssource]>0);
      assert(supernodes[ssink]>0);

            
      std::vector<default_color_type> color(num_vertices(currlcbg));
      std::vector<LGraph::edge_descriptor> pred(num_vertices(currlcbg));
      
      assert(num_edges(currlcbg)>0);

      std::set<LVertex> S_star;

      
      property_map < LGraph, vertex_index_t >::type
	idx = get(vertex_index, currlcbg);
      property_map < LGraph, vertex_distance_t >::type
	distance = get(vertex_distance, currlcbg); 
      capacity = get(edge_capacity, currlcbg);
      rev = get(edge_reverse, currlcbg);
      residual_capacity = get(edge_residual_capacity, currlcbg);
#ifdef CALCFLOW
      long flow = edmonds_karp_max_flow(currlcbg, ssource, ssink, capacity, residual_capacity, rev, &color[0], &pred[0]);
      ++cutcount;
      //kolmogorov is faster but its not clear how to find the disconnecting set 
      //long flow = kolmogorov_max_flow(currlcbg, capacity, residual_capacity, rev, &pred[0], &color[0],distance,idx,ssource,ssink);      

      //long flow = edmonds_karp_max_flow(currlcbg, ssource, ssink);
      //long flow = push_relabel_max_flow(currlcbg, ssource, ssink);
      //long flow = kolmogorov_max_flow(currlcbg, ssource, ssink);
      
      /*
      //Testing trimming graph of all but the current supernode source and sink
      typedef std::set<LVertex> SuperNodeMap;
      typedef filtered_graph<LGraph, 
	snode_efilter<LGraph>,snode_vfilter<LGraph> > FLGraph;

      snode_efilter<LGraph> efilter(&supernodes2,&currlcbg);
      snode_vfilter<LGraph> vfilter(&supernodes2);
      FLGraph filtlcbg(currlcbg, efilter, vfilter);
      supernodes2.erase(supernodes2.find(ssource));
      supernodes2.erase(supernodes2.find(ssink));
      long flow = kolmogorov_max_flow(filtlcbg, capacity, residual_capacity, rev, &pred[0], &color[0],distance,idx,ssource,ssink);
      */

      //Check flow since we may have already introduced a break 
      if(flow>0){
	assert(flow>0);
	graph_traits<LGraph>::out_edge_iterator ei, ei_end;
	graph_traits<LGraph>::vertex_iterator vi, vi_end;
	typedef color_traits<default_color_type> Color;
	for(tie(vi,vi_end) = vertices(currlcbg);vi!=vi_end;++vi){
	  if(color[*vi]!=Color::white()){
	    if(reusesupernodes || supernodes.find(*vi)==supernodes.end()){
	      S_star.insert(*vi);
	    }
	  }
	}
	for( std::set<LVertex>::iterator si = S_star.begin();si!=S_star.end();++si){
	  for(tie(ei,ei_end) = out_edges(*si,currlcbg);ei!=ei_end;++ei){
	    if(S_star.find(target(*ei,currlcbg))==S_star.end()){ 
	      if(reusesupernodes || supernodes.find(target(*ei,currlcbg))==supernodes.end()){
		disconnecting_set.push_back(*ei);
		#ifdef DEBUG
		std::cerr << "Disconnecting set " << get(vertex_name,currlcbg,source(*ei,currlcbg)) << "-" << get(vertex_name,currlcbg,target(*ei,currlcbg)) << std::endl; 
		put(edge_category,currlcbg,*ei,BLUE);
		#endif
	      }
	    }
	  }
	}
	#ifdef DEBUG
	std::cerr << " flow:" << flow << std::endl;
	#endif

#else
	int flow=0;
#endif
#ifdef DEBUG
	std::cerr << "Disconnecting set size:" << disconnecting_set.size() << std::endl; 
#endif

#ifdef PRINTFLOW
	//Write graph
	std::vector<int> ccvmap; //empty
	VertexSet maskedLCBs; //empty
#ifdef PRINTSEQS
	;
#else
	do_write_graphviz(currlcbg, std::string("gout.flow"+lexical_cast<std::string>(cutcount+filenumoffset)+".dot"),ccvmap,coordinates,maskedEdges,maskedLCBs,capacity,false);
	std::cerr << "Writing " << std::string("gout.flow"+lexical_cast<std::string>(cutcount+filenumoffset)+".dot") << std::endl;
#endif
#endif
	for(std::vector<LGraph::edge_descriptor>::iterator ei=disconnecting_set.begin();ei!=disconnecting_set.end();++ei){
	  boost::graph_traits<LGraph>::edge_descriptor e2;
	  //This edge may have been cut previously so
	  //first check if it is still present in the connectivity graph
	  tie(e2,found) = edge(source(*ei,currlcbg),target(*ei,currlcbg),currlcbg);
	  if(found){
	    LGraph::edge_descriptor maxe = *ei;
	    Vertex cuts = name2vertex[get(vertex_name,currlcbg,source(maxe,currlcbg))];
	    assert(vmap[source(maxe,currlcbg)]==cuts);
	    Vertex cutt = name2vertex[get(vertex_name,currlcbg,target(maxe,currlcbg))];
	    assert(vmap[target(maxe,currlcbg)]==cutt);

	    Edge e1;
	    tie(e1,found) = edge(cuts,cutt,fglcbsyn);
	    if(found){
	      assert(revcuts.find(std::make_pair(ssource,ssink))!=revcuts.end());
#ifdef DEBUG
		std::cerr << "cut " << get(vertex_name,currlcbg,source(maxe,currlcbg)) 
			  << "-" << get(vertex_name,currlcbg,target(maxe,currlcbg)) 
			  << " " << cuts << "-" << cutt
			  << std::endl;
		if((get(vertex_name,currlcbg,target(maxe,currlcbg)) == revcuts[std::make_pair(ssource,ssink)].first 
		    &&  get(vertex_name,currlcbg,source(maxe,currlcbg)) == revcuts[std::make_pair(ssource,ssink)].second)
		   ||
		   (get(vertex_name,currlcbg,source(maxe,currlcbg)) == revcuts[std::make_pair(ssource,ssink)].first 
		    &&  get(vertex_name,currlcbg,target(maxe,currlcbg)) == revcuts[std::make_pair(ssource,ssink)].second)
		   ){
		  std::cerr << "Split is trivial" << std::endl;
		}
		else{
		  std::cerr << "Split is non-local " << std::string("gout.flow"+lexical_cast<std::string>(cutcount+filenumoffset)) << " " << " supernodes:" << supernodes.size() << " " << "cut_set:" << disconnecting_set.size() << std::endl;
		}
		
#endif
	      maskedEdges.insert(std::make_pair(cuts,cutt)); 
	      put(edge_category,g,e1,BLUE);
	    }
	    else{
	      assert(revcuts.find(std::make_pair(ssource,ssink))!=revcuts.end());
	      tie(e1,found) = edge(cutt,cuts,fglcbsyn);
	      if(found){
#ifdef DEBUG
		std::cerr << "cut " << get(vertex_name,currlcbg,target(maxe,currlcbg)) 
			  << "-" << get(vertex_name,currlcbg,source(maxe,currlcbg))
			  << " " << cuts << "-" << cutt
			  << std::endl;
		if((get(vertex_name,currlcbg,target(maxe,currlcbg)) == revcuts[std::make_pair(ssource,ssink)].first 
		    &&  get(vertex_name,currlcbg,source(maxe,currlcbg)) == revcuts[std::make_pair(ssource,ssink)].second)
		   ||
		   (get(vertex_name,currlcbg,source(maxe,currlcbg)) == revcuts[std::make_pair(ssource,ssink)].first 
		    &&  get(vertex_name,currlcbg,target(maxe,currlcbg)) == revcuts[std::make_pair(ssource,ssink)].second)
		   ){
		  std::cerr << "Split is trivial" << std::endl;
		}
		else{
		  std::cerr << "Split is non-local " << std::string("gout.flow"+lexical_cast<std::string>(cutcount+filenumoffset)) << " " << " supernodes:" << supernodes.size() << " " << "cut_set:" << disconnecting_set.size() << std::endl;
		}
#endif
		maskedEdges.insert(std::make_pair(cutt,cuts));
		put(edge_category,g,e1,BLUE);
	      }
	    }
	    tie(e1,found) = edge(cutt,cuts,fglcbsyn);
	    assert(!found);
	    tie(e1,found) = edge(cuts,cutt,fglcbsyn);
	    assert(!found);
	    //Removing the edges from the graph can help short circuit future runs of the maxflow/mincut algorithm
	    //The check above for flow>0 ensures subsequent cuts are only considered if there is still 
	    //connectivity in the graph
	    //TODO, consider changing capacity to zero instead of removing for perf boost
	    //#ifdef DEBUG
	    //save edges so we can visualize
	    //#else
	    remove_edge(rev[maxe],currlcbg);
	    remove_edge(maxe,currlcbg);
	    //#endif
	  }
	  else{
	    //assert(false);
	  }
	}
#ifdef DEBUG
	//TESTING
	//Check to make sure the cuts eliminated all the flow
	//This is for testing only
	flow = kolmogorov_max_flow(currlcbg, capacity, residual_capacity, rev, &pred[0], &color[0],distance,idx,ssource,ssink);      
	std::cerr << "Remaining flow " << flow << std::endl;
	assert(flow==0);
#endif
	numcuts+=disconnecting_set.size();
#ifdef CALCFLOW
    }
    else{
      //Previous cut already broke flow between ssource-ssink. No further cuts needed
    }
#endif
      //supernodes2.insert(ssource);
      //supernodes2.insert(ssink);

      //Remove supernodes only if 
      //they are no longer referenced
      /*
	This doesn't work as expected. Can't get the vertex and
	associated edges to properly clear to property maps
	supernodes[ssource]--;
      supernodes[ssink]--;

      if(supernodes[ssource]==0){
	supernodes.erase(supernodes.find(ssource));
	clear_vertex(ssource,currlcbg);
	remove_vertex(ssource,currlcbg);
      }
      if(supernodes[ssink]==0){
	supernodes.erase(supernodes.find(ssink));
	clear_vertex(ssink,currlcbg);
	remove_vertex(ssink,currlcbg);
      }
      LGraph::edge_iterator ei,e_end;
      property_map < LGraph, edge_reverse_t >::type rev2 = get(edge_reverse, currlcbg);
      for(tie(ei, e_end) = edges(currlcbg); ei != e_end; ++ei) {
	std::cerr << *ei << std::endl;
	//This will segfault after removed nodes
	assert(rev2[rev2[*ei]]==*ei);
      }
      */
    }
  }
  return numcuts;
}



/*

template<typename TGraph, typename TFGraph, typename LCBGraph, typename VertexMap>
void createLCBGraph(TGraph & g, TFGraph & fglcbsyn, LCBGraph & currlcbg, LCB & lcb, VertexMap & vmap){
  int DEFAULT_CAP=1;
  property_map < LGraph, edge_capacity_t >::type
      capacity = get(edge_capacity, currlcbg);
  property_map < LGraph, edge_reverse_t >::type 
    rev = get(edge_reverse, currlcbg);
  property_map < LGraph, edge_residual_capacity_t >::type
    residual_capacity = get(edge_residual_capacity, currlcbg);
  
  std::map<VertexName, LVertex> currlcbv;
  std::map<VertexName, LVertex>::iterator pos;
  bool inserted;
  for(LCB::iterator vit = lcb.begin();vit!=lcb.end();++vit){
    Vertex v=*vit;
    VertexName sname = get(vertex_name,g,v);
    LVertex news;
    //
    //Insert vertex into currlcbg if needed
    tie(pos, inserted) = currlcbv.insert(std::make_pair(sname, LVertex()));
    if(inserted){
      news = add_vertex(sname,currlcbg);
      currlcbv[sname]=news;
      vmap[news]=v;
    }
    else{
      news = pos->second;
    }
    //
    //Add all edges for news
    //First make sure target vertex is part of currlcbg
    graph_traits<LCBSynFilterGraph>::out_edge_iterator out_i, out_end;
    for(tie(out_i, out_end) = out_edges(v, fglcbsyn); out_i != out_end; ++out_i){
      VertexName tname = get(vertex_name,g,target(*out_i,g));
      LVertex newt;
      tie(pos, inserted) = currlcbv.insert(std::make_pair(tname, LVertex()));
      if(inserted){
	newt = add_vertex(tname,currlcbg);
	currlcbv[tname]=newt;
	vmap[newt]=target(*out_i,g);
      }
      else{
	newt = pos->second;
      }
      //Now add the forward and reverse edges
      //and flow properties
      LGraph::edge_descriptor e1,e2;
      tie(e1, inserted) = edge(news,newt,currlcbg);
      if(!inserted){
	tie(e1, inserted) = edge(newt,news,currlcbg);
	assert(!inserted);
	tie(e1, inserted) = add_edge(news,newt,currlcbg);
	assert(inserted);
	tie(e2, inserted) = add_edge(newt,news,currlcbg);
	assert(inserted);
	//put(edge_reverse,currlcbg,e1,e2);
	//put(edge_reverse,currlcbg,e2,e1);
	rev[e1] = e2;
	assert(rev[e1]==e2);
	rev[e2] = e1;
	assert(rev[e2]==e1);
	capacity[e1]=DEFAULT_CAP;
	capacity[e2]=DEFAULT_CAP;
	residual_capacity[e1]=0;
	residual_capacity[e2]=0;
      }
    }
  }
}
template<typename TGraph,
	 typename TEdge>
void addFlowEdge(TGraph & g, 
		 TEdge & ss, 
		 TEdge & e,
		 rev,
		 residual_capacity){
  graph_traits < LGraph >::edge_descriptor e1,e2;
		tie(e1, inserted) = add_edge(ssink,*vit3,currlcbg);
		if(inserted){
		  //std::cerr << "Adding edge " << get(vertex_name,currlcbg,*vit3) << " --> sink:" << ssink << std::endl;
		  snodeedges++;
		  tie(e2, inserted) = add_edge(*vit3,ssink,currlcbg);
		  assert(inserted);
		  snodeedges++;
		  //put(edge_reverse,currlcbg,e1,e2);
		  //put(edge_reverse,currlcbg,e2,e1);
		  rev[e1] = e2;
		  assert(rev[e1]==e2);
		  rev[e2] = e1;
		  assert(rev[e2]==e1);
		  capacity[e1]=0;
		  capacity[e2]=std::numeric_limits<int>::max();
		  residual_capacity[e1]=0;
		  residual_capacity[e2]=0;
}

*/

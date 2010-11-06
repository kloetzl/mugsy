
//#########################
//Predicates
//
//isLabelMaxGap()
//Test if v1 and v2 are separated by < maxgap in label intersection(s1,s2)
bool isLabelMaxGap(Vertex v1, 
		   Vertex v2,
		   OrientedLabelSet & s1, 
		   OrientedLabelSet & s2,
		   VertexLabelIntervalMap & coordinates,
		   unsigned int maxgap,
		   SequenceGenomeMap & sequence2genome){

  OrientedLabelSet::iterator s1_it_end = s1.end();
  for(OrientedLabelSet::iterator s1_it=s1.begin();s1_it!=s1_it_end;++s1_it){
    Label seqidx = s1_it->first;
    //Label genomeidx = sequence2genome[seqidx];
    std::list<Vertex> sortedV;
    VertexIntervalMap currcoords;
    //only consider seqs present in both s1 and s2
    OrientedLabelSet::iterator s2_it = s2.find(*s1_it);
    if(s2_it != s2.end()){
      assert(seqidx==s2_it->first);
      sortedV.push_back(v1);
      sortedV.push_back(v2);
      currcoords.insert(std::make_pair(v1,
				       coordinates[std::make_pair(v1,seqidx)]));
      currcoords.insert(std::make_pair(v2,
				       coordinates[std::make_pair(v2,seqidx)]));

      //sort(sortedV.begin(),sortedV.end(),coordsorder(&coordinates,seqidx));
      //sortedV.sort(coordsorder(&coordinates,seqidx));
      sortedV.sort(coordsorder_vertex(&currcoords));

      int prevcoord=-1;
      int currstart,currend;
      std::list<Vertex>::iterator vit_end = sortedV.end();
      for(std::list<Vertex>::iterator vit = sortedV.begin();vit!=vit_end;++vit){
	assert(coordinates.find(std::make_pair(*vit,seqidx)) != coordinates.end());
	assert(coordinates.find(std::make_pair(*vit,seqidx))->second == currcoords.find(*vit)->second);
	tie(currstart,currend) = currcoords[*vit];
	if(prevcoord==-1){
	  assert(vit==sortedV.begin());
	}
	else{
	  //assert(coordinates.find(std::make_pair(*(vit-1),seqidx)) != coordinates.end());
	  //assert(coordinates.find(std::make_pair(*(vit-1),seqidx))->second.second == prevcoord);
	  assert(currstart<currend);
	  //assert(currstart>=prevcoord);
	  int dist = currstart-prevcoord;
	  if(dist>(int)maxgap){
	    return false;
	  }
	}
	prevcoord=currend;
      }
    }
  }
  return true;
}

//isLabelCollinear(s1,s2)
//
//Return true if there is no change in orientation between labels s1
//and s2.  This implies a collinear relationship between s1 and s2
//meaning there are no rearrangments between the labels in s1 and
//s2. Labels are comprised of a pair (seq,orient).  There are 2
//possibilities for a collinear relationship,

//S=intersection(seq(s1),seq(s2))
//return orient(s1 in S)==orient(s2 in S)
//OR
//revorient(s1 in S)==orient(s2 in S)

//In other words, this function checks both the stored orientation of
//s1 vs. s2 and rev(s1) vs s2 for all sequences shared between s1 and
//s2 and returns true if either comparison is collinear (ie, no change in
//orientation)
inline bool isLabelCollinear(OrientedLabelSet & s1, 
			     OrientedLabelSet & s2, 
			     SequenceGenomeMap & sequence2genome){
  //Implemented using bitmasks to store the presence of a sequence
  //and the orientation
  BitMask s1lmask,s2lmask;
  BitMask s1omask,s2omask,s1omaskrev;
  OrientedLabelSet::iterator s1_it_end = s1.end();
  for(OrientedLabelSet::iterator s1_it=s1.begin();s1_it!=s1_it_end;++s1_it){
    Label seqidx = s1_it->first;
    Label genomeidx = sequence2genome[seqidx];
    assert(genomeidx>=0&&genomeidx<BITMAX);
    s1lmask.set(genomeidx,1);
    //set omask for all + orientation
    if(s1_it->second == true){
      s1omask.set(genomeidx,1);
    }
  }
  s1omaskrev=s1omask;
  s1omaskrev.flip();
  OrientedLabelSet::iterator s2_it_end = s2.end();
  for(OrientedLabelSet::iterator s2_it=s2.begin();s2_it!=s2_it_end;++s2_it){
    Label seqidx = s2_it->first;
    Label genomeidx = sequence2genome[seqidx];
    assert(genomeidx>=0&&genomeidx<BITMAX);
    s2lmask.set(genomeidx,1);
    //set omask for all + orientation
    if(s2_it->second == true){
      s2omask.set(genomeidx,1);
    }
  }

  //Shared labels obtain by intersection using bitwise AND
  BitMask sharedlabels = (s1lmask&s2lmask);  

#if defined(V_DEBUG)
  BitMask difflabels = (s1lmask^s2lmask)&(s1lmask|s2lmask);  
  cout << "S1MASK: " << s1lmask << endl;
  cout << "S2MASK: " << s2lmask << endl;
  cout << "SHARED: " << sharedlabels << endl;
  cout << "DIFF: " << difflabels << endl;
  cout << "O1MASK: " << s1omask << endl;
  cout << "O2MASK: " << s2omask << endl;
  cout << "O1MASKa: " << (s1omask&sharedlabels) << endl;
  cout << "O2MASKa: " << (s2omask&sharedlabels) << endl;
#endif

  if((s1omask&sharedlabels)==(s2omask&sharedlabels) || 
     (s1omaskrev&sharedlabels)==(s2omask&sharedlabels)){
    return true;
  }
  else{
    return false;
  }
}

inline bool isLabelCollinearMask(BitMask & sharedlabels, BitMask & s1omask, BitMask & s2omask){
  if((s1omask&sharedlabels)==(s2omask&sharedlabels)){
    return true; 
  }
  else{
    BitMask s1omaskrev = s1omask;
    s1omaskrev.flip();
    if((s1omaskrev&sharedlabels)==(s2omask&sharedlabels)){
      return true;
    }
    else{
#ifdef DEBUG
      //std::cerr << (s1omask&sharedlabels) << std::endl
      //	  << (s2omask&sharedlabels) << std::endl << std::endl
      //	  << (s1omaskrev&sharedlabels) << std::endl
      //	  << (s2omask&sharedlabels) << std::endl;
#endif
      return false;
    }
  }
  assert(false);
}

//
//Return true if no LCB gaps > maxgap
template<typename TGraph> inline
bool checkLCBGaps(TGraph & g,
		  LCB & lcb,
		  std::vector<int> & ccvmap,
		  VertexLabelIntervalMap & coordinates,
		  unsigned int maxgap,
		  SequenceGenomeMap & sequence2genome){
  bool shortcircuit=true;

  bool badGap=false;
  int MINSPANLEN=0;
  LabelSet seqidxSet;
  std::set<Vertex> mmV;
  std::map<Label,std::set<Label> > seqspergenomeMap; //tracks the number of seqs per genome in an LCB
  std::map<Label,std::set<Label> >::iterator gpos;
  bool inserted;

  typename property_map <TGraph, vertex_orient_t >::type orientmap = get(vertex_orient,g);

  LCB::iterator lit_end = lcb.end();
  for(LCB::iterator lit=lcb.begin();lit!=lit_end;++lit){ 
    Vertex v = *lit;
    //OrientedLabelSet o1 = get(vertex_orient, g, v);
    OrientedLabelSet::iterator oit_end = orientmap[v].end();
    for(OrientedLabelSet::iterator oit=orientmap[v].begin();oit!=oit_end;++oit){ //all seqs on the vertex
      Label seqidx = oit->first;
      seqidxSet.insert(seqidx);
      tie(gpos, inserted) = seqspergenomeMap.insert(std::make_pair(sequence2genome[seqidx],std::set<Label>()));
      gpos->second.insert(seqidx);
      if(gpos->second.size()>1){
	return false;
      }
    }
  }

  LabelSet::iterator it2_end = seqidxSet.end();
  for(LabelSet::iterator it2 = seqidxSet.begin(); it2 != it2_end; ++it2){
    Label seqidx = *it2;
    //Label genomeidx = sequence2genome[seqidx];
    std::list<Vertex> sortedV;
    unsigned int spanlen=0;
    //std::map<pair<Vertex,Label>,Interval > currcoords;
    VertexIntervalMap currcoords;

    LCB::iterator lit_end = lcb.end();
    for(LCB::iterator lit=lcb.begin();lit!=lit_end;++lit){ 
      Vertex v = *lit;
      VertexLabelIntervalMap::iterator cit = coordinates.find(std::make_pair(v,seqidx));
      if(cit!=coordinates.end()){
	assert(cit->second.first<cit->second.second);
	sortedV.push_back(v);
	currcoords.insert(std::make_pair(v,
					 cit->second));
	spanlen = spanlen + get(vertex_len,g,v);
      }
    }
#ifdef DEBUG
    std::cerr << "checkLCBGaps seqidx: " << seqidx << " span length " << spanlen  << " MINSPANLEN " << MINSPANLEN << std::endl;
#endif
    if(spanlen >= MINSPANLEN){
      sortedV.sort(coordsorder_vertex(&currcoords));
      
      int prevcoord=-1;
      int currstart,currend;
      Vertex prevvertex=0;
      std::list<Vertex>::iterator vit_end = sortedV.end();
      for(std::list<Vertex>::iterator vit = sortedV.begin();vit!=vit_end;++vit){
	assert(coordinates.find(std::make_pair(*vit,seqidx)) != coordinates.end());
	assert(coordinates[std::make_pair(*vit,seqidx)] == currcoords[*vit]);
	tie(currstart,currend) = currcoords[*vit];
	if(prevcoord==-1){
	  assert(vit==sortedV.begin());
	  prevvertex=*vit;
	}
	else{
	  assert(coordinates.find(std::make_pair(prevvertex,seqidx)) != coordinates.end());
	  assert(coordinates.find(std::make_pair(prevvertex,seqidx))->second.second == prevcoord);
	  assert(currstart<currend);
	  int dist = currstart-prevcoord;
#ifdef DEBUG
	    std::cerr << "Checking dist:" << dist << " > " << maxgap
		      << " between V:" <<  get(vertex_name,g,prevvertex) << " " << currstart << "-" << currend 
		      << " and V:" <<  get(vertex_name,g,*vit)
		      << " on seqidx:" << seqidx << std::endl;
#endif

	  if(dist>(int)maxgap){
	    badGap=true;
#ifdef DEBUG
	    std::cerr << "Large gap dist:" << dist << " > " << maxgap
		      << " between V:" <<  get(vertex_name,g,prevvertex) 
		      << " and V:" <<  get(vertex_name,g,*vit)
		      << " on seqidx:" << seqidx << std::endl;
#endif
	    if(shortcircuit){
	      return !badGap;
	    }
	    else{
	      mmV.insert(*vit);
	    }
	    //std::cerr << "NO SHORT CIRCUIT" << std::endl;
	  }
	}
	prevvertex=*vit;
	prevcoord=currend;
      }
    }
    else{
      //std::cerr << "Skipping check of seqidx: " << seqidx << " span length " << spanlen  << " MINSPANLEN " << MINSPANLEN << std::endl;
    }
#ifdef DEBUG
    if(mmV.size()>0){
      std::cerr << "Num bad vertices:" << mmV.size()  
		<< " LCB size:" << lcb.size() << std::endl;
    }
#endif
  }
  return !badGap;
}

//
//Bitmask implementation
inline bool checkLCBOrient(BitMask & lcbl1,
			   BitMask & lcbl2,
			   BitMask & lcbo1,
			   BitMask & lcbo2){
#ifdef DEBUG
    std::cerr << "LCBOrient l1:" << lcbl1 << " l2:" << lcbl2 << std::endl
	      << "LCBOrient o1:" << lcbo1 << " o2:" << lcbo2 << std::endl;
#endif
  BitMask sharedlabels = (lcbl1&lcbl2);
#ifdef DEBUG
    std::cerr << "LCBOrient shared:" << sharedlabels << std::endl;
#endif
  //Make sure this function is working symmetrically
  assert(isLabelCollinearMask(sharedlabels,lcbo1,lcbo2)==isLabelCollinearMask(sharedlabels,lcbo2,lcbo1));
  if(isLabelCollinearMask(sharedlabels,lcbo1,lcbo2)){
#ifdef DEBUG
      std::cerr << "Match" << std::endl;
#endif
    return true;
  }
  else{
#ifdef DEBUG
      std::cerr << "MisMatch" << std::endl;
#endif
    return false;
  }
}

template<typename TLCBOrientMap>
bool checkLCBOrient(TLCBOrientMap & lcborientmap,
		    typename TLCBOrientMap::key_type &lcbidx1, 
		    typename TLCBOrientMap::key_type &lcbidx2,
		    BitMask &longlabelmask){
  BitMask t1=lcborientmap[lcbidx1].first&longlabelmask;
  BitMask t2=lcborientmap[lcbidx2].first&longlabelmask;
  BitMask t3=lcborientmap[lcbidx1].second&longlabelmask;
  BitMask t4=lcborientmap[lcbidx2].second&longlabelmask;
  return checkLCBOrient(t1,t2,t3,t4);
}

//
//Return true if no mismatch between Vertex orientations within an LCB
//Currently implemented with bitmasks
//TODO
//Improve perfomance
//Profiling shows this check is primary performance bottleneck
//First attempt above was to cache orientation of the LCB once instead of checking all pairs
//for consisency
//Surprisingly, this does not appear to improve performance?  
template<typename TGraph> inline
bool checkLCBOrient(TGraph & g,
		    LCB & lcb,
		    BitMask &lcbmask1,
		    BitMask &lcbmask2,
		    BitMask &longlabelmask){
  

  LCB::iterator it,it2,it_end,it2_end;

  //Use orientmap so that we can pass by reference using lvalue []
  typename property_map < TGraph, vertex_orient_t >::type orientmap = get(vertex_orient,g);
  typename property_map < TGraph, vertex_orientmask_t >::type orientmaskmap = get(vertex_orientmask,g);
  typename property_map < TGraph, vertex_vlabelmask_t >::type labelmaskmap = get(vertex_vlabelmask,g);
  it_end = lcb.end();
  it2_end = lcb.end();
  BitMask l1 = lcbmask1&longlabelmask;
  BitMask l2 = lcbmask2&longlabelmask;
  for(it = lcb.begin();it!=it_end;++it){
    for(it2 = it+1;it2!=it2_end;++it2){
      if(checkLCBOrient(l1,//lcbmask1&longlabelmask,
			l2,//lcbmask2&longlabelmask,
			orientmaskmap[*it]&longlabelmask,
			orientmaskmap[*it2]&longlabelmask)){
      }
      else{
#ifdef DEBUG
	std::cerr << "SAM Mismatch " << *it << "-" << *it2 << std::endl;
#endif	
	return false;
      }
    }
  }
  return true;
}
template<typename TGraph> inline
bool checkLCBOrient(TGraph & g,
		    LCB & lcb,
		    SequenceGenomeMap & sequence2genome){

  BitMask longlabelmask;
  longlabelmask.set();
  return checkLCBOrient(g,lcb,longlabelmask,sequence2genome);
}

template<typename TGraph> inline
bool checkLCBOrient(TGraph & g,
		    LCB & lcb,
		    BitMask & longlabelmask,
		    SequenceGenomeMap & sequence2genome){
  
  LCB::iterator it,it2,it_end,it2_end;

  //Use orientmap so that we can pass by reference using lvalue []
  typename property_map < TGraph, vertex_orient_t >::type orientmap = get(vertex_orient,g);
  typename property_map < TGraph, vertex_orientmask_t >::type orientmaskmap = get(vertex_orientmask,g);
  typename property_map < TGraph, vertex_vlabelmask_t >::type labelmaskmap = get(vertex_vlabelmask,g);
  typename property_map < TGraph, vertex_label_t >::type labelset = get(vertex_label,g);
  typename property_map < TGraph, vertex_len_t >::type lenmap = get(vertex_len,g);
  it_end = lcb.end();
  it2_end = lcb.end();
  for(it = lcb.begin();it!=it_end;++it){
    BitMask o1 = orientmaskmap[*it];
    BitMask l1 = labelmaskmap[*it];
    for(it2 = it+1;it2!=it2_end;++it2){
      BitMask sharedlabels = ((l1&labelmaskmap[*it2])&longlabelmask);
      if(isLabelCollinearMask(sharedlabels,o1,orientmaskmap[*it2])){
      }
      else{
#ifdef DEBUG
	std::cerr << "SAM Mismatch " << *it << "-" << *it2 << std::endl;
#endif
	return false;
      }
    }
  }
  return true;
}

template<typename TGraph, typename TLCBOrientMap> inline
bool checkLCBOrient_old(TGraph & g,
		    LCB & lcb,
		    int lcbidx,
		    TLCBOrientMap & lcborientmap){
  bool shortcircuit=true;
  //Check for label orientation mismatches
  bool mmOrient=false;
  LCB::iterator it,it2,it_end,it2_end;
  EdgeSet mmV;
  //Use orientmap so that we can pass by reference using lvalue []
  typename property_map < TGraph, vertex_orient_t >::type orientmap = get(vertex_orient,g);
  typename property_map < TGraph, vertex_orientmask_t >::type orientmaskmap = get(vertex_orientmask,g);
  typename property_map < TGraph, vertex_vlabelmask_t >::type labelmaskmap = get(vertex_vlabelmask,g);
  it_end = lcb.end();
  it2_end = lcb.end();
  for(it = lcb.begin();it!=it_end;++it){
    BitMask o1 = orientmaskmap[*it];
    BitMask l1 = labelmaskmap[*it];
    for(it2 = it+1;it2!=it2_end;++it2){
      BitMask sharedlabels = (l1&labelmaskmap[*it2]);
#ifdef DEBUG
	std::cerr << "LCBOrientS l1:" << l1 << " l2:" << labelmaskmap[*it2] << std::endl
		  << "LCBOrientS o1:" << o1 << " o2:" << orientmaskmap[*it2] << std::endl;
	std::cerr << "LCBOrientS shared:" << sharedlabels << std::endl;
#endif
      //assert(checkLCBOrient(lcborientmap[lcbidx].second,lcborientmap[lcbidx].second,o1,orientmaskmap[*it2]) 
      //== isLabelCollinearMask(sharedlabels,o1,orientmaskmap[*it2]));
      if(isLabelCollinearMask(sharedlabels,o1,orientmaskmap[*it2])){
#ifdef DEBUG
	  std::cerr << "Match" << std::endl;
#endif
      }
      else{
#ifdef DEBUG
	  std::cerr << "MisMatch" << std::endl;
	  mmV.insert(std::make_pair(*it,*it2));
#endif
	mmOrient=true;
	if(shortcircuit){
	  return !mmOrient;
	}
      }
    }
  }
#ifdef DEBUG
    for(EdgeSet::iterator it = mmV.begin();it!=mmV.end();++it){
      OrientedLabelSet l1 = get(vertex_orient, g, it->first);
      OrientedLabelSet l2 = get(vertex_orient, g, it->second);
      Edge e1;
      bool found;
      
      tie(e1,found) = edge(it->first,it->second,g);
      if(found){
	//assert(get(edge_category,g,e1)==PURPLE);
      }
      tie(e1,found) = edge(it->first,it->second,g);
      if(found){
	//assert(get(edge_category,g,e1)==PURPLE);
      }
      std::cerr << "Orient Mismatch V1:" <<  get(vertex_name,g,it->first) << " L1:";
      printlabel(l1);
      std::cerr <<  " V2:" << get(vertex_name,g,it->second) << " L2:";
      printlabel(l2); 
      std::cerr << std::endl;
    }
#endif
    return !mmOrient;
}

template<typename TLabel>
bool sameLabel(TLabel v1, TLabel v2, TLabel e){
  if(v1==v2){
    if(v1==e){
      return true;
    }
    else{
      return false;
    }
  }
  else{
    return false;
  }
}

template<typename TLabel>
bool sameOrient(TLabel o1, TLabel o2, TLabel v1){
  BitMask revo=o1;
  revo.flip();
  revo=revo&v1;
  //Must logical AND with labelmask v1
  if(o1==o2 || o2==revo){
    return true;
  }
  else{
    return false;
  }
}

//
//TODO
//First determine the orientation of the sequences in the LCB
//Then orient each vertex in the LCB
//assignLCBOrient(LCB)
// BitMask s1omaskrev = s1omask;
//  s1omaskrev.flip();
//  if((s1omask&sharedlabels)==(s2omask&sharedlabels) || 
//     (s1omaskrev&sharedlabels)==(s2omask&sharedlabels)){
//    return true;
//  }
//  else{
//    return false;
//Need to assign LCBs a true orientation for each seq
//This is not as simple as always setting a bit for every + seq
//because blocks within an LCB can be reversed
//Need to mark each block and whether it is flipped wrt to the LCB orientation
//Attempt to mark this flipped state when adding the block to the LCB2
template<typename TGraph>
inline std::pair<BitMask,BitMask> setLCBOrient(TGraph & g, 
					       LCB & lcb, 
					       std::vector<Vertex> & badV,
					       SequenceGenomeMap & sequence2genome){
  BitMask labelmask;
  BitMask orientmask;
  setLCBOrient(g,labelmask,orientmask,lcb,badV,sequence2genome);
  return std::make_pair(labelmask,orientmask);
}

template<typename TGraph>
inline void setLCBOrient(TGraph & g, 
			 BitMask & labelmask, 
			 BitMask & orientmask, 
			 LCB & lcb, 
			 std::vector<Vertex> & badV,
			 SequenceGenomeMap & sequence2genome){
#ifdef DEBUG
    std::cerr << "Setting orientation for lcb with " << lcb.size() << " vertices" << std::endl;
#endif
  badV.clear();
  assert(!orientmask.any());
  assert(!labelmask.any());
  typename property_map < TGraph, vertex_vlabelmask_t >::type lmaskmap = get(vertex_vlabelmask,g);
  typename property_map < TGraph, vertex_orientmask_t >::type omaskmap = get(vertex_orientmask,g);
  typename property_map < TGraph, vertex_len_t >::type lenmap = get(vertex_len,g);
  typename property_map < TGraph, vertex_label_t >::type labelset = get(vertex_label,g);
  //Save label mask and determine evaluation order for orientmask
  std::list<Vertex> sortedV;
  std::map<Label,int> spans;
  std::map<Label,int>::iterator sit;
  //bool found;

  std::map<Vertex,int> vorientmatchcount;
  std::map<Vertex,int>::iterator mit;

  LCB::iterator it,it_end,it2,it2_end;
  it_end = lcb.end();
  it2_end = lcb.end();


  BitMask omask;
  BitMask lmask;
  BitMask sharedlabels;
  //Foreach vertex in the LCB, sum the number of bases with compatible orientation labeling
  for(it = lcb.begin();it!=it_end;++it){
    assert(lmaskmap[*it].any());//at least one seq > MINSPANLEN
    labelmask = labelmask|lmaskmap[*it];
    sortedV.push_back(*it);
    //Check num bp from other vertices compatible within the LCBs 
    omask = omaskmap[*it];
    lmask = lmaskmap[*it];
    int len = lenmap[*it];
    mit = vorientmatchcount.find(*it);
    for(it2 = lcb.begin();it2!=it2_end;++it2){
      if(it2!=it){
	sharedlabels = (lmask&lmaskmap[*it2]);
	if(isLabelCollinearMask(sharedlabels,omask,omaskmap[*it2])){
	  //Update count of bp
	  if(mit==vorientmatchcount.end()){
	    vorientmatchcount[*it]=len;
	    mit = vorientmatchcount.find(*it);
	  }
	  else{
	    mit->second = mit->second+len;
	  }
	}
      }
    }
  }
  //Determine orientation compatible with most bp in the LCB
  //
  //Sort by vorientmatchcount so that we consider most frequent orientations first
  //before alternative orientations
  //TODO, confirm this is actually working
  sortedV.sort(matchmaporder(&vorientmatchcount));
#ifdef DEBUG
    std::cerr << "LCB labelmask :" << labelmask << std::endl;
#endif
  BitMask currlabelmask;

  //Realize speedup if < 64 genomes
#if BITMAX > 64
  boost::unordered_set<BitMask,bitsethasher_string> altOrients;
#else
  boost::unordered_set<BitMask,bitsethasher_ulong> altOrients;
#endif

  std::set<Vertex> badVS;
  std::list<Vertex>::iterator svit,svit_end;
  svit_end = sortedV.end();
  for(svit=sortedV.begin();svit!=svit_end;++svit){
    if((omaskmap[*svit]&currlabelmask)  //Vertex orients for seqs in LCB
       ==(orientmask&lmaskmap[*svit])){ //LCB orients on current vertex
      orientmask = orientmask|omaskmap[*svit];
    }
    else{
      BitMask s1omaskrev = omaskmap[*svit];
      s1omaskrev.flip();
      s1omaskrev = s1omaskrev&lmaskmap[*svit];
      if((s1omaskrev&currlabelmask)
	 ==(orientmask&lmaskmap[*svit])){
	orientmask = orientmask|s1omaskrev;
      }
      else{
	//Alternative orientation
	//Sorting on bp above ensures this alternative is congruent with 
	//fewer bp than at least on other alternative orientation for the LCB
	badV.push_back(*svit);
#ifdef DEBUG
	  badVS.insert(*svit);
#endif
	altOrients.insert((omaskmap[*svit]&labelmask));
	
	//(1)Simple breakpoints
	//A simple breakpoint occues when a pair of vertices in the block have incompatible orientations
	//For example, consider an LCB with vertices 1(a+b+) and 2(a+b-)

	//(2)Compound orientation breakpoints
	//A compound breakpoint occurs when a combination of
	//vertices has an orientation that is incompatible with 
	//other vertices in a block
	//For example, consider an LCB with vertices 1(a+b+,d-) 2(b+c-) 3(a+b+c-) 4(c-,d+)
	//All pairwise comparisons are congruent in orientation if reversals are allowed
	//But not all multiway combinations are congruent because 
	//      1+2 requires (a+b+c-d-)
	//while 3+4 requires (a+b+c-d+)
	//producing an incompatibility in the orientations of d
	//when considered in the context of the other
	//vertices in the block

	//Save a map v->totallcbsize congruent with orientmask(v)
	//Sort by totallcbsize in decreasing order
	//Build LCB mask in sorted order 
	//This ensures orientations congruent with the most bp are considered first
	//Each alternative orientation encountered is a compound breakpoint

	//Count length/number of vertices congruent with this alternative orientation
#ifdef DEBUG
	  std::cerr << "Mismatched orient for vertex " << get(vertex_name,g,*svit) << " with vertex/lcb/shared labels " <<std::endl 
		    << " vertex_l     :" << lmaskmap[*svit] << std::endl
		    << " lcb_l        :" << currlabelmask << std::endl
		    << " shared_l     :" << (lmaskmap[*svit]&currlabelmask) << std::endl
		    << " lcb_o        :" << orientmask << std::endl
		    << " vertex_o     :" << (omaskmap[*svit]&currlabelmask) << std::endl
		    << " vertex_ro    :" << (s1omaskrev&currlabelmask) << std::endl
		    << " shared_o     :" << (orientmask&lmaskmap[*svit]) << std::endl;
	//assert(false); 
	  std::cerr << "Num bp matching orient len:" << vorientmatchcount[*svit] << std::endl;
#endif
      }
    }
    //Update label mask
    currlabelmask = currlabelmask|lmaskmap[*svit];
  }
#ifdef DEBUG
    std::cerr << "LCB orientmask:" << orientmask << std::endl;
#endif
  //Check that orientmask is only set on member sequences in the LCB
  assert((orientmask&labelmask)==orientmask);
  /*
  //TESTING
  if(0 && DEBUG){
    for(it = lcb.begin();it!=it_end;++it){
      //Label for vertex contains strict subset of seqs as LCB
      assert(((lmaskmap[*it]&labelmask))==lmaskmap[*it]);
      assert(((omaskmap[*it]&labelmask))==((omaskmap[*it]&lmaskmap[*it]&labelmask)));
      assert(((omaskmap[*it]&labelmask))==((omaskmap[*it]&lmaskmap[*it])));
      
      //Check that vertex omask matches lcb omask for member sequences on the vertex
      if((omaskmap[*it]&labelmask)==(orientmask&lmaskmap[*it])){
	if(badVS.find(*it)!=badVS.end()){
	  assert(false);
	}
      }
      else{
	BitMask s1omaskrev = omaskmap[*it];
	s1omaskrev.flip();
	s1omaskrev = s1omaskrev&lmaskmap[*it];
	if((s1omaskrev&labelmask)==(orientmask&lmaskmap[*it])){
	  if(badVS.find(*it)!=badVS.end()){
	    assert(false);
	  }
	}
	else{
	  if(badVS.find(*it)==badVS.end()){
	    std::cerr << "Vertex " << get(vertex_name,g,*it) << " does not match LCB orientation " << std::endl
		      << "vertex_o: " << (omaskmap[*it]&labelmask) << std::endl
		      << "vertex_ro:" << (s1omaskrev&labelmask) << std::endl
		      << "lcb_o:    " << (orientmask&lmaskmap[*it]) << std::endl;
	    
	    assert(false);
	  }
	}
      }
    }
  }
  */
  if(badV.size()>0){
#ifdef DEBUG
    std::cerr << "LCB has " << badV.size()  << "/" << lcb.size() 
	      << " misoriented vertices. Max num alternative orients " << altOrients.size() << std::endl;
#endif
  }
#ifdef DEBUG
    std::cerr << "Setting orientation done" << std::endl;
#endif
}

template<typename TGraph, typename TLCBMap, typename TComponentMap>
void setLCBOrient(TGraph & g, TLCBMap & lcborientmap, TComponentMap & componentMap, SequenceGenomeMap & sequence2genome){
#ifdef DEBUG
    std::cerr <<"Resetting lcborientmap" << std::endl;
#endif
  lcborientmap.clear();
  typename property_map < TGraph, vertex_vlabelmask_t >::type lmaskmap = get(vertex_vlabelmask,g);
  typename property_map < TGraph, vertex_orientmask_t >::type omaskmap = get(vertex_orientmask,g);

  int lcbidx=0;
  std::vector<Vertex> badV;
  typename TComponentMap::iterator lit,lit_end;
  lit_end=componentMap.end();
  for(lit = componentMap.begin();lit!=lit_end;++lit){
    lcborientmap[lcbidx] = setLCBOrient(g,*lit,badV,sequence2genome);
    //At least one genome label must be set
    assert(lcborientmap[lcbidx].first.any());
    lcbidx++;
  }
}

template<typename TLCB,
	 typename TLenMap,
	 typename TLabelMap>
BitMask setSpanMask(TLCB & lcb,
		    TLenMap & lenmap,
		    TLabelMap & labelset,
		    SequenceGenomeMap & sequence2genome,
		    int MINSPANLEN=0){
   
  std::map<Label,int> spans;
  std::map<Label,int>::iterator sit;
  LCB::iterator it,it_end,it2,it2_end;
  LabelSet::iterator lit,lit_end;
  BitMask longlabelmask;
  it_end = lcb.end();
  it2_end = lcb.end();

  //Only consider sequences with > MINSPANLEN
  for(it = lcb.begin();it!=it_end;++it){
    //assert(labelset.find(*it)!=labelset.end());
    //assert(lenmap.find(*it)!=lenmap.end());
    lit_end = labelset[*it].end();
    for(lit = labelset[*it].begin();lit != lit_end;++lit){
      assert(sequence2genome.find(*lit)!=sequence2genome.end());
      Label l = sequence2genome[*lit]; 
      sit = spans.find(l);
      if(sit==spans.end()){
	spans.insert(std::make_pair(l,lenmap[*it]));
      }
      else{
	int prev = sit->second;
	sit->second = prev + lenmap[*it];
	assert(spans[l] == (prev + lenmap[*it]));
      }
    }
  }
  for(std::map<Label,int>::iterator sit = spans.begin();sit!=spans.end();++sit){
    if(sit->second>(int)MINSPANLEN){
      longlabelmask.set(sit->first);
    }
  }
  return longlabelmask;
}



//####################
//LCB reporting functions
//

unsigned int get_LCB_length(LCB & lcb, 
			    property_map < LCBSynFilterGraph, vertex_orient_t >::type orientmap,
			    property_map < LCBSynFilterGraph, vertex_len_t >::type lenmap,
			    VertexLabelIntervalMap & coordinates,
			    LCBLabelIntervalMap & lcbcoords,
			    unsigned int lcbidx, 
			    int & totallen,
			    SequenceGenomeMap & sequence2genome,
			    int minlength=0
			    ){
  std::map<unsigned int,unsigned int> mincoordsbyseq;
  std::map<unsigned int,unsigned int> maxcoordsbyseq;
  std::map<Label,int> spans;
  OrientedLabelSet alllabel;
  LCB::iterator it;
  for(it = lcb.begin();it!=lcb.end();++it){
    OrientedLabelSet::iterator it2_end = orientmap[*it].end();
    for(OrientedLabelSet::iterator it2 = orientmap[*it].begin();it2!=it2_end;++it2){
      assert(it2->first>=0);
      Label seqidx = it2->first;
      Label genomeidx = sequence2genome[seqidx];
      if(spans.find(genomeidx)==spans.end()){
	spans.insert(std::make_pair(genomeidx,lenmap[*it]));
      }
      else{
	spans[genomeidx] = spans[genomeidx] + lenmap[*it];
      }
    }
  }
  for(it = lcb.begin();it!=lcb.end();++it){
    OrientedLabelSet::iterator it2_end = orientmap[*it].end();
    for(OrientedLabelSet::iterator it2 = orientmap[*it].begin();it2!=it2_end;++it2){
      assert(it2->first>=0);
      Label seqidx = it2->first;
      Label genomeidx = sequence2genome[seqidx];
      assert(spans.find(genomeidx)!=spans.end());
      if(spans[genomeidx] >= minlength){
	if(mincoordsbyseq.find(seqidx)==mincoordsbyseq.end()){
	  mincoordsbyseq[seqidx] = std::numeric_limits<unsigned int>::max();
	}
	if(maxcoordsbyseq.find(seqidx)==maxcoordsbyseq.end()){
	  maxcoordsbyseq[seqidx] = std::numeric_limits<unsigned int>::min();
	}
	
	alllabel.insert(OrientedLabel(seqidx,true));
	VertexLabelIntervalMap::iterator vit =  coordinates.find(std::make_pair(*it,seqidx));
	assert(vit!=coordinates.end());
	assert(vit->second.second>vit->second.first);
	totallen+=abs(vit->second.second-vit->second.first);	
	mincoordsbyseq[seqidx] = ((unsigned int)vit->second.first<mincoordsbyseq[seqidx]) ? vit->second.first : mincoordsbyseq[seqidx];
	maxcoordsbyseq[seqidx] = ((unsigned int)vit->second.second>maxcoordsbyseq[seqidx]) ? vit->second.second : maxcoordsbyseq[seqidx];
	assert(mincoordsbyseq[seqidx] != std::numeric_limits<unsigned int>::max());
	assert(maxcoordsbyseq[seqidx] != std::numeric_limits<unsigned int>::min());
	//std::cerr << "seq:" << seqidx << " " << *it << " " << maxcoordsbyseq[seqidx] << " " <<  mincoordsbyseq[seqidx] << std::endl;
      }
    }
  }
  if(alllabel.size()==0){
    return 0;
  }
  //Save spanning coords for lcb in lcbcoords
  OrientedLabelSet::iterator it2_end = alllabel.end();
  for(OrientedLabelSet::iterator it2 = alllabel.begin();it2!=it2_end;++it2){
    Label seqidx = it2->first;
    Label genomeidx = sequence2genome[seqidx];
    assert(spans[genomeidx] >= minlength);

    LCBLabelIntervalMap::iterator it=lcbcoords.find(std::make_pair(lcbidx,seqidx));
    if(it!=lcbcoords.end()){
      lcbcoords.erase(it);
    }
    assert(mincoordsbyseq[seqidx]<maxcoordsbyseq[seqidx]);
    lcbcoords.insert(std::make_pair(std::make_pair(lcbidx,seqidx),std::make_pair(mincoordsbyseq[seqidx],maxcoordsbyseq[seqidx])));
  }
  //Maximum span of all seqs
  unsigned int maxminlen = 0;
  OrientedLabelSet::iterator it_end = alllabel.end();
  for(OrientedLabelSet::iterator it = alllabel.begin();it!=it_end;++it){
    Label seqidx = it->first;
    Label genomeidx = sequence2genome[seqidx];
    assert(spans[genomeidx] >= minlength);
    assert(maxcoordsbyseq[seqidx]>mincoordsbyseq[seqidx]);
    if(maxcoordsbyseq[seqidx]>mincoordsbyseq[seqidx]){
      unsigned int len = (unsigned int)(maxcoordsbyseq[seqidx] - mincoordsbyseq[seqidx]);
      //std::cerr << "seq:" << seqidx << " " << maxcoordsbyseq[seqidx] << " " <<  mincoordsbyseq[seqidx] << " " << len << std::endl;
      maxminlen = len>maxminlen ? len : maxminlen;
    }
    else{
      //std::cerr << "??seq:" << seqidx << " " << maxcoordsbyseq[seqidx] << " " <<  mincoordsbyseq[seqidx] << std::endl;
      assert(false);
    }
  }
  assert(maxminlen>0);
  //std::cerr << maxminlen << std::endl;
  return maxminlen;
}


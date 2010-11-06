 /*==========================================================================
                SeqAn - The Library for Sequence Analysis
                          http://www.seqan.de 
 ============================================================================
  Copyright (C) 2007

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.

 ============================================================================
  $Id: graph_align_tcoffee_io.h 1919 2008-05-02 15:54:46Z rausch@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_GRAPH_ALIGN_TCOFFEE_IO_H
#define SEQAN_HEADER_GRAPH_ALIGN_TCOFFEE_IO_H

#include <assert.h>

namespace SEQAN_NAMESPACE_MAIN
{

/////////////////////////////////////////////////////////////////////////////
// Input and Output of an alignment graph, Tree Reading im Newick Format
//////////////////////////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////////////////////////
// Tags
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

/**
.Tag.Alignment Graph Format.value.TCoffeeLib:
	T-Coffee library format to read and write an alignment graph.
*/

struct TCoffeeLib_;
typedef Tag<TCoffeeLib_> const TCoffeeLib;

//////////////////////////////////////////////////////////////////////////////

/**
.Tag.Alignment Graph Format.value.BlastLib:
	A blast library for matches for an alignment graph.
*/

struct BlastLib_;
typedef Tag<BlastLib_> const BlastLib;

//////////////////////////////////////////////////////////////////////////////

/**
.Tag.Alignment Graph Format.value.MummerLib:
	A mummer library for matches for an alignment graph.
*/

struct MummerLib_;
typedef Tag<MummerLib_> const MummerLib;


//////////////////////////////////////////////////////////////////////////////

/**
.Tag.Alignment Graph Format.value.NewickFormat:
	NewickFormat format to write a guide tree.
*/

struct NewickFormat_;
typedef Tag<NewickFormat_> const NewickFormat;



/////////////////////////////////////////////////////////////////////////////
// T-Coffee Library Reading / Writing
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////



template<typename TFile, typename TFragment, typename TSpec, typename TScoreValue, typename TSpec2, typename TNames>
void 
read(TFile & file,
	 String<TFragment, TSpec>& matches,
	 String<TScoreValue, TSpec2>& scores,
	 TNames const& names,
	 TCoffeeLib) 
{
	SEQAN_CHECKPOINT
	typedef typename Size<TFragment>::Type TSize;
	typedef typename Value<TFile>::Type TValue;
	typedef typename Value<TNames>::Type TName;
	
	
	// Map the names to slots
	typedef std::map<TName, TSize> TNameToPosition;
	TNameToPosition namePosMap;
	for(TSize i = 0;i<length(names);++i) namePosMap.insert(std::make_pair(names[i], i));
	
	// Remember the correct spots
	TValue c;
	if (_streamEOF(file)) return;
	else c = _streamGet(file);
	_parse_skipLine(file, c);
	TSize nseq = (TSize) _parse_readNumber(file, c);
	_parse_skipLine(file, c);

	// Read sequences
	typedef String<TSize> TMapping;
	TMapping posMap;
	resize(posMap, nseq);
	for(TSize i=0; i<nseq; ++i) {
		TName myName;
		_parse_readIdentifier(file, myName, c);
		value(posMap, i) = namePosMap.find(myName)->second;
		_parse_skipLine(file, c);
	}

	bool seq1ToN = false;
	if (_streamEOF(file)) return;
	
	typedef std::pair<std::pair<TSize, TSize>, TScoreValue> TResiduePair;
	typedef std::set<TResiduePair> TResiduePairSet;
	String<TResiduePairSet> resPair;
	resize(resPair, nseq * nseq);	
	TSize seq1 = 0;
	TSize seq2 = 0;
	bool firstPass = true;
	while (!_streamEOF(file)) {
		_parse_skipWhitespace(file,c);
		if (_streamEOF(file)) break;
		if (c == '#') {
			c = _streamGet(file);
			_parse_skipWhitespace(file,c);
			seq1 = _parse_readNumber(file, c);
			_parse_skipWhitespace(file,c);
			seq2 = _parse_readNumber(file, c);
			if (firstPass) {
				firstPass = false;
				if ((seq1 != 0) && (seq2 != 0)) seq1ToN = true;
			}
			if (seq1ToN) {
				--seq1;
				--seq2;
			}
			seq1 = value(posMap, seq1);
			seq2 = value(posMap, seq2);
		} else if (c == '!') {
			_parse_skipLine(file, c);
		} else {
			TSize res1 = _parse_readNumber(file, c);
			_parse_skipWhitespace(file,c);
			TSize res2 = _parse_readNumber(file, c);
			_parse_skipWhitespace(file,c);
			TSize weight = _parse_readNumber(file, c);
			_parse_skipLine(file,c);

			if (seq1 < seq2) {
				TSize index = seq1 * nseq + seq2;
				resPair[index].insert(std::make_pair(std::make_pair(--res1,--res2), weight));
			} else {
				TSize index = seq2 * nseq + seq1;
				resPair[index].insert(std::make_pair(std::make_pair(--res2,--res1), weight));
			}		
		}
	}
	for(TSize i = 0; i<length(resPair); ++i) {
		if (resPair[i].empty()) continue;
		TSize seq1 = i / nseq;
		TSize seq2 = i % nseq;
		//std::cout << "#" << seq1 << ',' << seq2 << std::endl;
		typename TResiduePairSet::const_iterator pos = resPair[i].begin();
		typename TResiduePairSet::const_iterator posEnd = resPair[i].end();
		TSize startMatch1 = pos->first.first;
		TSize startMatch2 = pos->first.second;
		TScoreValue carg = pos->second;
		TSize len = 1;
		//std::cout << pos->first.first << ',' << pos->first.second << ',' << pos->second << std::endl;
		++pos;
		while(pos != posEnd) {
			if ((startMatch1 + len == pos->first.first) &&
				(startMatch2 + len == pos->first.second) &&
				(carg / (TScoreValue) len == pos->second)) {
					carg += pos->second;
					++len;
			} else {
			  //SVA
			  //appendValue(matches, TFragment(seq1, startMatch1, seq2, startMatch2, len,len));
			  appendValue(matches, TFragment(seq1, startMatch1, seq2, startMatch2, len));
				appendValue(scores, carg);
				startMatch1 = pos->first.first;
				startMatch2 = pos->first.second;
				carg = pos->second;
				len = 1;
			}
			//std::cout << pos->first.first << ',' << pos->first.second << ',' << pos->second << std::endl;
			++pos;
		}
		//SVA
		//appendValue(matches, TFragment(seq1, startMatch1, seq2, startMatch2, len,len));
		appendValue(matches, TFragment(seq1, startMatch1, seq2, startMatch2, len));
		appendValue(scores, carg);
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TString, typename TSpec, typename TNames>
void 
read(TFile & file,
	 StringSet<TString, TSpec>& oriStr,
	 TNames& names,
	 TCoffeeLib) 
{
	SEQAN_CHECKPOINT
	typedef typename Size<TNames>::Type TSize;
	typedef typename Position<TFile>::Type TPosition;
	typedef typename Value<TFile>::Type TValue;
	
	TValue c;
	if (_streamEOF(file)) return;
	else c = _streamGet(file);

	// Ignore first line
	_parse_skipLine(file, c);
	
	// Read number of sequences
	TSize nSeq = (TSize) _parse_readNumber(file, c);
	resize(oriStr, nSeq);
	_parse_skipLine(file, c);

	// Read sequences
	for(TSize i=0; i<nSeq; ++i) {
		appendValue(names, _parse_readIdentifier(file, c));
		_parse_skipWhitespace(file, c);
		_parse_readNumber(file, c);
		_parse_skipWhitespace(file, c);
		_parse_readSequenceData(file,c,oriStr[i]);
		_parse_skipLine(file, c);
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TStringSet, typename TCargo, typename TSpec, typename TNames>
void write(TFile & file, 
		   Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
		   TNames& names,
		   TCoffeeLib) 
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename Value<TStringSet>::Type TString;
	typedef typename Size<TStringSet>::Type TSize;
	
	typedef std::pair<std::pair<TSize, TSize>, TCargo> TResiduePair;
	typedef std::set<TResiduePair> TResiduePairSet;
	TSize nseq = length(stringSet(g));
	String<TResiduePairSet> resPair;
	resize(resPair, nseq * nseq);
	
	typedef typename Iterator<TGraph, EdgeIterator>::Type TIter;
	TIter it(g);
	for(;!atEnd(it);++it) {
		TVertexDescriptor sV = sourceVertex(it);
		TVertexDescriptor tV = targetVertex(it);
		TSize fragLen = fragmentLength(g,sV);
		TSize fragPos1 = fragmentBegin(g,sV);
		TSize fragPos2 = fragmentBegin(g,tV);
		TSize seq1 = sequenceId(g,sV);
		TSize seq2 = sequenceId(g,tV);
		TCargo my_carg =  getCargo(*it);
		if (my_carg <= 0) my_carg = 1;
		else my_carg = (TCargo) ((double) my_carg / (double) fragLen);
		for(TSize i = 0; i<fragLen; ++i) {
			resPair[seq1 * nseq + seq2].insert(std::make_pair(std::make_pair(fragPos1 + i, fragPos2 + i), my_carg));
		}
	}

	_streamWrite(file, "! TC_LIB_FORMAT_01\n");
	TSize len = length(getStringSet(g));
	_streamPutInt(file, len);
	_streamPut(file, '\n');
	for(TSize i=0;i<len;++i) {
		_streamWrite(file, names[i]);
		_streamPut(file, ' ');
		TString str = value(getStringSet(g), i);
		_streamPutInt(file, length(str));
		_streamPut(file, ' ');
		_streamWrite(file, str);
		_streamPut(file, '\n');
	}

	for(TSize i = 0; i<length(resPair); ++i) {
		if (resPair[i].empty()) continue;
		TSize seq1 = i / nseq;
		TSize seq2 = i % nseq;
		_streamPut(file, '#');
		_streamPutInt(file, seq1 + 1);
		_streamPut(file, ' ');
		_streamPutInt(file, seq2 + 1);
		_streamPut(file, '\n');	
		typename TResiduePairSet::const_iterator pos = resPair[i].begin();
		typename TResiduePairSet::const_iterator posEnd = resPair[i].end();
		while(pos != posEnd) {
			_streamPutInt(file, pos->first.first + 1);
			_streamPut(file, ' ');
			_streamPutInt(file, pos->first.second + 1);
			_streamPut(file, ' ');
			_streamPutInt(file, pos->second);
			_streamPut(file, '\n');	
			++pos;
		}
	}
	_streamWrite(file, "! SEQ_1_TO_N");
	_streamPut(file, '\n');
}


/////////////////////////////////////////////////////////////////////////////
// FastaAlign Reading
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TString, typename TSpec, typename TNames>
void 
read(TFile & file,
	 StringSet<TString, TSpec>& oriStr,
	 TNames& names,
	 FastaAlign) 
{
	SEQAN_CHECKPOINT
	typedef typename Size<TNames>::Type TSize;
	typedef TSize TWord;
	typedef typename Value<TFile>::Type TValue;
	
	TValue c;
	if (_streamEOF(file)) return;
	else c = _streamGet(file);

	// Read sequences
	TString seq;
	while(!_streamEOF(file)) {
		_parse_skipWhitespace(file, c);
		if (_streamEOF(file)) break;
		if (c == '>') {
			if (length(seq)) {
				appendValue(oriStr, seq);
				clear(seq);
			}
			c = _streamGet(file);
			appendValue(names, _parse_readIdentifier(file, c));
			_parse_skipLine(file, c);
		} else if ((c == '-') || (c == '.') || (c == '\n') || (c == '\r')) {
			c = _streamGet(file);
		} else {
			appendValue(seq, c);
			c = _streamGet(file);
		}
	}
	if (length(seq)) appendValue(oriStr, seq);
}


template<typename TValue, typename TSpec2, typename TFragment, typename TSpec, typename TScores, typename TSize>
void 
_collectSegmentMatches(String<TValue, TSpec2> const& mat,
		       String<TFragment, TSpec>& matches,
		       TScores& scores,
		       TSize nseq) 
{
	SEQAN_CHECKPOINT
	TSize len = length(mat) / nseq;
	TValue gapChar = gapValue<TValue>();

	// Create the anchor graph
	typedef String<TFragment, TSpec> TFragmentString;
	typedef typename Iterator<TFragmentString>::Type TFragmentStringIter;
	typedef std::pair<TSize, TSize> TResiduePair;
	typedef std::set<TResiduePair> TResiduePairSet;
	String<TResiduePairSet> resPair;
	resize(resPair, nseq * nseq);	
	//SVA accumulate all the runs of ungapped columns
	for(TSize seq1 = 0; seq1 < nseq - 1; ++seq1) {
		for(TSize seq2 = seq1 + 1; seq2 < nseq; ++seq2) {
			TSize index = seq1 * nseq + seq2;
			TSize offset1 = 0;
			TSize offset2 = 0;
			for(TSize col = 0; col<len; ++col) {
				if (value(mat, seq1 * len + col) != gapChar) {
					if (value(mat, seq2 * len + col) != gapChar) {
						resPair[index].insert(std::make_pair(offset1, offset2));
						++offset1;
						++offset2;
					} else ++offset1;
				} else if (value(mat, seq2 * len + col) != gapChar) ++offset2;
			}
		}
	}
	//SVA append runs to matches
	for(TSize i = 0; i<length(resPair); ++i) {
		if (resPair[i].empty()) continue;
		TSize seq1 = i / nseq;
		TSize seq2 = i % nseq;
		typename TResiduePairSet::const_iterator pos = resPair[i].begin();
		typename TResiduePairSet::const_iterator posEnd = resPair[i].end();
		TSize startMatch1 = pos->first;
		TSize startMatch2 = pos->second;
		TSize len = 1;
		++pos;
		while(pos != posEnd) {
			if ((startMatch1 + len == pos->first) && (startMatch2 + len == pos->second)) ++len;
			else {
			  //SVA
			  //appendValue(matches, TFragment(seq1, startMatch1, seq2, startMatch2, len,len));
			  appendValue(matches, TFragment(seq1, startMatch1, seq2, startMatch2, len));
			  appendValue(scores, len);
			  startMatch1 = pos->first;
			  startMatch2 = pos->second;
			  len = 1;
			}
			++pos;
		}
		//SVA
		//appendValue(matches, TFragment(seq1, startMatch1, seq2, startMatch2, len, len));
		appendValue(matches, TFragment(seq1, startMatch1, seq2, startMatch2, len));
		appendValue(scores, len);
	}
}

template<typename TValue, typename TSpec2, typename TFragment, typename TSpec, typename TScores, typename TSize, typename TStringSet>
  void 
  _collectSegmentMatchesOffset(String<TValue, TSpec2> const& mat,
			       String<TFragment, TSpec>& matches,
			       TScores& scores,
			       TSize nseq,
			       String<TSize>& coords,
			       String<bool>& orients,
			       TStringSet const& strSet,
			       TSize inseq1,
			       TSize inseq2
			       )
{
  SEQAN_CHECKPOINT
    
  TSize flen = length(mat) / nseq;
  TValue gapChar = gapValue<TValue>();
  
  // Create the anchor graph
  typedef String<TFragment, TSpec> TFragmentString;
  typedef typename Iterator<TFragmentString>::Type TFragmentStringIter;
  typedef std::pair<TSize, TSize> TResiduePair;
  typedef std::set<TResiduePair> TResiduePairSet;
  String<TResiduePairSet> resPair;
  resize(resPair, 1);//nseq * nseq);	
  
  //SVA accumulate all the runs of ungapped columns
  /*
  for(TSize seq1 = 0; seq1 < nseq; ++seq1) {
    if(seq1 == inseq1){// || seq1 == inseq2){//hack for pairwise only
      //for(TSize seq2 = seq1 + 1; seq2 < nseq; ++seq2) {
      for(TSize seq2 = 0; seq2 < nseq; ++seq2) {
	if(seq2 == inseq2){// || seq2 == inseq1){//hack for pairwise only
	  TSize index = seq1 * nseq + seq2;
	  TSize offset1 = 0;
	  TSize offset2 = 0;
	  if(orients[inseq2])
	    ;//offset2 = length(value(strSet, inseq2))-1;
	  for(TSize col = 0; col<len; ++col) {
	    if (value(mat, seq1 * len + col) != gapChar) {
	      if (value(mat, seq2 * len + col) != gapChar) {
		resPair[index].insert(std::make_pair(offset1, offset2));
		++offset1;
		if(orients[inseq2])
		  ++offset2;//--offset2;
		else
		  ++offset2;
	      } else ++offset1;
	    } else if (value(mat, seq2 * len + col) != gapChar) 
	      if(orients[inseq2])
		++offset2;//--offset2;
	      else
		++offset2;
	  }
	}
      }
    }
  }
  */
  TSize index = 0;//inseq1 * nseq + inseq2;
  TSize offset1 = 0;
  TSize offset2 = 0;
  for(TSize col = 0; col<flen; ++col) {
    if (value(mat, inseq1 * flen + col) != gapChar) {
      if (value(mat, inseq2 * flen + col) != gapChar) {
	resPair[index].insert(std::make_pair(offset1, offset2));
	++offset1;
	if(orients[inseq2])
	  ++offset2;//--offset2;
	else
	  ++offset2;
      } else ++offset1;
    } else if (value(mat, inseq2 * flen + col) != gapChar) 
      if(orients[inseq2])
	++offset2;//--offset2;
      else
	++offset2;
  }
  /*
  //SVA append runs to matches
  for(TSize i = 0; i<length(resPair); ++i) {
    if (resPair[i].empty()) continue;
    TSize seq1 = i / nseq;
    TSize seq2 = i % nseq;
    if(seq1 != inseq1){// && seq1 != inseq2){
      std::cerr << "Bad input seqs" << std::endl;
      exit(1);
    }
    if(seq2 != inseq2){// && seq2 != inseq2){
      std::cerr << "Bad input seqs" << std::endl;
      exit(1);
    }
  }
  */ 
  //for(TSize i = 0; i<length(resPair); ++i) {
  //if (resPair[i].empty()) continue;

  //TSize seq1 = inseq1;
  //TSize seq2 = inseq2;
  typename TResiduePairSet::const_iterator pos = resPair[0].begin();
  typename TResiduePairSet::const_iterator posEnd = resPair[0].end();
  TSize startMatch1 = pos->first;
  TSize startMatch2 = pos->second;
  //Should be zero 
  if(startMatch1 !=0 && startMatch2 != 0){
    std::cerr << startMatch1 << ":" << startMatch2 << std::endl;
    exit(1);
  }
  TSize len = 1;
  ++pos;
  while(pos != posEnd) {
    if ((startMatch1 + len == pos->first) && (startMatch2 + len == pos->second)) ++len;
    else{
      //std::cout << "Len: " << len << " " << startMatch1 + len << "==" << pos->first << " " 
      //<< startMatch2 + len << "==" <<  pos->second << std::endl;
      //std::cout << "S1: " << inseq1 << " " << length(value(strSet, inseq1)) << " " << coords[inseq1] << ":" << startMatch1+coords[inseq1] 
      //<< " o:" << orients[inseq1] << " len:" << len << std::endl;
      //std::cout << "S2: " << inseq2 << " " << length(value(strSet, inseq2)) << " " << coords[inseq2] << ":" << startMatch2+coords[inseq2] 
      //<< " o:" << orients[inseq2] << " len:" << len << std::endl;
      bool rev=false;
      if(orients[inseq1] == orients[inseq2]){
	if(orients[inseq1]){
	  std::cerr << "Bad relative orientation for seq1,seq2 "<< orients[inseq1] << ":" << orients[inseq2] << std::endl;
	  exit(1);
	}
      }
      else{
	if(orients[inseq1]){
	  rev = true;
	  //currently die, alternatively swap inseq1 and inseq2
	  std::cerr << "Unsupported orientation '-' for seq1" << std::endl;
	  exit(1);
	}
	else{
	  if(orients[inseq2]){
	    rev = true;
	  }
	}
      } 
      typedef Fragment<TSize, ExactReversableFragment<TSpec> > TFragmentRev;
      if(rev){
#ifdef DEBUGGING
	std::cout << "REVMATCH" << std::endl;
	std::cout << "S1: " << inseq1 << " " << length(value(strSet, inseq1)) << " " << coords[inseq1] << ":" << startMatch1+coords[inseq1] 
		  << " o:" << orients[inseq1] << " len:" << len << std::endl;
	std::cout << "S2: " << inseq2 << " " << length(value(strSet, inseq2)) 
		  << " startMatch2:" << startMatch2 << " " << coords[inseq2] << ":" << length(value(strSet, inseq2))-(coords[inseq2]+len)-startMatch2
		  << " o:" << orients[inseq2] << " len:" << len << std::endl;
#endif
	//__appendNewMatch(matches, scores, inseq1, inseq2, 
	//startMatch1+coords[inseq1], (length(value(strSet, inseq2)) - (startMatch2+coords[inseq2]+len)),len,rev);
	//Gap placed exactly at the wrong end
	//__appendNewMatch(matches, scores, inseq1, inseq2, 
	//startMatch1+coords[inseq1], length(value(strSet, inseq2))-(coords[inseq2]-startMatch2), len,rev);
	//__appendNewMatch(matches, scores, inseq1, inseq2, 
	//startMatch1+coords[inseq1], coords[inseq2]-len-startMatch2, len,rev);
#ifdef REVCOORD
	__appendNewMatch(matches, scores, inseq1, inseq2, 
			 startMatch1+coords[inseq1], coords[inseq2]+startMatch2, len,rev);
#else
	assert(coords[inseq1]<length(value(strSet, inseq1)));
	assert(coords[inseq2]<length(value(strSet, inseq2)));
	assert((int)coords[inseq1]>=0);
	assert((int)coords[inseq2]>=0);
	assert((int)(startMatch1+coords[inseq1])>=0);
	assert((int)(startMatch2+coords[inseq2])>=0);
	assert((int)(startMatch1+coords[inseq1]+len)<=(int)length(value(strSet, inseq1)));
	assert(length(value(strSet, inseq2))-coords[inseq2]-startMatch2-len<=length(value(strSet, inseq2)));
	__appendNewMatch(matches, scores, inseq1, inseq2, 
			 startMatch1+coords[inseq1], length(value(strSet, inseq2))-coords[inseq2]-startMatch2-len, len,rev);
	
#endif
      }
      else{
	assert(startMatch1+coords[inseq1]<length(value(strSet, inseq1)));
	assert(startMatch2+coords[inseq2]<length(value(strSet, inseq2)));
	//assert(startMatch1+coords[inseq1]>=0);
	//assert(startMatch2+coords[inseq2]>=0);
	assert(startMatch1+coords[inseq1]+len<=length(value(strSet, inseq1)));
	assert(startMatch2+coords[inseq2]+len<=length(value(strSet, inseq2)));
	__appendNewMatch(matches, scores, inseq1, inseq2, 
			 startMatch1+coords[inseq1], startMatch2+coords[inseq2], len,rev);
      }
      startMatch1 = pos->first;
      startMatch2 = pos->second;
      len = 1;
    }
    ++pos;
  }
  //std::cout << "Len: " << len << " startmatch1:" << startMatch1 << " startmatch2:"  << startMatch2 << std::endl;
  //std::cout << "S1: " << inseq1 << " " << length(value(strSet, inseq1)) << " " << coords[inseq1] << ":" << startMatch1+coords[inseq1] 
  //<< " o:" << orients[inseq1] << " len:" << len << std::endl;
  //std::cout << "S2: " << inseq2 << " " << length(value(strSet, inseq2)) << " " << coords[inseq2] << ":" << startMatch2+coords[inseq2] 
  //<< " o:" << orients[inseq2] << " len:" << len << std::endl;
  
  bool rev=false;
  if(orients[inseq1] == orients[inseq2]){
    if(orients[inseq1]){
      std::cerr << "Bad relative orientation for seq1,seq2 "<< orients[inseq1] << ":" << orients[inseq2] << std::endl;
      exit(1);
    }
  }
  else{
    if(orients[inseq1]){
      rev = true;
      std::cerr << "Unsupported orientation '-' for seq1" << std::endl;
      exit(1);
    }
    else{
      if(orients[inseq2]){
	rev = true;
      }
    }
  } 
  if(rev){
#ifdef DEBUGGING
    std::cout << "S1: " << inseq1 << " " << length(value(strSet, inseq1)) << " " << coords[inseq1] << ":" << startMatch1+coords[inseq1] 
	      << " o:" << orients[inseq1] << " len:" << len << std::endl;
    std::cout << "S2: " << inseq2 << " " << length(value(strSet, inseq2)) << " startMatch2:" << startMatch2 << " " << coords[inseq2] << ":"  <<  length(value(strSet, inseq2))-coords[inseq2]+len-startMatch2
      
	      << " o:" << orients[inseq2] << " len:" << len << std::endl;
#endif
    assert(startMatch1+coords[inseq1]<length(value(strSet, inseq1)));
    assert(startMatch1+coords[inseq1]+len<=length(value(strSet, inseq1)));
    //assert(length(value(strSet, inseq2)) - startMatch2+coords[inseq2] >= 0);
    //assert(length(value(strSet, inseq2)) - (startMatch2+coords[inseq2]+len) >= 0);
    //assert(length(value(strSet, inseq2)) - (startMatch2+coords[inseq2]+len)<=length(value(strSet, inseq2)));
    
    //Original
    //__appendNewMatch(matches, scores, inseq1, inseq2, 
    //startMatch1+coords[inseq1], (length(value(strSet, inseq2)) - (startMatch2+coords[inseq2]+len)), len,rev);
    //I think this is correct but progressive align doesn't work right with this one
    //__appendNewMatch(matches, scores, inseq1, inseq2, 
    //startMatch1+coords[inseq1], coords[inseq2]-len-startMatch2, len,rev);
    //__appendNewMatch(matches, scores, inseq1, 
    //estartMatch1+coords[inseq1], length(value(strSet, inseq2))-(coords[inseq2]-startMatch2), len,rev);
#ifdef REVCOORD
    __appendNewMatch(matches, scores, inseq1, inseq2, 
		     startMatch1+coords[inseq1], coords[inseq2]+startMatch2, len,rev);
#else
    assert((int)(startMatch1+coords[inseq1])>=0);
    assert((int)(startMatch2+coords[inseq2])>=0);
    assert(startMatch1+coords[inseq1]+len<=length(value(strSet, inseq1)));
    assert(length(value(strSet, inseq2))-coords[inseq2]-startMatch2<=length(value(strSet, inseq2)));
    __appendNewMatch(matches, scores, inseq1, inseq2, 
		     startMatch1+coords[inseq1], length(value(strSet, inseq2))-coords[inseq2]-startMatch2-len, len,rev);
#endif
  }
  else{
#ifdef DEBUGGING
    std::cout << "S1: " << inseq1 << " " << length(value(strSet, inseq1)) << " " << coords[inseq1] << ":" << startMatch1+coords[inseq1] 
	      << " o:" << orients[inseq1] << " len:" << len << std::endl;
    std::cout << "S2: " << inseq2 << " " << length(value(strSet, inseq2)) << " " << coords[inseq2] << ":"  <<  startMatch2+coords[inseq2]
      
	      << " o:" << orients[inseq2] << " len:" << len << std::endl;
#endif
    assert(inseq1<length(coords));
    assert(inseq2<length(coords));
    assert(inseq1<length(orients));
    assert(inseq2<length(orients));
    assert(startMatch1+coords[inseq1]<length(value(strSet, inseq1)));
    assert(startMatch2+coords[inseq2]<length(value(strSet, inseq2)));
    //assert(startMatch1+coords[inseq1]>=0);
    //assert(startMatch2+coords[inseq2]>=0);
    assert(startMatch1+coords[inseq1]+len<=length(value(strSet, inseq1)));
    assert(startMatch2+coords[inseq2]+len<=length(value(strSet, inseq2)));
    __appendNewMatch(matches, scores, inseq1, inseq2, 
		     startMatch1+coords[inseq1], startMatch2+coords[inseq2], len,rev);
  }
}


//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TFragment, typename TSpec, typename TScoreValue, typename TSpec2, typename TNames>
void 
read(TFile & file,
	 String<TFragment, TSpec>& matches,
	 String<TScoreValue, TSpec2>& scores,
	 TNames const& origNames,
	 FastaAlign) 
{
	SEQAN_CHECKPOINT
	typedef typename Size<TNames>::Type TSize;
	typedef typename Value<TFile>::Type TValue;
	typedef typename Value<TNames>::Type TName;
	
	TValue c;
	if (_streamEOF(file)) return;
	else c = _streamGet(file);

	// Read sequences
	String<TValue> mat;
	TNames names;
	TName nextSeq;
	while(!_streamEOF(file)) {
		_parse_skipWhitespace(file, c);
		if (_streamEOF(file)) break;
		if (c == '>') {
			c = _streamGet(file);
			clear(nextSeq);
			_parse_readIdentifier(file, nextSeq, c);
			appendValue(names, nextSeq);
			_parse_skipLine(file, c);
		} else if ((c == '\n') || (c == '\r')) {
			c = _streamGet(file);
		} else {
			appendValue(mat, c);
			c = _streamGet(file);
		}
	}
	// Reorder rows according to names order
	String<TValue> finalMat = mat;
	TSize nseq = length(names);
	TSize len = length(mat) / nseq;
	for(TSize i = 0; i<nseq; ++i) {
		if (value(names, i) == value(origNames, i)) continue;
		else {
			for(TSize j = 0; j<length(origNames); ++j) {
				if (value(names, i) != value(origNames, j)) continue;
				// Copy the whole row
				infix(finalMat, j * len, j * len + len) = infix(mat, i * len, i*len + len);
				break;
			}
		}
	}
	clear(mat);

	// Collect the segment matches
	_collectSegmentMatches(finalMat, matches, scores, nseq); 
}


 template<typename TFile, typename TFragment, typename TSpec, typename TScoreValue, typename TSpec2, typename TStringSet, typename TNames>
  void 
  read(TFile & file,
       String<TFragment, TSpec>& matches,
       String<TScoreValue, TSpec2>& scores,
       TStringSet const& strSet,
       TNames const& origNames,
       MultiFastaAlign) 
{
  SEQAN_CHECKPOINT
  typedef typename Size<TNames>::Type TSize;
  typedef typename Value<TFile>::Type TValue;
  typedef typename Value<TNames>::Type TName;
  
  TValue c;
  TValue gapChar = gapValue<TValue>();
  if (_streamEOF(file)) return;
  else c = _streamGet(file);
  
  // Map the names to slots
  typedef std::map<TName, TSize> TNameToPosition;
  TNameToPosition namePosMap;
  
  //TSize seq;
  TName currSeq;

  // Read sequences
  String<TValue> mat;
  String<TSize> coords;
  String<bool> orients;
  TNames names;

  TSize beg,flen,seqlen;
  bool rev;

  std::map<TName,int> name2pos;
  assert(length(origNames)>0);
  for(unsigned int j = 0; j<length(origNames); ++j) {
    assert(value(origNames,j)==origNames[j]);
  }
  for(TSize j = 0; j<length(origNames); ++j) {  
    name2pos[value(origNames,j)] = j;
  }
	
  //typedef std::map<TName, std::pair<TSize,TSize> > TNameToCoords;
  //typedef std::map<TName, bool > TNameToOrient;
  //TNameToCoords coordsMap;
  //TNameToOrient orientMap;

  beg=0;
  flen=0;
  seqlen=0;
  rev=false;
#ifdef SEQAN_PROFILE
	  std::cerr << "Parsing pairwise alignments ";
#endif
  unsigned int numAligns=0;

  unsigned int ungappedlen=0;
  unsigned int reported_ungappedlen=0;

  while(!_streamEOF(file)) {
    _parse_skipWhitespace(file, c);
    if (_streamEOF(file)) break;
    if (c == '>') {
      c = _streamGet(file);
      clear(currSeq);
      _parse_readIdentifier(file, currSeq, c);
      appendValue(names, currSeq);
      _parse_skipWhitespace(file, c);
      //Start coord
      beg = _parse_readNumber(file, c);
      _parse_skipWhitespace(file, c);
      //Len coord
      flen = _parse_readNumber(file, c);
      _parse_skipWhitespace(file, c);
      //Read orientation
      if(c == '+'){
	rev = false;
      }
      else{
	if(c == '-'){
	  rev = true;
	}
	else{
	  std::cerr << "Unknown orient string:" << c << std::endl;
	}
      }
      //gobble orient
      _parse_readNumber(file, c);
      _parse_skipWhitespace(file, c);
      seqlen = _parse_readNumber(file, c);
      assert(seqlen > 0);
      if(rev){
#ifdef REVCOORD
	//Coords relative to the the reverse strand 
	appendValue(coords,beg);
#else
	//Coords relative to the the reverse strand 
	//std::cout << "beg:"<< beg << " flen:" << flen << " seqlen:" << seqlen << std::endl;
	//assert(seqlen-(beg+flen)>=0);
	//assert(seqlen-(beg+flen)<seqlen);
	//assert(seqlen-(beg+flen)+flen<=seqlen);

	assert((int)beg>=0);
	assert(beg<seqlen);
	assert(beg+flen<=seqlen);
	//Coords relative to the leading strand
	//appendValue(coords,seqlen-(beg+flen));
	appendValue(coords,beg);
#endif
      }
      else{
	appendValue(coords,beg);
      }
      appendValue(orients,rev);
      reported_ungappedlen+=flen;
      _parse_skipLine(file, c);
    } else if ((c == '\n') || (c == '\r')) {
      c = _streamGet(file);
    } else if (c == '='){
      //end of alignment
      // Reorder rows according to names order

      String<TSize> finalcoords = coords;
      String<bool> finalorients = orients;
      resize(finalcoords,length(origNames));
      resize(finalorients,length(origNames));

      TSize nseq = length(names);
      if(nseq>2){
	std::cerr << "Only pairwise matches currently supported nseq:" << nseq << std::endl;
	exit(1);
      }
      assert(nseq==2);
      TSize seq1=0,seq2=0;
      bool isSeq1=false;
      bool isSeq2=false;
      bool isSeq1rev=false,isSeq2rev=false;
      if(nseq<=0){
	std::cerr << "No sequences" << std::endl;
	assert(false);
      }
      if(length(mat) % nseq != 0){
	std::cerr << "Bad length of input matrix:" << length(mat) 
		  << ". Not divisible by number of sequences:" << nseq 
		  << std::endl;
	assert(false);
      }

      TSize len = length(mat) / nseq;
      String<TValue> finalMat = mat;
      resize(finalMat,length(origNames)*len);
      if(name2pos.find(value(names,0))!=name2pos.end()){
	seq1 = name2pos[value(names,0)];
	isSeq1=true;
	isSeq2rev = orients[0];
	assert(coords[0]<length(value(strSet, seq1)));
	infix(finalMat, seq1 * len, seq1 * len + len) = infix(mat, 0,len);
      }
      else{
	isSeq1=false;
	//assert(false);
      }

      //      infix(finalMat, j * len, j * len + len) = infix(mat, i * len, i * len + len);


      if(name2pos.find(value(names,1))!=name2pos.end()){
	seq2 = name2pos[value(names,1)];
	isSeq2=true;
	isSeq2rev = orients[1];
	assert(coords[1]<length(value(strSet, seq2)));
	infix(finalMat, seq2 * len, seq2 * len + len) = infix(mat, len, len + len);
      }
      else{
	isSeq2=false;
	//assert(false);
      }

      
      /*
      for(TSize i = 0; i<nseq; ++i) {
	for(TSize j = 0; j<length(origNames); ++j) {
	  if (value(names, i) == value(origNames, j)){
	    //std::cout << "Orig " << value(origNames,j) << " " << length(value(strSet, j))
	    //<< " localindex " << i 
	    //<< " finalindex " << j << " finalmatpos:" << j*len << "-" << j*len+len 
	    //<< " localmatrix len:" << length(mat) << std::endl;
	    if(!isSeq1){
	      seq1 = j;
	      isSeq1 = true;
	      isSeq1rev = orients[i];
	      assert(coords[i]<length(value(strSet, seq1)));
	    }
	    else{
	      if(!isSeq2){
		seq2 = j;
		isSeq2 = true;
		isSeq2rev = orients[i];
		assert(coords[i]<length(value(strSet, seq2)));
	      }
	      else{
		std::cerr << "Only pairwise matches currently supported" << std::endl;
		assert(false);
	      }
	    }
	    assert(j<length(finalcoords));
	    assert(j<length(finalorients));
	    finalcoords[j] = coords[i];
	    finalorients[j] = orients[i];
	    if(i * len+len>length(mat)){
	      std::cerr << "Bad value:" << i * len+len << " for matrix of length:" << length(mat) << std::endl;
	      assert(false);
	    }
	    if(j * len+len>length(finalMat)){
	      std::cerr << "Bad value:" << j * len+len << " for matrix of length:" << length(finalMat) << std::endl;
	      assert(false);
	    }
	    infix(finalMat, j * len, j * len + len) = infix(mat, i * len, i * len + len);
	    break;
	  }
	}
      }
      */

      if(isSeq1 && isSeq2){
	assert(isSeq1 && isSeq2);
	assert(len+len<=length(mat));
	assert(len+len<=length(finalMat));
	
	if(reported_ungappedlen != ungappedlen){
	  std::cerr << "Bad length of matrix per sequence:" << ungappedlen << " expecting " << reported_ungappedlen << std::endl;
	  assert(false);
	}
	assert(seq1<length(finalcoords));
	assert(seq2<length(finalcoords));
	
	finalcoords[seq1] = coords[0];
	finalorients[seq1] = orients[0];
	finalcoords[seq2] = coords[1];
	finalorients[seq2] = orients[1];
	//std::cerr << "Input seqs " << origNames[seq1] << " " << origNames[seq2] << std::endl;
	//if(length(names)!=length(finalcoords) || length(names)!=length(finalorients)){
	assert(idToPosition(strSet, seq1)==seq1);
	assert(idToPosition(strSet, seq2)==seq2);
	if(length(origNames)!=length(finalcoords) || length(origNames)!=length(finalorients)){
	  std::cerr << "Bad length of coords,orients" << std::endl;
	  assert(false);
	}
	assert(!isSeq1rev);
	unsigned int ali_counter = length(matches);
	if(!isSeq1rev && !isSeq2rev){
	  // Collect the segment matches
	  _collectSegmentMatchesOffset(finalMat, matches, scores, length(origNames), 
				       finalcoords, finalorients, 
				       strSet,
				       seq1, seq2);
	}
	else{
	  //finalcoords here for finalcoords[seq2] should be relative to the end of the match
	  _collectSegmentMatchesOffset(finalMat, matches, scores, length(origNames), 
				       finalcoords, finalorients, 
				       strSet,
				       seq1, seq2);
	  //Perform some simple checks for correctness of seq ids
	  typedef String<Fragment<> > TFragmentString;
	  typedef String<Fragment<> > TAlignmentString;
	  typedef typename Iterator<TAlignmentString,Standard>::Type TAliIterator;
	  TAliIterator ali_it = begin(matches,Standard())+ali_counter;
	  TAliIterator ali_end = end(matches,Standard());
	  unsigned int seq_i_id,begin_,end_;
	  while(ali_it != ali_end)
	    {
	      Fragment<> segment = *ali_it;
	      seq_i_id = sequenceId(segment,0);
	      begin_ = fragmentBegin(segment,seq_i_id);
	      end_ = begin_ + fragmentLength(segment,seq_i_id);

	      assert(seq_i_id==seq1);
	      assert(seq_i_id==idToPosition(strSet,seq_i_id));

	      seq_i_id = sequenceId(segment,1);
	      begin_ = fragmentBegin(segment,seq_i_id);
	      end_ = begin_ + fragmentLength(segment,seq_i_id);

	      assert(seq_i_id==seq2);
	      assert(seq_i_id==idToPosition(strSet,seq_i_id));

	      ++ali_counter;
	      ++ali_it;
	    }
	}
	if(++numAligns%1000==0){
#ifdef SEQAN_PROFILE
	  std::cerr << numAligns << " ";
	  std::cerr.flush();
#endif
	}
      }
      else{
	std::cerr << "Unable to find pairwise matches "  
		  << isSeq1 << ":" << seq1 << " name:" << value(names,0) 
		  << " -- " 
		  << isSeq2 << ":" << seq2 << " name:" << value(names,1) 
		  << std::endl;
	//assert(false);
      }

      clear(mat);
      clear(coords);
      clear(orients);
      clear(names);
      ungappedlen=0;
      reported_ungappedlen=0;

      _parse_skipLine(file, c);
    }
    else{
      if(c != gapChar){
	ungappedlen++;
      }
      appendValue(mat, c);
      c = _streamGet(file);
    }
  }
#ifdef DEBUGGING
  std::cout << std::endl;
#endif
  //assert(length(matches)>0); zero matches is also ok
  return;
}

/////////////////////////////////////////////////////////////////////////////
// BLAST Library Reading / Writing
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TSizeSpec, typename TSpec1, typename TSpec2, typename TSize>
inline void 
__includeFragment(String<Fragment<TSizeSpec, ExactReversableFragment<TSpec1> >, TSpec2>& matches, 
				  TSize seq1Id, 
				  TSize beg1, 
				  TSize seq2Id, 
				  TSize beg2, 
				  TSize len, 
				  bool reversed)
{
	SEQAN_CHECKPOINT
	  typedef Fragment<TSizeSpec, ExactReversableFragment<TSpec1> > TFragment;
	  
	//appendValue(matches, TFragment(seq1Id, beg1, seq2Id, beg2, len, len, reversed));
	appendValue(matches, TFragment(seq1Id, beg1, seq2Id, beg2, len, reversed));
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSizeSpec, typename TSpec1, typename TSpec2, typename TSize>
inline void 
  __includeFragment(String<Fragment<TSizeSpec, ExactFragment<TSpec1> >, TSpec2>& matches, 
				  TSize seq1Id, 
				  TSize beg1, 
				  TSize seq2Id, 
				  TSize beg2, 
				  TSize len, 
				  bool reversed)
{
	SEQAN_CHECKPOINT
	  typedef Fragment<TSizeSpec, ExactFragment<TSpec1> > TFragment;
	//if (!reversed) appendValue(matches, TFragment(seq1Id, beg1, seq2Id, beg2, len,len));
	if (!reversed) appendValue(matches, TFragment(seq1Id, beg1, seq2Id, beg2, len));
}

//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TFragment, typename TSpec1, typename TScoreValue, typename TSpec2, typename TNames>
inline void 
read(TFile & file,
	 String<TFragment, TSpec1>& matches,
	 String<TScoreValue, TSpec2>& scores,
	 TNames& names,
	 BlastLib) 
{
	SEQAN_CHECKPOINT
	typedef typename Size<TNames>::Type TSize;
	typedef typename Value<TFile>::Type TValue;
	typedef typename Value<TNames>::Type TName;
	
	// Map the names to slots
	typedef std::map<TName, TSize> TNameToPosition;
	TNameToPosition namePosMap;
	for(TSize i = 0;i<length(names);++i) namePosMap.insert(std::make_pair(names[i], i));
	
	// Read the Blast file
	TValue c;
	if (_streamEOF(file)) return;
	else c = _streamGet(file);

	TName seq1;
	TName seq2;
	while (!_streamEOF(file)) {
		clear(seq1);
		clear(seq2);
		_parse_skipWhitespace(file, c);
		seq1 = _parse_readIdentifier(file, c);
		_parse_skipWhitespace(file, c);
		seq2 = _parse_readIdentifier(file, c);
		if (seq1 == seq2) {
			_parse_skipLine(file, c);
			continue;
		}
		TSize seq1Id = namePosMap[seq1];
		TSize seq2Id = namePosMap[seq2];
		_parse_skipWhitespace(file, c);
		_parse_readDouble(file, c);
		_parse_skipWhitespace(file, c);
		_parse_readNumber(file, c);
		_parse_skipWhitespace(file, c);
		_parse_readNumber(file, c);
		_parse_skipWhitespace(file, c);
		_parse_readNumber(file, c);
		_parse_skipWhitespace(file, c);
		TSize beg1 = _parse_readNumber(file, c);
		_parse_skipWhitespace(file, c);
		TSize end1 = _parse_readNumber(file, c);
		TSize len = end1 - beg1 + 1;
		_parse_skipWhitespace(file, c);
		TSize beg2 = _parse_readNumber(file, c);
		_parse_skipWhitespace(file, c);
		TSize end2 = _parse_readNumber(file, c);
		_parse_skipWhitespace(file, c);
		_parse_readIdentifier(file, c);
		_parse_skipWhitespace(file, c);
		TScoreValue rawScore = (TScoreValue) _parse_readDouble(file, c);

		bool reversed = false;
		if (beg1 > end1) { TSize tmp = beg1; beg1 = end1; end1 = tmp; reversed = !reversed; }
		if (beg2 > end2) { TSize tmp = beg2; beg2 = end2; end2 = tmp; reversed = !reversed; }
		//// Debug code
		//std::cout << seq1Id << ',' << beg1 << ',' << seq2Id << ',' << beg2 << ',' << len << std::endl;
		//std::cout << infix(strSet[seq1Id], beg1, beg1+len) << std::endl;
		//std::cout << infix(strSet[seq2Id], beg2, beg2+len) << std::endl;

		__includeFragment(matches, seq1Id, --beg1, seq2Id, --beg2, len, reversed);
		appendValue(scores, rawScore);
		_parse_skipLine(file, c);
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TStringSet, typename TCargo, typename TSpec, typename TNames, typename TEdgeMap>
void write(TFile & file, 
		   Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
		   TNames& names,
		   TEdgeMap& edgeMap,
		   BlastLib) 
{
	SEQAN_CHECKPOINT
	typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename Value<TStringSet>::Type TString;
	typedef typename Size<TStringSet>::Type TSize;
	
	TStringSet& str = stringSet(g);
	
	typedef typename Iterator<TGraph, EdgeIterator>::Type TIter;
	TIter it(g);
	for(;!atEnd(it);++it) {
		TVertexDescriptor sV = sourceVertex(it);
		TVertexDescriptor tV = targetVertex(it);
		if (sequenceId(g,sV) < sequenceId(g,tV)) {
			TVertexDescriptor tmp = sV;
			sV = tV;
			tV = tmp;
		}
		TSize fragLen = fragmentLength(g,sV);
		TSize fragPos1 = fragmentBegin(g,sV);
		TSize fragPos2 = fragmentBegin(g,tV);
		TSize seq1 = idToPosition(str, sequenceId(g,sV));
		TSize seq2 = idToPosition(str, sequenceId(g,tV));
		TCargo my_carg =  getCargo(*it);
		_streamWrite(file, names[seq1]);
		_streamPut(file, '\t');	
		_streamWrite(file, names[seq2]);
		_streamPut(file, '\t');	
		_streamPutInt(file, 0);
		_streamPut(file, '\t');	
		_streamPutInt(file, fragLen);
		_streamPut(file, '\t');	
		_streamPutInt(file, 0);
		_streamPut(file, '\t');	
		_streamPutInt(file, 0);
		_streamPut(file, '\t');	
		_streamPutInt(file, fragPos1+1);
		_streamPut(file, '\t');	
		_streamPutInt(file, fragPos1 + fragLen);
		_streamPut(file, '\t');	
		if (!property(edgeMap, *it)) {
			_streamPutInt(file, fragPos2+1);
			_streamPut(file, '\t');	
			_streamPutInt(file, fragPos2 + fragLen);
			_streamPut(file, '\t');	
		} else {
			_streamPutInt(file, fragPos2 + fragLen);
			_streamPut(file, '\t');		
			_streamPutInt(file, fragPos2+1);
			_streamPut(file, '\t');	
		}
		_streamPutInt(file, 0);
		_streamPut(file, '\t');	
		_streamPutInt(file, my_carg);
		_streamPut(file, '\n');
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TStringSet, typename TCargo, typename TSpec, typename TNames>
void write(TFile & file, 
		   Graph<Alignment<TStringSet, TCargo, TSpec> > const& g,
		   TNames& names,
		   BlastLib) 
{
	SEQAN_CHECKPOINT
	String<bool> edgeMap;
	fill(edgeMap, getIdUpperBound(_getEdgeIdManager(g)), false);
	write(file, g, names, edgeMap, BlastLib());
}



/////////////////////////////////////////////////////////////////////////////
// MUMMER Format Reading
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

template<typename TPos, typename TSpec2, typename TSpec1, typename TScores, typename TId, typename TSize>
inline void 
  __appendNewMatch(String<Fragment<TPos, ExactReversableFragment<TSpec2> >, TSpec1>& matches,
		   TScores& scores,
				 TId seq1Id,
				 TId seq2Id,
				 TSize beg1,
				 TSize beg2,
				 TSize len,
				 bool reversed) 
{
	typedef Fragment<TPos, ExactReversableFragment<TSpec2> > TFragment;
	appendValue(matches, TFragment(seq1Id, beg1, seq2Id, beg2, len, reversed));
	appendValue(scores, len);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TPos, typename TSpec2, typename TSpec1, typename TScores, typename TId, typename TSize>
inline void 
  __appendNewMatch(String<Fragment<TPos, ExactFragment<TSpec2> >, TSpec1>& matches,
				 TScores& scores,
				 TId seq1Id,
				 TId seq2Id,
				 TSize beg1,
				 TSize beg2,
				 TSize len,
				 bool reversed) 
{
	typedef Fragment<TPos, Fragment<TSpec2> > TFragment;
	appendValue(matches, TFragment(seq1Id, beg1, seq2Id, beg2, len,reversed));
	appendValue(scores, len);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TFragment, typename TSpec1, typename TScoreValue, typename TSpec2, typename TStringSet, typename TNames>
inline void 
read(TFile & file,
	 String<TFragment, TSpec1>& matches,
	 String<TScoreValue, TSpec2>& scores,
	 TStringSet const& strSet,
	 TNames const& names,
	 MummerLib) 
{
	SEQAN_CHECKPOINT
	typedef typename Size<TNames>::Type TSize;
	typedef typename Value<TFile>::Type TValue;
	typedef typename Value<TNames>::Type TName;
	
	// Map the names to slots
	typedef std::map<TName, TSize> TNameToPosition;
	TNameToPosition namePosMap;
	for(TSize i = 0;i<length(names);++i) namePosMap.insert(std::make_pair(value(names, i), i));
	
	// Read the Mummer file
	TValue c;
	if (_streamEOF(file)) return;
	else c = _streamGet(file);
	TName seq1;
	TName seq2;
	TSize seq1Id = 0;
	TSize seq2Id = 0;
	bool reversed = false;
	while (!_streamEOF(file)) {
		if (c == '>') {
			c = _streamGet(file);
			_parse_skipWhitespace(file, c);
			seq1 = _parse_readIdentifier(file, c);
			 seq1Id = namePosMap[seq1];
			 _parse_skipWhitespace(file, c);
			 if (c == 'R') {
				 reversed = true;
				 _parse_skipLine(file, c);
			 } else reversed = false;
		} else {
			_parse_skipWhitespace(file, c);
			if (_streamEOF(file)) {
				break;
			}
			seq2 = _parse_readIdentifier(file, c);
			
			
			seq2Id = namePosMap[seq2];
			_parse_skipWhitespace(file, c);
			TSize beg2 = _parse_readNumber(file, c);
			_parse_skipWhitespace(file, c);
			TSize beg1 = _parse_readNumber(file, c);
			_parse_skipWhitespace(file, c);
			TSize len = _parse_readNumber(file, c);
			_parse_skipLine(file, c);
			if (seq1Id == seq2Id) continue;
			if (!reversed) __appendNewMatch(matches, scores, seq1Id, seq2Id, --beg1, --beg2, len, reversed);
			else __appendNewMatch(matches, scores, seq1Id, seq2Id, (length(value(strSet, seq1Id)) - (--beg1 + len)), --beg2, len, reversed);
		}
	}
}

/////////////////////////////////////////////////////////////////////////////
// Newick Format
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TCargo, typename TSpec, typename TNames>
void 
read(TFile & file,
	 Graph<Tree<TCargo, TSpec> >& guideTree,
	 TNames& names,
	 NewickFormat) 
{
	SEQAN_CHECKPOINT
	typedef Graph<Tree<TCargo, TSpec> > TGuideTree;
	typedef typename VertexDescriptor<TGuideTree>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGuideTree>::Type TEdgeDescriptor;
	typedef typename Size<TGuideTree>::Type TSize;
	typedef typename Id<TGuideTree>::Type TId;
	typedef typename Position<TFile>::Type TPosition;
	typedef typename Value<TFile>::Type TValue;
	typedef typename Value<TNames>::Type TName;
	TVertexDescriptor nilVertex = getNil<TVertexDescriptor>();


	if (length(names) < 3) {
		TVertexDescriptor v1 = addVertex(guideTree);
		TVertexDescriptor v2 = addVertex(guideTree);
		TVertexDescriptor internalVertex = addVertex(guideTree);
		addEdge(guideTree, internalVertex, v1, (TCargo) 1 * SEQAN_DISTANCE_UNITY);
		addEdge(guideTree, internalVertex, v2, (TCargo) 1 * SEQAN_DISTANCE_UNITY);
		assignRoot(guideTree, internalVertex);
		return;
	}


	typedef std::map<TName, TId> TNameToId;
	TNameToId nameToId;
	for(TId i=0; i<length(names);++i) {
		addVertex(guideTree);	// Create the sequence vertices
		nameToId.insert(std::make_pair(names[i], i));
	}

	TValue c;
	if (_streamEOF(file)) return;
	else c = _streamGet(file);

	TVertexDescriptor lastVertex = nilVertex;
	TVertexDescriptor lastChild = nilVertex;
	while (!_streamEOF(file)) {
		if (c=='(') {
			if (lastVertex == nilVertex) {
				lastVertex = addVertex(guideTree);
				assignRoot(guideTree, lastVertex);
			} else {
				TVertexDescriptor ch = addChild(guideTree, lastVertex);
				lastVertex = ch;
			}
			c = _streamGet(file);
			_parse_skipWhitespace(file, c);
		} else if (c==')') {
			if (!isRoot(guideTree, lastVertex)) {
				lastChild = lastVertex;
				lastVertex = parentVertex(guideTree, lastVertex);
			} else {
				lastChild = lastVertex;
				lastVertex = nilVertex;
			}
			c = _streamGet(file);
			_parse_skipWhitespace(file, c);
		} else if (c==',') {
			c = _streamGet(file);
			_parse_skipWhitespace(file, c);
		} else if (c==':') {
			c = _streamGet(file);
			cargo(findEdge(guideTree, lastVertex, lastChild)) = (TCargo) (_parse_readDouble(file,c) * SEQAN_DISTANCE_UNITY);
		} else if (c==';') {
			c = _streamGet(file);
			_parse_skipWhitespace(file, c);
		} else {
			TName tmp = _parse_readIdentifier(file, c);
			//std::cout << tmp << std::endl;
			if (lastVertex == nilVertex) {
				// Tree is rooted at a leaf
				// Create artificial root node
				lastVertex = length(names);
				assignRoot(guideTree, addVertex(guideTree));
				addEdge(guideTree, getRoot(guideTree), lastVertex);
				addEdge(guideTree, getRoot(guideTree), nameToId[tmp]);
			} else {
				addEdge(guideTree, lastVertex, nameToId[tmp]);
			}
			lastChild = nameToId[tmp];
		}
	}
	
	// Root the tree if necessary 
	if (outDegree(guideTree, lastChild) > 2) {
		TVertexDescriptor myRoot = addVertex(guideTree);
		assignRoot(guideTree, myRoot);
		typedef typename Iterator<TGuideTree, OutEdgeIterator>::Type TOutEdgeIterator;
		TOutEdgeIterator it(guideTree, lastChild);
		goNext(it); goNext(it);
		TVertexDescriptor tV = targetVertex(it);
		TCargo c = cargo(*it);
		removeEdge(guideTree, lastChild, tV);
		addEdge(guideTree, myRoot, tV, (TCargo) (c / 2));
		addEdge(guideTree, myRoot, lastChild, (TCargo) (c / 2));
	}

	//std::fstream strm1; // Alignment graph as dot
	//strm1.open("tree23.dot", std::ios_base::out | std::ios_base::trunc);
	//write(strm1,guideTree,DotDrawing());
	//strm1.close();

	//std::cout << guideTree << std::endl;
}


//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec, typename TNames, typename TNewickString, typename TVertexDescriptor>
void 
_buildNewickString(Graph<Tree<TCargo, TSpec> >& guideTree,
				   TNames& names,
				   TNewickString& str,
				   TVertexDescriptor v,
				   bool collapseRoot) 
{
	typedef Graph<Tree<TCargo, TSpec> > TGuideTree;
	typedef typename EdgeDescriptor<TGuideTree>::Type TEdgeDescriptor;
	typedef typename Iterator<TGuideTree, AdjacencyIterator>::Type TAdjIter;
	if (isLeaf(guideTree, v)) {
		append(str, names[v], Generous());
	} else {
		if ((collapseRoot) && (isRoot(guideTree, v))) {
			TAdjIter adjIterRoot(guideTree, v);
			TVertexDescriptor v1 = *adjIterRoot;
			goNext(adjIterRoot);
			TVertexDescriptor v2 = *adjIterRoot;
			if (isLeaf(guideTree, v2)) {
				TVertexDescriptor tmp = v1;
				v1 = v2;
				v2 = tmp;
			}
			String<char> subStr;
			appendValue(subStr, '(', Generous());
			_buildNewickString(guideTree, names, subStr, v1, collapseRoot);
			TEdgeDescriptor e1 = findEdge(guideTree, v, v1);
			TEdgeDescriptor e2 = findEdge(guideTree, v, v2);
			appendValue(subStr, ':', Generous());
			::std::ostringstream weight;
			weight << (cargo(e1) + cargo(e2));
			append(subStr, weight.str().c_str(), Generous());
			TAdjIter adjIter(guideTree, v2);
			for(;!atEnd(adjIter); goNext(adjIter)) {
				appendValue(subStr, ',', Generous());
				_buildNewickString(guideTree, names, subStr, *adjIter, collapseRoot);
				TEdgeDescriptor e = findEdge(guideTree, v2, *adjIter);
				appendValue(subStr, ':', Generous());
				::std::ostringstream weight;
				weight << cargo(e);
				append(subStr, weight.str().c_str(), Generous());
			}
			appendValue(subStr, ')', Generous());
			append(str, subStr, Generous());
		} else {
			TAdjIter adjIter(guideTree, v);
			String<char> subStr;
			appendValue(subStr, '(', Generous());
			_buildNewickString(guideTree, names, subStr, *adjIter, collapseRoot);
			TEdgeDescriptor e = findEdge(guideTree, v, *adjIter);
			appendValue(subStr, ':', Generous());
			::std::ostringstream weight;
			weight << cargo(e);
			append(subStr, weight.str().c_str(), Generous());
			goNext(adjIter);
			for(;!atEnd(adjIter); goNext(adjIter)) {
				appendValue(subStr, ',', Generous());
				_buildNewickString(guideTree, names, subStr, *adjIter, collapseRoot);
				TEdgeDescriptor e = findEdge(guideTree, v, *adjIter);
				appendValue(subStr, ':', Generous());
				::std::ostringstream weight;
				weight << cargo(e);
				append(subStr, weight.str().c_str(), Generous());
			}
			appendValue(subStr, ')', Generous());
			append(str, subStr, Generous());
		}
	}
}


//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TCargo, typename TSpec, typename TNames>
void 
write(TFile & file,
	  Graph<Tree<TCargo, TSpec> >& guideTree,
	  TNames& names,
	  bool collapseRoot,
	  NewickFormat) 
{
	typedef Graph<Tree<TCargo, TSpec> > TGuideTree;
	typedef typename Size<TGuideTree>::Type TSize;
	typedef typename VertexDescriptor<TGuideTree>::Type TVertexDescriptor;
	
	String<char> myNewickString;
	_buildNewickString(guideTree, names, myNewickString, getRoot(guideTree), collapseRoot);
	_streamWrite(file, myNewickString);
	_streamPut(file, ';');	
	_streamPut(file, '\n');	
}

//////////////////////////////////////////////////////////////////////////////

template<typename TFile, typename TCargo, typename TSpec, typename TNames>
void 
write(TFile & file,
	  Graph<Tree<TCargo, TSpec> >& guideTree,
	  TNames& names,
	  NewickFormat) 
{
	write(file,guideTree, names, false, NewickFormat());
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...

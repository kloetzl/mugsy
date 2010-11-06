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
  $Id: graph_align_tcoffee_msa.h 1911 2008-05-02 09:28:04Z rausch@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_GRAPH_ALIGN_TCOFFEE_MSA_H
#define SEQAN_HEADER_GRAPH_ALIGN_TCOFFEE_MSA_H

namespace SEQAN_NAMESPACE_MAIN
{

template<typename TAlphabet, typename TScore>
struct MsaOptions {
public:
	// Rescore segment matches after refinement
	bool rescore;

	// Output format
	// 0: Fasta
	// 1: Msf
	unsigned int outputFormat;

	// Scoring object
	TScore sc;

	// All methods to compute a guide tree
	// 0: Neighbor-joining
	// 1: UPGMA single linkage
	// 2: UPGMA complete linkage
	// 3: UPGMA average linkage
	// 4: UPGMA weighted average linkage
	unsigned int build;

	// All methods to compute segment matches
	// 0: global alignment
	// 1: local alignments
	// 2: overlap alignments
	// 3: longest common subsequence
	String<unsigned int> method;

	// Various input and output file names
	String<std::string> alnfiles;		// External alignment files
	String<std::string> libfiles;		// T-Coffee library files
	String<std::string> blastfiles;		// Blast match files
	String<std::string> mummerfiles;	// MUMmer files
	std::string outfile;				// Output file name
	std::string seqfile;				// Sequence file name
	std::string infile;					// Alignment file for alignment evaluation
	std::string treefile;				// Guide tree file

  //MUGSY options
  std::string distance;				// 
  std::string minlength;				// 
  std::string refine;				// 
  std::string unique;				// 
  std::string blockfile;
  std::string segmentation;
  std::string duplications;                   
  std::string allownestedlcbs;

  unsigned int anchorwin;

  //TESTING, disabled options
  //unsigned int partition; 
  //afloat possharedcutoff;
  //unsigned int posscorewindow;
  //unsigned int filter;

	// Initialization
	MsaOptions() : rescore(true), outputFormat(0), build(0) {}
};

//////////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TScore>
inline void
evaluateAlignment(MsaOptions<TAlphabet, TScore> const& msaOpt) {
	typedef typename Value<TScore>::Type TScoreValue;
	typedef String<TAlphabet> TSequence;
	typedef typename Size<TSequence>::Type TSize;
	StringSet<TSequence, Owner<> > origStrSet;
	StringSet<String<char> > names;

	// Read the sequences
	std::fstream strm;
	strm.open(msaOpt.infile.c_str(), ::std::ios_base::in | ::std::ios_base::binary);
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
	strm_lib.open(msaOpt.infile.c_str(), ::std::ios_base::in | ::std::ios_base::binary);
	read(strm_lib,matches, scores, names, FastaAlign());	
	strm_lib.close();

	// Build the alignment graph
	//SVA make type int to support negative edge weights
	typedef Graph<Alignment<TDepSequenceSet, int> > TGraph;
	//typedef Graph<Alignment<TDepSequenceSet, TSize> > TGraph;
	TGraph g(strSet);
	buildAlignmentGraph(matches, g, FrequencyCounting() );

	// Print the scoring information
	TScoreValue gop = msaOpt.sc.data_gap_open;
	TScoreValue gex = msaOpt.sc.data_gap_extend;
	std::cout << "Scoring parameters:" << std::endl;
	std::cout << "*Gap opening: " << gop << std::endl;
	std::cout << "*Gap extension: " << gex << std::endl;
	std::cout << "*Scoring matrix: " << std::endl;
	TSize alphSize = ValueSize<TAlphabet>::VALUE;
	std::cout << "   ";
	for(TSize col = 0; col<alphSize; ++col) std::cout << TAlphabet(col) << ',';
	std::cout << std::endl;
	for(TSize row = 0; row<alphSize; ++row) {
		for(TSize col = 0; col<alphSize; ++col) {
			if (col == 0) std::cout << TAlphabet(row) << ": ";
			std::cout << score(msaOpt.sc, TAlphabet(row), TAlphabet(col));
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
	String<char> mat;
	if (convertAlignment(g, mat)) {
		TScoreValue alignScore = alignmentEvaluation(g, msaOpt.sc, numGapEx, numGap, numPairs, pairCount, alignLen);
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
	} else {
		std::cout << "No valid alignment!" << std::endl;
	}
}


//////////////////////////////////////////////////////////////////////////////

template<typename TStrSpec, typename TSpec, typename TList, typename TScore, typename TSegmentMatches, typename TScores>
inline void 
__appendSegmentMatches(StringSet<String<AminoAcid, TStrSpec>, Dependent<TSpec> > const& str,
					   TList const& pList,
					   TScore const&,
					   TSegmentMatches& matches,
					   TScores& scores)
{
	Blosum62 local_score(-1,-8);
	appendSegmentMatches(str, pList, local_score, matches, scores, LocalPairwise_Library() );
}

//////////////////////////////////////////////////////////////////////////////

template<typename TValue, typename TStrSpec, typename TSpec, typename TList, typename TScore, typename TSegmentMatches, typename TScores>
inline void 
__appendSegmentMatches(StringSet<String<TValue, TStrSpec>, Dependent<TSpec> > const& str,
					   TList const& pList,
					   TScore const& score_type,
					   TSegmentMatches& matches,
					   TScores& scores)
{
	appendSegmentMatches(str, pList, score_type, matches, scores, LocalPairwise_Library() );
}


//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TCargo, typename TSpec, typename TStringSet1, typename TNames, typename TAlphabet, typename TScore>
inline void
globalMsaAlignment(Graph<Alignment<TStringSet, TCargo, TSpec> >& gAlign, 
				   TStringSet1& sequenceSet,
				   TNames& sequenceNames,
				   MsaOptions<TAlphabet, TScore> const& msaOpt)
{
	typedef typename Value<TScore>::Type TScoreValue;
	typedef typename Size<TStringSet>::Type TSize;
	typedef Graph<Alignment<TStringSet, TSize> > TGraph;
	typedef typename Id<TGraph>::Type TId;
	typedef double TDistanceValue;
	
	
	// Initialize alignment object
	clear(gAlign);
	assignStringSet(gAlign, sequenceSet);

	// Some alignment constants
	TStringSet& seqSet = stringSet(gAlign);
	TSize nSeq = length(seqSet);
	TSize threshold = 30;

	// Select all possible pairs for global and local alignments
	String<TSize> pList;
	selectPairs(seqSet, pList);

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
		for(;begIt != begItEnd; goNext(begIt)) {
			std::ifstream strm_lib;
			strm_lib.open((*begIt).c_str(), ::std::ios_base::in | ::std::ios_base::binary);
			read(strm_lib, matches, scores, sequenceNames, FastaAlign());
			strm_lib.close();
		}
	}

	// Include computed segment matches
	if (!empty(msaOpt.method)) {
		typedef typename Iterator<String<unsigned int>, Standard>::Type TIter;
		TIter begIt = begin(msaOpt.method, Standard() );
		TIter begItEnd = end(msaOpt.method, Standard() );
		for(;begIt != begItEnd; goNext(begIt)) {
			if (*begIt == 0) appendSegmentMatches(seqSet, pList, msaOpt.sc, matches, scores, distanceMatrix, GlobalPairwise_Library() );
			else if (*begIt == 1) __appendSegmentMatches(seqSet, pList, msaOpt.sc, matches, scores);
			else if (*begIt == 2) {
				Nothing noth;
				appendSegmentMatches(seqSet, pList, msaOpt.sc, matches, scores, noth, AlignConfig<true,true,true, true>(), GlobalPairwise_Library() );
			}
			else if (*begIt == 3) appendSegmentMatches(seqSet, pList, matches, scores, Lcs_Library() );
		}	
	}

	// Include a T-Coffee library
	if (!empty(msaOpt.libfiles)) {
		typedef typename Iterator<String<std::string>, Standard>::Type TIter;
		TIter begIt = begin(msaOpt.libfiles, Standard() );
		TIter begItEnd = end(msaOpt.libfiles, Standard() );
		for(;begIt != begItEnd; goNext(begIt)) {
			std::ifstream strm_lib;
			strm_lib.open((*begIt).c_str(), std::ios_base::in | std::ios_base::binary);
			read(strm_lib, matches, scores, sequenceNames, TCoffeeLib());
			strm_lib.close();
		}
	}

	// Include MUMmer segment matches
	if (!empty(msaOpt.mummerfiles)) {
		typedef typename Iterator<String<std::string>, Standard>::Type TIter;
		TIter begIt = begin(msaOpt.mummerfiles, Standard() );
		TIter begItEnd = end(msaOpt.mummerfiles, Standard() );
		for(;begIt != begItEnd; goNext(begIt)) {
			std::ifstream strm_lib;
			strm_lib.open((*begIt).c_str(), std::ios_base::in | std::ios_base::binary);
			read(strm_lib, matches, scores, seqSet, sequenceNames, MummerLib());		
			strm_lib.close();
		}
	}

	// Include BLAST segment matches
	if (!empty(msaOpt.blastfiles)) {
		typedef typename Iterator<String<std::string>, Standard>::Type TIter;
		TIter begIt = begin(msaOpt.blastfiles, Standard() );
		TIter begItEnd = end(msaOpt.blastfiles, Standard() );
		for(;begIt != begItEnd; goNext(begIt)) {
			std::ifstream strm_lib;
			strm_lib.open((*begIt).c_str(), std::ios_base::in | std::ios_base::binary);
			read(strm_lib, matches, scores, sequenceNames, BlastLib());
			strm_lib.close();
		}
	}

	// Use these segment matches for the initial alignment graph
	TGraph g(seqSet);
	if (!msaOpt.rescore) buildAlignmentGraph(matches, scores, g, FractionalScore() );
	else buildAlignmentGraph(matches, scores, g, msaOpt.sc, ReScore() );
	clear(matches);
	clear(scores);

	// Guide tree
	Graph<Tree<TDistanceValue> > guideTree;
	if (!empty(msaOpt.treefile)) {
		std::fstream strm_tree;
		strm_tree.open(msaOpt.treefile.c_str(), ::std::ios_base::in | ::std::ios_base::binary);
		read(strm_tree, guideTree, sequenceNames, NewickFormat());	// Read newick tree
		strm_tree.close();
	} else {
		// Check if we have a valid distance matrix
		if (empty(distanceMatrix)) getDistanceMatrix(g, distanceMatrix, KmerDistance());
		if (msaOpt.build == 0) njTree(distanceMatrix, guideTree);
		else if (msaOpt.build == 1) upgmaTree(distanceMatrix, guideTree, UpgmaMin());
		else if (msaOpt.build == 2) upgmaTree(distanceMatrix, guideTree, UpgmaMax());
		else if (msaOpt.build == 3) upgmaTree(distanceMatrix, guideTree, UpgmaAvg());
		else if (msaOpt.build == 4) upgmaTree(distanceMatrix, guideTree, UpgmaWeightAvg());
	}
	clear(distanceMatrix);
		
	// Triplet extension
	if (nSeq < threshold) tripletLibraryExtension(g);
	else tripletLibraryExtension(g, guideTree, threshold / 2);

	// Progressive Alignment
	progressiveAlignment(g, guideTree, gAlign);

	clear(guideTree);
	clear(g);

	//TStringSet& str = stringSet(gAlign);
	//for(TSize i = 0; i<length(str);++i) {
	//	for(TSize j=0;j<length(str[i]);++j) {
	//		std::cout << ordValue(str[i][j]) << ',';
	//	}
	//	std::cout << std::endl;
	//}
}


//////////////////////////////////////////////////////////////////////////////


template<typename TStringSet, typename TCargo, typename TSpec, typename TScore>
inline void
globalMsaAlignment(Graph<Alignment<TStringSet, TCargo, TSpec> >& gAlign,
				   TScore const& scoreObject)
{
	typedef typename Value<TStringSet>::Type TString;
	typedef typename Value<TString>::Type TAlphabet;
	TStringSet sequenceSet = stringSet(gAlign);
	String<String<char> > sequenceNames;
	fill(sequenceNames, length(sequenceSet), String<char>("tmpName"));
	MsaOptions<AminoAcid, TScore> msaOpt;
	msaOpt.sc = scoreObject;
	appendValue(msaOpt.method, 0);  // Global pairwise
	appendValue(msaOpt.method, 1);	// Local pairwise
	globalMsaAlignment(gAlign, sequenceSet, sequenceNames, msaOpt);
}

//////////////////////////////////////////////////////////////////////////////


template<typename TSource, typename TSpec, typename TScore>
inline void
globalMsaAlignment(Align<TSource, TSpec>& align,
				   TScore const& scoreObject)
{
	typedef StringSet<TSource, Dependent<> > TStringSet;
	TStringSet sequenceSet = stringSet(align);
	Graph<Alignment<TStringSet, void, WithoutEdgeId> > gAlign(sequenceSet);
	globalMsaAlignment(gAlign, scoreObject);
	convertAlignment(gAlign, align);
}


//////////////////////////////////////////////////////////////////////////////
// Just two testing functions
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TStringSet, typename TMatches>
void
_debugMatches(TStringSet& str, 
			  TMatches& matches)
{
	typedef typename Id<TStringSet>::Type TId;
	typedef typename Size<TStringSet>::Type TSize;

	// Print all the matches
	std::cout << "The sequences:" << std::endl;
	for(TSize i = 0;i<length(str);++i) {
		std::cout << positionToId(str,i) << ':' << str[i] << std::endl;
	}
	std::cout << "The matches:" << std::endl;
	for(TSize i = 0;i<length(matches);++i) {
		TId tmp_id1 = sequenceId(matches[i],0);
		std::cout << tmp_id1 << ',' << fragmentBegin(matches[i],tmp_id1) << ',';
		for(TSize j = fragmentBegin(matches[i],tmp_id1); j < fragmentBegin(matches[i],tmp_id1) + fragmentLength(matches[i],tmp_id1); ++j) {
			std::cout << str[idToPosition(str, tmp_id1)][j];
		}
		TId tmp_id2 = sequenceId(matches[i],1);
		std::cout << ',' <<	tmp_id2 << ',' << fragmentBegin(matches[i],tmp_id2) << ',';
		for(TSize j = fragmentBegin(matches[i],tmp_id2); j < fragmentBegin(matches[i],tmp_id2) + fragmentLength(matches[i],tmp_id2); ++j) {
			std::cout << str[idToPosition(str, tmp_id2)][j];
		}
		std::cout << std::endl;

		SEQAN_TASSERT(sequenceId(matches[i],0) != sequenceId(matches[i],1))
		SEQAN_TASSERT(fragmentBegin(matches[i],tmp_id1) < length(str[idToPosition(str, tmp_id1)]))
		SEQAN_TASSERT(fragmentBegin(matches[i],tmp_id1) + fragmentLength(matches[i],tmp_id1) <= length(str[idToPosition(str, tmp_id1)]))
		SEQAN_TASSERT(fragmentBegin(matches[i],tmp_id2) < length(str[idToPosition(str, tmp_id2)]))
		SEQAN_TASSERT(fragmentBegin(matches[i],tmp_id2) + fragmentLength(matches[i],tmp_id2) <= length(str[idToPosition(str, tmp_id2)]))
		SEQAN_TASSERT(fragmentLength(matches[i],tmp_id2) == fragmentLength(matches[i],tmp_id1))
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TGraph>
void
_debugRefinedMatches(TGraph& g)
{
	typedef typename Id<TGraph>::Type TId;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
	
	std::cout << "Refined matches" << std::endl;
	TEdgeIterator it_tmp(g);
	for(;!atEnd(it_tmp);++it_tmp) {
		TId id1 = sequenceId(g,sourceVertex(it_tmp));
		TId id2 = sequenceId(g,targetVertex(it_tmp));
		std::cout << id1 << ',' << fragmentBegin(g,sourceVertex(it_tmp)) << ',';
		std::cout << label(g,sourceVertex(it_tmp));
		std::cout << ',' <<	id2 << ',' << fragmentBegin(g,targetVertex(it_tmp)) << ',';
		std::cout << label(g,targetVertex(it_tmp));
		std::cout << " (" << cargo(*it_tmp) << ")";
		std::cout << std::endl;	

		SEQAN_TASSERT(sequenceId(g,sourceVertex(it_tmp)) != sequenceId(g,targetVertex(it_tmp)))
		SEQAN_TASSERT(fragmentBegin(g,sourceVertex(it_tmp)) < length((stringSet(g))[idToPosition((stringSet(g)), id1)]))
		SEQAN_TASSERT(fragmentBegin(g,sourceVertex(it_tmp)) + fragmentLength(g,sourceVertex(it_tmp)) <= length((stringSet(g))[idToPosition((stringSet(g)), id1)]))
		SEQAN_TASSERT(fragmentBegin(g,targetVertex(it_tmp)) < length((stringSet(g))[idToPosition((stringSet(g)), id2)]))
		SEQAN_TASSERT(fragmentBegin(g,targetVertex(it_tmp)) + fragmentLength(g,targetVertex(it_tmp)) <= length((stringSet(g))[idToPosition((stringSet(g)), id2)]))
		SEQAN_TASSERT(fragmentLength(g,sourceVertex(it_tmp)) == fragmentLength(g,targetVertex(it_tmp)))

	}
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...


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
==========================================================================*/

#ifndef SEQAN_HEADER_SEQAN_TCOFFEE_H
#define SEQAN_HEADER_SEQAN_TCOFFEE_H

#include <seqan/graph_msa.h>

#ifdef SEQAN_PROFILE
SEQAN_PROTIMESTART(__myProfileTime); // Profiling
#endif

using namespace seqan;

//////////////////////////////////////////////////////////////////////////////////

template<typename TConfigOptions, typename TScore, typename TAlphabet>
inline int
evaluateAlignment(TConfigOptions const& cfgOpt, TScore const& scType, TAlphabet) {
	typedef typename Value<TScore>::Type TScoreValue;
	typedef String<TAlphabet> TSequence;
	typedef typename Size<TSequence>::Type TSize;
	typedef String<char> TName;
	StringSet<TSequence, Owner<> > origStrSet;
	StringSet<TName> names;

	// Read the sequences
	std::fstream strm;
	strm.open(toCString(value(cfgOpt, "infile")), std::ios_base::in | std::ios_base::binary);
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
	strm_lib.open(toCString(value(cfgOpt, "infile")), std::ios_base::in | std::ios_base::binary);
	read(strm_lib,matches, scores, names, FastaAlign());	
	strm_lib.close();

	// Build the alignment graph
	typedef Graph<Alignment<TDepSequenceSet, TSize> > TGraph;
	TGraph g(strSet);
	buildAlignmentGraph(matches, g, FrequencyCounting() );

	// Print the scoring information
	TScoreValue gop = scoreGapOpen(scType);
	TScoreValue gex = scoreGap(scType);
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
			std::cout << score(scType, TAlphabet(row), TAlphabet(col));
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

	return 0;

}

//////////////////////////////////////////////////////////////////////////////////

template<typename TConfigOptions>
inline int
evaluateAlignment(TConfigOptions const& cfgOpt, AminoAcid) 
{
	if (!length(value(cfgOpt, "matrix"))) {
		typedef typename Value<Blosum62>::Type TScoreValue;
		Blosum62 sc(-1 * _stringToNumber<TScoreValue>(value(cfgOpt, "gex")) , -1 * _stringToNumber<TScoreValue>(value(cfgOpt, "gop")));
		return evaluateAlignment(cfgOpt, sc, AminoAcid() );
	} else {
		typedef Score<int, ScoreMatrix<> > TScore;
		typedef typename Value<TScore>::Type TScoreValue;
		TScore sc;
		loadScoreMatrix(sc, value(cfgOpt, "matrix"));
		sc.data_gap_extend = -1 * _stringToNumber<TScoreValue>(value(cfgOpt, "gex")); 
		sc.data_gap_open = -1 * _stringToNumber<TScoreValue>(value(cfgOpt, "gop"));
		return evaluateAlignment(cfgOpt, sc, AminoAcid() );
	}
}

//////////////////////////////////////////////////////////////////////////////////

template<typename TConfigOptions>
inline int
evaluateAlignment(TConfigOptions const& cfgOpt, Dna5) 
{
	typedef Score<int> TScore;
	typedef typename Value<TScore>::Type TScoreValue;
	TScore scType(_stringToNumber<TScoreValue>(value(cfgOpt, "msc")),_stringToNumber<TScoreValue>(value(cfgOpt, "mmsc")),-1 * _stringToNumber<TScoreValue>(value(cfgOpt, "gex")),-1 * _stringToNumber<TScoreValue>(value(cfgOpt, "gop")));
	return evaluateAlignment(cfgOpt, scType, Dna5() );
}

//////////////////////////////////////////////////////////////////////////////////

template<typename TConfigOptions>
inline int
evaluateAlignment(TConfigOptions const& cfgOpt, Rna5) 
{
	typedef Score<int> TScore;
	typedef typename Value<TScore>::Type TScoreValue;
	TScore scType(_stringToNumber<TScoreValue>(value(cfgOpt, "msc")),_stringToNumber<TScoreValue>(value(cfgOpt, "mmsc")),-1 * _stringToNumber<TScoreValue>(value(cfgOpt, "gex")),-1 * _stringToNumber<TScoreValue>(value(cfgOpt, "gop")));
	return evaluateAlignment(cfgOpt, scType, Rna5() );
}


//////////////////////////////////////////////////////////////////////////////

template<typename TString, typename TSpec, typename TName, typename TSpec2, typename TConfigOptions, typename TAlignmentGraph>
inline void
predefinedGlobalAlignment(StringSet<TString, TSpec> const& seqSet,
						  StringSet<TName, TSpec2> const& nameSet,
						  TConfigOptions& cfgOpt,
						  TAlignmentGraph& gOut)
{
	typedef StringSet<TString, Dependent<> > TStringSet;
	typedef typename Size<TStringSet>::Type TSize;
	typedef Graph<Alignment<TStringSet, TSize> > TGraph;
	typedef typename Id<TGraph>::Type TId;
	TSize nSeq = length(seqSet);
	TSize threshold = 30;

	// Select all pairs
	String<Pair<TId, TId> > pList;
	selectPairs(seqSet, pList);

	// Set-up a distance matrix
	typedef String<double> TDistanceMatrix;
	typedef typename Value<TDistanceMatrix>::Type TDistanceValue;
	TDistanceMatrix distanceMatrix;

	// Set-up alignment scoring matrices
	typedef Blosum62 TScore;
	typedef typename Value<TScore>::Type TScoreValue;
	TScore score_type_global(-1,-11);
	TScore score_type_local(-2,-8);

	// Containers for segment matches and corresponding scores 
	typedef String<Fragment<> > TFragmentString;
	TFragmentString matches;
	typedef String<TScoreValue> TScoreValues;
	TScoreValues scores;

	// Compute segment matches from global pairwise alignments
	appendSegmentMatches(seqSet, pList, score_type_global, matches, scores, distanceMatrix, GlobalPairwise_Library() );
	
	// Append segment matches from local pairwise alignments
	appendSegmentMatches(seqSet, pList, score_type_local, matches, scores, LocalPairwise_Library() );

	// Select a subset of pairs for end-gaps free
	typedef std::multimap<TDistanceValue, Pair<TId, TId> > TBestPairs;
	TBestPairs bestPairs;
	TDistanceValue dist = 0;
	for(TSize i=0;i<nSeq-1;++i) {
		for(TSize j=i+1;j<nSeq;++j) {
			TDistanceValue d = value(distanceMatrix, i*nSeq+j);
			bestPairs.insert(std::make_pair(d, Pair<TId, TId>(positionToId(seqSet, i),positionToId(seqSet, j)))); 
			dist+=d;
		}
	}
	TSize numPairs = nSeq * (nSeq - 1) / 2;
	dist /= numPairs;
	String<Pair<TId, TId> > pListLocal;
	pListLocal = pList;
	if (nSeq > threshold) {
		clear(pListLocal);
		typename TBestPairs::reverse_iterator itBestR = bestPairs.rbegin();
		typename TBestPairs::iterator itBestF = bestPairs.begin();
		TSize limit = 4 * nSeq;
		for(TSize i=0;i<limit;++i, ++itBestR, ++itBestF) {
			appendValue(pListLocal, itBestF->second);
			appendValue(pListLocal, itBestR->second);
		}
	}

	// Append segment matches from end-gaps free alignments
	Nothing noth;
	if (dist > 0.75) {
		Blosum30 sT(-4,-20);
		AlignConfig<true,true,true,true> ac;
		appendSegmentMatches(seqSet, pListLocal, sT, matches, scores, noth, ac, GlobalPairwise_Library() );
	} else if (dist < 0.50) {
		Blosum80 sT(-2, -12);
		AlignConfig<true,true,true,true> ac;
		appendSegmentMatches(seqSet, pListLocal, sT, matches, scores, noth, ac, GlobalPairwise_Library() );
	} else {
		Blosum62 sT(-3,-14);
		AlignConfig<true,true,true,true> ac;
		appendSegmentMatches(seqSet, pListLocal, sT, matches, scores, noth, ac, GlobalPairwise_Library() );
	}

	// Include external segment matches
	if (length(value(cfgOpt, "blast"))) {
		std::fstream strm_lib;
		strm_lib.open(toCString(value(cfgOpt, "blast")), std::ios_base::in | std::ios_base::binary);
		read(strm_lib, matches, scores, nameSet, BlastLib());
		strm_lib.close();
	}

	// Score the matches
	scoreMatches(seqSet, score_type_global, matches, scores);
	
	// Use these segment matches for the initial alignment graph
	TGraph g(seqSet);
	buildAlignmentGraph(matches, scores, g, FrequencyCounting() );
	
	// Guide tree
	Graph<Tree<double> > guideTree;
	if (length(value(cfgOpt, "usetree"))) {
		std::fstream strm_tree;
		strm_tree.open(toCString(value(cfgOpt, "usetree")), std::ios_base::in | std::ios_base::binary);
		read(strm_tree,guideTree,nameSet,NewickFormat());	// Read newick tree
		strm_tree.close();
	} else {
		//upgmaTree(distanceMatrix, guideTree);
		slowNjTree(distanceMatrix, guideTree);
	}

	// Triplet extension and progressive alignment
	if (nSeq < threshold) {
		// Full triplet...
		tripletLibraryExtension(g);

		// ... and normal progressive alignment with guide tree
		progressiveAlignment(g, guideTree, gOut);
	} else {
		// Triplet only on groups of sequences
		progressiveAlignment(g, guideTree, gOut, threshold);
	}

	clear(guideTree);
	clear(distanceMatrix);
	clear(g);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TString, typename TSpec, typename TName, typename TSpec2, typename TScoreMatrix, typename TConfigOptions, typename TAlignmentGraph>
inline void
globalAlignment(StringSet<TString, TSpec> const& seqSet,
				StringSet<TName, TSpec2> const& nameSet,
				TConfigOptions const& cfgOpt,
				TScoreMatrix const& scType,
				TAlignmentGraph& gOut)
{
	typedef typename Value<TScoreMatrix>::Type TScoreValue;
	typedef StringSet<TString, Dependent<> > TStringSet;
	typedef typename Size<TStringSet>::Type TSize;
	typedef Graph<Alignment<TStringSet, TSize> > TGraph;
	typedef typename Id<TGraph>::Type TId;
	TSize nSeq = length(seqSet);
	TSize threshold = 30;

	// Select all pairs
	String<Pair<TId, TId> > pList;
	selectPairs(seqSet, pList);

	// Set-up a distance matrix
	typedef String<double> TDistanceMatrix;
	typedef typename Value<TDistanceMatrix>::Type TDistanceValue;
	TDistanceMatrix distanceMatrix;

	// Containers for segment matches and corresponding scores 
	typedef String<Fragment<> > TFragmentString;
	TFragmentString matches;
	typedef String<TScoreValue> TScoreValues;
	TScoreValues scores;

	// Include segment matches from other alignments
	if (length(value(cfgOpt, "aln"))) {
#ifdef SEQAN_PROFILE
		std::cout << "Parsing external alignment files:" << std::endl;
#endif
		String<char> alnFiles = value(cfgOpt, "aln");
		String<char> alignmentFile;
		for(TSize i = 0; i<=length(alnFiles); ++i) {
			if ((i == length(alnFiles) || (value(alnFiles, i) == ','))) {		
#ifdef SEQAN_PROFILE
				std::cout << "*Alignment file: " << alignmentFile << std::endl;
#endif
				std::stringstream input;
				input << alignmentFile;
				std::ifstream strm_lib;
				strm_lib.open(input.str().c_str(), std::ios_base::in | std::ios_base::binary);
				read(strm_lib, matches, scores, nameSet, FastaAlign());
				strm_lib.close();
				clear(alignmentFile);
			} else {
				if ((value(alnFiles, i) != ' ') && (value(alnFiles, i) != '\t')) appendValue(alignmentFile, value(alnFiles, i));
			}
		}
#ifdef SEQAN_PROFILE
		std::cout << "External segment matches done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;
#endif
	}

	// Include computed segment matches
	if (length(value(cfgOpt, "method"))) {
#ifdef SEQAN_PROFILE
		std::cout << "Computing segment matches:" << std::endl;
#endif
		String<char> methodNames = value(cfgOpt, "method");
		String<char> currentMethod;
		for(TSize i = 0; i<=length(methodNames); ++i) {
			if ((i == length(methodNames) || (value(methodNames, i) == ','))) {		
#ifdef SEQAN_PROFILE
				std::cout << "*Method: " << currentMethod << std::endl;
#endif
				if (currentMethod == "global") {
					// Compute segment matches from global pairwise alignments
					appendSegmentMatches(seqSet, pList, scType, matches, scores, distanceMatrix, GlobalPairwise_Library() );
				} else if (currentMethod == "local") {
					// Compute segment matches from local pairwise alignments
					appendSegmentMatches(seqSet, pList, scType, matches, scores, LocalPairwise_Library() );
					std::cout << "Hi " << std::endl;
				} else if (currentMethod == "overlap") {
					// Compute segment matches from overlap alignments
					Nothing noth;
					AlignConfig<true,true,false,false> ac;
					appendSegmentMatches(seqSet, pList, scType, matches, scores, noth, ac, GlobalPairwise_Library() );
				} else if (currentMethod == "lcs") {
					// Compute segment matches from the longest common subsequence
					appendSegmentMatches(seqSet, matches, scores, Lcs_Library() );
				} else {
#ifdef SEQAN_PROFILE
					std::cout << "*Unknown method!!!" << std::endl;
#endif
				}
#ifdef SEQAN_PROFILE
				std::cout << "*Done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;
#endif
				clear(currentMethod);
			} else {
				if ((value(methodNames, i) != ' ') && (value(methodNames, i) != '\t')) appendValue(currentMethod, value(methodNames, i));
			}
		}	
	}

	// Include a T-Coffee library
	if (length(value(cfgOpt, "lib"))) {
#ifdef SEQAN_PROFILE
		std::cout << "Parsing a T-Coffee Library:" << std::endl;
#endif
		String<char> libNames = value(cfgOpt, "lib");
		String<char> currentLib;
		for(TSize i = 0; i<=length(libNames); ++i) {
			if ((i == length(libNames) || (value(libNames, i) == ','))) {
#ifdef SEQAN_PROFILE
				std::cout << "*T-Coffee library: " << currentLib << std::endl;
#endif
				std::stringstream input;
				input << currentLib;
				std::ifstream strm_lib;
				strm_lib.open(input.str().c_str(), std::ios_base::in | std::ios_base::binary);
				read(strm_lib, matches, scores, nameSet, TCoffeeLib());
				strm_lib.close();
				clear(currentLib);
			} else {
				if ((value(libNames, i) != ' ') && (value(libNames, i) != '\t')) appendValue(currentLib, value(libNames, i));
			}
		}
#ifdef SEQAN_PROFILE
		std::cout << "Parsing done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;
#endif
	}
	
	// Include MUMmer segment matches
	if (length(value(cfgOpt, "mummer"))) {
#ifdef SEQAN_PROFILE
		std::cout << "Parsing MUMmer segment matches:" << std::endl;
#endif
		String<char> mummerFiles = value(cfgOpt, "mummer");
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
				read(strm_lib, matches, scores, seqSet, nameSet, MummerLib());		
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

	// Include BLAST segment matches
	if (length(value(cfgOpt, "blast"))) {
#ifdef SEQAN_PROFILE
		std::cout << "Parsing BLAST segment matches:" << std::endl;
#endif
		String<char> blastFiles = value(cfgOpt, "blast");
		String<char> currentBlastFile;
		for(TSize i = 0; i<=length(blastFiles); ++i) {
			if ((i == length(blastFiles) || (value(blastFiles, i) == ','))) {		
#ifdef SEQAN_PROFILE
				std::cout << "*BLAST file: " << currentBlastFile << std::endl;
#endif
				std::stringstream input;
				input << currentBlastFile;
				std::ifstream strm_lib;
				strm_lib.open(input.str().c_str(), std::ios_base::in | std::ios_base::binary);
				read(strm_lib, matches, scores, nameSet, BlastLib());
				strm_lib.close();
				clear(currentBlastFile);
			} else {
				if ((value(blastFiles, i) != ' ') && (value(blastFiles, i) != '\t')) appendValue(currentBlastFile, value(blastFiles, i));
			}
		}
#ifdef SEQAN_PROFILE
		std::cout << "Parsing done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;
#endif
	}

#ifdef SEQAN_PROFILE
	std::cout << "Total number of segment matches: " << length(matches) << std::endl;
#endif
	// Score the matches
	if ((value(cfgOpt, "rescore") == "true") && (length(value(cfgOpt, "seq")))) {
#ifdef SEQAN_PROFILE	
		std::cout << "Scoring method: Re-Score" << std::endl;
#endif
		scoreMatches(seqSet, scType, matches, scores);
#ifdef SEQAN_PROFILE
		std::cout << "Scoring done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;
#endif
	}
	
	// Use these segment matches for the initial alignment graph
	TGraph g(seqSet);
#ifdef SEQAN_PROFILE
	std::cout << "Construction of alignment graph: FractionalScore" << std::endl;
#endif
	// Build the alignment graph
	// Make a dependent StringSet
	typedef Dna5 TAlphabet;
	typedef StringSet<TString, TSpec> TSequenceSet;
	typedef String<TAlphabet> TSequence;
	typedef StringSet<TSequence, Dependent<> > TDepSequenceSet;
	//TDepSequenceSet seqSet2(seqSet);
	//typedef Graph<Alignment<TSequenceSet, TSize> > TGraph2;
	//TGraph2 g2(seqSet);
	//buildAlignmentGraph(matches, g2, FrequencyCounting() );
	buildAlignmentGraph(matches, scores, g, FractionalScore() );
	//#############
	//Scoring code
	String<char> align;
	std::cout << convertAlignment(g2,align) << std::endl;
	TScoreValue gop = scoreGapOpen(scType);
	TScoreValue gex = scoreGap(scType);
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
	    std::cout << score(scType, TAlphabet(row), TAlphabet(col));
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
	TScoreValue alignScore = alignmentEvaluation(g2, scType, numGapEx, numGap, numPairs, pairCount, alignLen);
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
	//########3
	clear(matches);
	clear(scores);

	std::fstream strm;
	strm.open("sam.out", std::ios_base::out | std::ios_base::trunc);
	write(strm,g,nameSet,Raw());
	strm.close();
#ifdef SEQAN_PROFILE
	std::cout << "Alignment graph construction done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;
#endif

	// Guide tree
	Graph<Tree<double> > guideTree;
	if (length(value(cfgOpt, "usetree"))) {
#ifdef SEQAN_PROFILE
		std::cout << "Guide Tree: " << value(cfgOpt, "usetree") << std::endl;
#endif
		std::fstream strm_tree;
		strm_tree.open(toCString(value(cfgOpt, "usetree")), std::ios_base::in | std::ios_base::binary);
		read(strm_tree,guideTree,nameSet,NewickFormat());	// Read newick tree
		strm_tree.close();
	} else {
#ifdef SEQAN_PROFILE
		std::cout << "Guide Tree: Neighbor Joining" << std::endl;
#endif
		// Check if we have a valid distance matrix
		if (empty(distanceMatrix)) getDistanceMatrix(g, distanceMatrix, LibraryDistance());
		slowNjTree(distanceMatrix, guideTree);
	}
#ifdef SEQAN_PROFILE
	std::cout << "Guide tree done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;
#endif
	
	// Triplet extension
	if (nSeq < threshold) tripletLibraryExtension(g);
	else groupBasedTripletExtension(g, guideTree, threshold / 2);
#ifdef SEQAN_PROFILE
	std::cout << "Triplet extension done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;
#endif

	// Progressive Alignment
	progressiveAlignment(g, guideTree, gOut);
#ifdef SEQAN_PROFILE
	std::cout << "Progressive alignment done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;
#endif

	clear(guideTree);
	clear(distanceMatrix);
	clear(g);
#ifdef SEQAN_PROFILE
	std::cout << "Clean-up done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;
#endif
}

//////////////////////////////////////////////////////////////////////////////

template<typename TString, typename TSpec, typename TName, typename TSpec2, typename TConfigOptions, typename TAlignmentGraph>
inline void
customizedMsaAlignment(StringSet<TString, TSpec> const& seqSet,
					   StringSet<TName, TSpec2> const& nameSet,
					   TConfigOptions& cfgOpt,
					   TAlignmentGraph& gOut,
					   AminoAcid)
{
#ifdef SEQAN_PROFILE
	std::cout << "Scoring: " << std::endl;
#endif
	if (!length(value(cfgOpt, "matrix"))) {
		typedef typename Value<Blosum62>::Type TScoreValue;
		TScoreValue gopening = -1 * _stringToNumber<TScoreValue>(value(cfgOpt, "gop")); 
		TScoreValue gextend = -1 * _stringToNumber<TScoreValue>(value(cfgOpt, "gex")); 
#ifdef SEQAN_PROFILE
		std::cout << "*Gap open penalty: " << gopening << std::endl;
		std::cout << "*Gap extension penalty: " << gextend << std::endl;
		std::cout << "*Scoring Matrix: Blosum62" << std::endl;
#endif
		Blosum62 sc(gextend , gopening);
		globalAlignment(seqSet, nameSet, cfgOpt, sc, gOut);
	} else {
		typedef Score<int, ScoreMatrix<> > TScore;
		typedef typename Value<TScore>::Type TScoreValue;
		TScore sc;
		TScoreValue gopening = -1 * _stringToNumber<TScoreValue>(value(cfgOpt, "gop")); 
		TScoreValue gextend = -1 * _stringToNumber<TScoreValue>(value(cfgOpt, "gex")); 
#ifdef SEQAN_PROFILE
		std::cout << "*Gap open penalty: " << gopening << std::endl;
		std::cout << "*Gap extension penalty: " << gextend << std::endl;
		std::cout << "*Scoring Matrix: " << value(cfgOpt, "matrix") << std::endl;
#endif
		loadScoreMatrix(sc, value(cfgOpt, "matrix"));
		sc.data_gap_extend = (TScoreValue) (gextend);
		sc.data_gap_open = (TScoreValue) (gopening);
		globalAlignment(seqSet, nameSet, cfgOpt, sc, gOut);
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TString, typename TSpec, typename TName, typename TSpec2, typename TConfigOptions, typename TAlignmentGraph>
inline void
customizedMsaAlignment(StringSet<TString, TSpec> const& seqSet,
					   StringSet<TName, TSpec2> const& nameSet,
					   TConfigOptions& cfgOpt,
					   TAlignmentGraph& gOut,
					   Dna5)
{
	typedef Score<int> TScore;
	typedef typename Value<TScore>::Type TScoreValue;
	TScoreValue gop = -1 * _stringToNumber<TScoreValue>(value(cfgOpt, "gop")); 
	TScoreValue gex = -1 * _stringToNumber<TScoreValue>(value(cfgOpt, "gex")); 
	TScoreValue match = _stringToNumber<TScoreValue>(value(cfgOpt, "msc")); 
	TScoreValue mismatch = _stringToNumber<TScoreValue>(value(cfgOpt, "mmsc")); 
#ifdef SEQAN_PROFILE
	std::cout << "Scoring: " << std::endl;
	std::cout << "*Gap open penalty: " << gop << std::endl;
	std::cout << "*Gap extension penalty: " << gex << std::endl;
	std::cout << "*Match score: " << match << std::endl;
	std::cout << "*Mismatch score: " << mismatch << std::endl;
#endif
	TScore scType(match,mismatch,gex,gop);
	globalAlignment(seqSet, nameSet, cfgOpt, scType, gOut);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TString, typename TSpec, typename TName, typename TSpec2, typename TConfigOptions, typename TAlignmentGraph>
inline void
customizedMsaAlignment(StringSet<TString, TSpec> const& seqSet,
					   StringSet<TName, TSpec2> const& nameSet,
					   TConfigOptions& cfgOpt,
					   TAlignmentGraph& gOut,
					   Rna5)
{
	typedef Score<int> TScore;
	typedef typename Value<TScore>::Type TScoreValue;
	TScoreValue gop = -1 * _stringToNumber<TScoreValue>(value(cfgOpt, "gop")); 
	TScoreValue gex = -1 * _stringToNumber<TScoreValue>(value(cfgOpt, "gex")); 
	TScoreValue match = _stringToNumber<TScoreValue>(value(cfgOpt, "msc")); 
	TScoreValue mismatch = _stringToNumber<TScoreValue>(value(cfgOpt, "mmsc")); 
#ifdef SEQAN_PROFILE
	std::cout << "Scoring: " << std::endl;
	std::cout << "*Gap open penalty: " << gop << std::endl;
	std::cout << "*Gap extension penalty: " << gex << std::endl;
	std::cout << "*Match score: " << match << std::endl;
	std::cout << "*Mismatch score: " << mismatch << std::endl;
#endif
	TScore scType(match,mismatch,gex,gop);
	globalAlignment(seqSet, nameSet, cfgOpt, scType, gOut);
}

//////////////////////////////////////////////////////////////////////////////////

template<typename TConfigOptions, typename TAlphabet>
inline int
customizedMsaAlignment(TConfigOptions const& cfgOpt, TAlphabet) {

	//////////////////////////////////////////////////////////////////////////////
	// Read the sequences
	//////////////////////////////////////////////////////////////////////////////
	typedef String<TAlphabet> TSequence;
	typedef String<char> TName;
	StringSet<TSequence, Owner<> > origStrSet;
	StringSet<TName> names;

	// Read the sequences from a T-Coffee library file or a regular sequence file in FASTA format
	if ((!length(value(cfgOpt, "seq"))) && (length(value(cfgOpt, "lib")))) {
		std::fstream strm;
		strm.open(toCString(value(cfgOpt, "lib")), std::ios_base::in | std::ios_base::binary);
		read(strm,origStrSet,names,TCoffeeLib());	
		strm.close();			
	} else if (length(value(cfgOpt, "seq"))) {
		_loadSequences(value(cfgOpt, "seq"), origStrSet, names);
	} else {
		std::cerr << valueHelp(cfgOpt) << std::endl;
		return -1;
	}
	typedef StringSet<TSequence, Dependent<> > TDepSequenceSet;
	typedef typename Size<TDepSequenceSet>::Type TSize;
	TDepSequenceSet strSet(origStrSet);
#ifdef SEQAN_PROFILE
	std::cout << "Number of sequences: " << length(strSet) << std::endl;
	std::cout << "Import of sequences done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;
#endif

	//////////////////////////////////////////////////////////////////////////////
	// Alignment of the sequences
	//////////////////////////////////////////////////////////////////////////////
	typedef Graph<Alignment<TDepSequenceSet, TSize> > TGraph;
	Graph<Alignment<TDepSequenceSet, void, WithoutEdgeId> > gOut(strSet);
	
	// Predefined alignment pipeline
	//predefinedGlobalAlignment(strSet, names, cfgOpt, gOut);
	// Custom Alignment
	customizedMsaAlignment(strSet, names, cfgOpt, gOut, TAlphabet() );

	//////////////////////////////////////////////////////////////////////////////
	// Alignment output
	//////////////////////////////////////////////////////////////////////////////

	if (value(cfgOpt, "output") == "fasta") {
		std::fstream strm;
		strm.open(toCString(value(cfgOpt, "outfile")), std::ios_base::out | std::ios_base::trunc);
		write(strm,gOut,names,FastaFormat());
		strm.close();
	} else if (value(cfgOpt, "output") == "msf") {
		std::fstream strm;
		strm.open(toCString(value(cfgOpt, "outfile")), std::ios_base::out | std::ios_base::trunc);
		write(strm,gOut,names,MsfFormat());
		strm.close();
	}
#ifdef SEQAN_PROFILE
	std::cout << "Output done: " << SEQAN_PROTIMEUPDATE(__myProfileTime) << " seconds" << std::endl;
#endif

	return 0;
}

#endif


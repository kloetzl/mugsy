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
  $Id: graph_consensus_library.h 1809 2008-03-31 12:57:59Z rausch@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_GRAPH_CONSENSUS_LIBRARY_H
#define SEQAN_HEADER_GRAPH_CONSENSUS_LIBRARY_H

namespace SEQAN_NAMESPACE_MAIN
{



//////////////////////////////////////////////////////////////////////////////

template<typename TFragment, typename TSpec1, typename TStringSet, typename TSize>
inline void 
getAlignmentStatistics(String<TFragment, TSpec1>& matches,
					   TStringSet& str,
					   TSize& from,
					   TSize& matchLength,
					   TSize& overlapLength,
					   TSize& alignLength)
{
	SEQAN_CHECKPOINT
	getAlignmentStatistics(matches, str, (TSize) from, (TSize) length(matches), matchLength, overlapLength, alignLength, typename Value<TStringSet>::Type() );
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSize>
inline void 
__getAlignmentStatistics(Nothing&,
						 TSize,
						 TSize,
						 TSize,
						 double,
						 double)
{
	SEQAN_CHECKPOINT
}

//////////////////////////////////////////////////////////////////////////////

template<typename TValue, typename TSpec, typename TSize>
inline void 
__getAlignmentStatistics(String<TValue, TSpec>& dist,
						 TSize i,
						 TSize j,
						 TSize nseq,
						 double matchLen,
						 double quality)
{
	SEQAN_CHECKPOINT
	TValue sim = matchLen * matchLen * quality;
	if (sim > 10000) sim = 10000;
	assignValue(dist, i*nseq + j, 10000 - sim);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TCargo, typename TSpec, typename TSize>
inline void 
__getAlignmentStatistics(Graph<Undirected<TCargo, TSpec> >& dist,
						 TSize i,
						 TSize j,
						 TSize,
						 double matchLen,
						 double quality)
{
	SEQAN_CHECKPOINT
	TCargo sim = matchLen * matchLen * quality;
	if (sim > 10000) sim = 10000;
	addEdge(dist, i, j, 10000 - sim);
}

//////////////////////////////////////////////////////////////////////////////
// Layout-based pair selection
//////////////////////////////////////////////////////////////////////////////

template<typename TSize>
struct _LessPair :
	public ::std::unary_function<Pair<TSize, TSize>, bool>
{
	inline bool 
	operator() (Pair<TSize, TSize> const& a1, Pair<TSize, TSize> const& a2) const {
		if (a1.i1 == a2.i1) return (a1.i2 < a2.i2);
		else return (a1.i1 < a2.i1);
	}
};

template<typename TSize>
struct _LessTripel :
	public ::std::unary_function<Pair<TSize, Triple<TSize, TSize, TSize> >, bool>
{
	inline bool 
	operator() (Pair<TSize, Triple<TSize, TSize, TSize> > const& a1, Pair<TSize, Triple<TSize, TSize, TSize> > const& a2) {
		return (a1.i1 < a2.i1);
	}
};


//////////////////////////////////////////////////////////////////////////////

template<typename TString, typename TSpec, typename TBegEndPos, typename TSize, typename TPairList, typename TPos, typename TSpec2>
inline void 
selectPairs(StringSet<TString, TSpec> const& str,
			TBegEndPos const& begEndPos,
			TSize bandwidth,							
			TPairList& pList,
			String<Pair<TPos, TPos>, TSpec2>& dList)
{
	SEQAN_CHECKPOINT
	typedef String<Pair<TPos, TPos>, TSpec2>  TDistanceList;
	typedef StringSet<TString, TSpec> TStringSet;
	typedef Pair<TPos, TPos> TDiagPair;
	typedef typename Value<TPairList>::Type TPair;
	typedef typename Iterator<TPairList>::Type TPairIter;
	typedef typename Iterator<TBegEndPos>::Type TBegEndIter;

	// Workaround for strange celera behaviour (just for contained reads)
#ifdef CELERA_OFFSET
	TSize contained_offset=200;
#else
	TSize contained_offset=0;
#endif
	
	// Sort the reads by their first index position
	TBegEndIter begEndIt = begin(begEndPos);
	TBegEndIter begEndItEnd = end(begEndPos);
	typedef Triple<TSize, TSize, TSize> TInfo;
	typedef String<Pair<TSize, TInfo> > TPosIndexList;
	typedef typename Iterator<TPosIndexList>::Type TPosIter;
	TPosIndexList posIndex;
	resize(posIndex, length(begEndPos));
	TPosIter posIndexIt = begin(posIndex);
	TPosIter posIndexItEnd = end(posIndex);
	for(TSize index = 0;begEndIt != begEndItEnd; goNext(begEndIt), goNext(posIndexIt), ++index) {
		TSize pos1 = (value(begEndIt)).i1;
		TSize pos2 = (value(begEndIt)).i2;
		if (pos1 < pos2) value(posIndexIt) = Pair<TSize, TInfo>(pos1,TInfo(index, pos1, pos2));
		else value(posIndexIt) = Pair<TSize, TInfo>(pos2,TInfo(index, pos1, pos2));
	}
	std::sort(begin(posIndex, Standard() ), end(posIndex, Standard() ), _LessTripel<TSize>() );

	// The expected overlap by a pair of reads (represented by its index)
	typedef String<Pair<TSize, TSize> > TOverlapIndexList;
	TOverlapIndexList ovlIndex;
	TSize pairLen = 0;  // Pair Counter
	TPos const initialRadius = (bandwidth + 1) / 2;
	TPos const lengthDivider = 5;	// Overlap / 2^lengthDivider is added to the radius

	// Find all overlapping reads
	TSize nseq = length(str);
	TDistanceList preDList;
	TPairList prePList;
	reserve(preDList, nseq * 40);
	reserve(prePList, nseq * 40);
	reserve(ovlIndex, nseq * 40);
	posIndexIt = begin(posIndex);
	posIndexItEnd = end(posIndex);
	for(;posIndexIt != posIndexItEnd; goNext(posIndexIt)) {
		TSize index1 = ((value(posIndexIt)).i2).i1;
		TPos posIi1 = ((value(posIndexIt)).i2).i2;
		TPos posIi2 = ((value(posIndexIt)).i2).i3;
		TSize lenI = 0;
		bool forwardI = true;
		if (posIi1 < posIi2) lenI = posIi2 - posIi1;
		else {
			lenI = posIi1 - posIi2;
			forwardI = false;
		}
		TPosIter posIndexIt2 = posIndexIt;
		goNext(posIndexIt2);
		for(;posIndexIt2 != posIndexItEnd; goNext(posIndexIt2)) {
			if ((value(posIndexIt)).i1 + lenI <= (value(posIndexIt2)).i1) break;
			TSize index2 = ((value(posIndexIt2)).i2).i1;
			TPos posJi1 = ((value(posIndexIt2)).i2).i2;
			TPos posJi2 = ((value(posIndexIt2)).i2).i3;

			// Diagonal boundaries of the band
			// Initialization values are used if one read is contained in the other 
			TSize lenJ = 0;
			bool forwardJ = true;
			if (posJi1 < posJi2) lenJ = posJi2 - posJi1;
			else {
				forwardJ = false;
				lenJ = posJi1 - posJi2;
			}
			TPos diagLow = -1 * (TPos) lenJ;
			TPos diagHigh = (TPos) lenI;
			TPos radius = initialRadius;		// Increased by overlap length

			// Read orientations
			if (forwardI) {
				// 1) Forward - Forward
				if (forwardJ) {
					if ((posJi2 < posIi2) && (posJi1 < posIi1)) {
						TPos offset = (posIi1 - posJi1);
						radius += (posJi2 - posJi1 - offset) >> lengthDivider;
						if (-1 * offset - radius > diagLow) diagLow = -1 * offset - radius;
						if (-1 * offset + radius < diagHigh) diagHigh = -1 * offset + radius;
					} else if ((posJi1 > posIi1) && (posJi2 > posIi2)) {
						TPos offset = (posJi1 - posIi1);
						radius += (posIi2 - posIi1 - offset) >> lengthDivider;
						if (offset + radius < diagHigh) diagHigh = offset + radius;
						if (offset - radius > diagLow) diagLow = offset - radius;
					} else {
						radius += contained_offset;
						if (posIi1 < posJi1) {
							TPos offset = (posJi1 - posIi1);
							radius += (posJi2 - posJi1) >> lengthDivider;
							if (offset + radius < diagHigh) diagHigh = offset + radius;
							if (offset - radius > diagLow) diagLow = offset - radius;
						} else {
							TPos offset = (posIi1 - posJi1);
							radius += (posIi2 - posIi1) >> lengthDivider;
							if (-1 * offset + radius < diagHigh) diagHigh = -1 * offset + radius;
							if (-1 * offset - radius > diagLow) diagLow = -1 * offset - radius;
						}
					}
				} else { // 2) Forward - Reverse
					if ((posJi1 < posIi2) && (posJi2 < posIi1)) {
						TPos offset = (posIi1 - posJi2);
						radius += (posJi1 - posJi2 - offset) >> lengthDivider;
						if (-1 * offset + radius < diagHigh) diagHigh = -1 * offset + radius;
						if (-1 * offset - radius > diagLow) diagLow = -1 * offset - radius;
					} else if ((posJi2 > posIi1) && (posJi1 > posIi2)) {
						TPos offset = (posJi2 - posIi1);
						radius += (posIi2 - posIi1 - offset) >> lengthDivider;
						if (offset + radius < diagHigh) diagHigh = offset + radius;
						if (offset - radius > diagLow) diagLow = offset - radius;
					} else {
						radius += contained_offset;  
						if (posIi1 < posJi2) {
							TPos offset = (posJi2 - posIi1);
							radius += (posJi1 - posJi2) >> lengthDivider;
							if (offset + radius < diagHigh) diagHigh = offset + radius;
							if (offset - radius > diagLow) diagLow = offset - radius;
						} else {
							TPos offset = (posIi1 - posJi2);
							radius += (posIi2 - posIi1) >> lengthDivider;
							if (-1 * offset + radius < diagHigh) diagHigh = -1 * offset + radius;
							if (-1 * offset - radius > diagLow) diagLow = -1 * offset - radius;
						}	
					}
				}
			} else { 
				// 3) Reverse - Forward
				if (forwardJ) {
					if ((posIi1 > posJi2) && (posIi2 > posJi1)) {
						TPos offset = (posIi2 - posJi1);
						radius += (posJi2 - posJi1 - offset) >> lengthDivider;
						if (-1 * offset - radius > diagLow) diagLow = -1 * offset - radius;
						if (-1 * offset + radius < diagHigh) diagHigh = -1 * offset + radius;
					} else if ((posJi1 > posIi2) && (posJi2 > posIi1)) {
						TPos offset = (posJi1 - posIi2);
						radius += (posIi1 - posIi2 - offset) >> lengthDivider;
						if (offset + radius < diagHigh) diagHigh = offset + radius;
						if (offset - radius > diagLow) diagLow = offset - radius;
					} else {
						radius += contained_offset;  
						if (posIi2 < posJi1) {
							TPos offset = (posJi1 - posIi2);
							radius += (posJi2 - posJi1) >> lengthDivider;
							if (offset + radius < diagHigh) diagHigh = offset + radius;
							if (offset - radius > diagLow) diagLow = offset - radius;
						} else {
							TPos offset = (posIi2 - posJi1);
							radius += (posIi1 - posIi2) >> lengthDivider;
							if (-1 * offset + radius < diagHigh) diagHigh = -1 * offset + radius;
							if (-1 * offset - radius > diagLow) diagLow = -1 * offset - radius;
						}
					}
				} else { // 4) Reverse - Reverse
					if ((posJi1 < posIi1) && (posJi2 < posIi2)) {
						TPos offset = (posIi2 - posJi2);
						radius += (posJi1 - posJi2 - offset) >> lengthDivider;
						if (-1 * offset + radius < diagHigh) diagHigh = -1 * offset + radius;
						if (-1 * offset - radius > diagLow) diagLow = -1 * offset - radius;
					} else if ((posJi2 > posIi2) && (posJi1 > posIi1)) {
						TPos offset = (posJi2 - posIi2);
						radius += (posIi1 - posIi2 - offset) >> lengthDivider;
						if (offset + radius < diagHigh) diagHigh = offset + radius;
						if (offset - radius > diagLow) diagLow = offset - radius;
					} else {
						radius += contained_offset;  
						if (posIi2 < posJi2) {
							TPos offset = (posJi2 - posIi2);
							radius += (posJi1 - posJi2) >> lengthDivider;
							if (offset + radius < diagHigh) diagHigh = offset + radius;
							if (offset - radius > diagLow) diagLow = offset - radius;
						} else {
							TPos offset = (posIi2 - posJi2);
							radius += (posIi1 - posIi2) >> lengthDivider;
							if (-1 * offset + radius < diagHigh) diagHigh = -1 * offset + radius;
							if (-1 * offset - radius > diagLow) diagLow = -1 * offset - radius;
						}
					}
				}
			}

			// Append this pair of reads
			if (index1 < index2) {
				appendValue(prePList, TPair(positionToId(str, index1), positionToId(str, index2)));
				appendValue(preDList, TDiagPair(diagLow, diagHigh));
			} else {
				appendValue(prePList, TPair(positionToId(str, index2), positionToId(str, index1)));
				appendValue(preDList, TDiagPair(-1 * diagHigh, -1 * diagLow));
			}
			// Estimate the overlap quality
			TPos avgDiag = (diagLow + diagHigh) / 2;
			if (avgDiag < 0) avgDiag *= -1;
			appendValue(ovlIndex, Pair<TSize, TSize>((TSize) (avgDiag), pairLen)); 
			++pairLen;
		}
	}

	// Sort the pairs, better expected overlaps come first
	std::sort(begin(ovlIndex, Standard() ), end(ovlIndex, Standard() ), _LessPair<TSize>() );
	typedef typename Iterator<TOverlapIndexList>::Type TOVLIter;
	TOVLIter itOvl = begin(ovlIndex);
	TOVLIter itOvlEnd = end(ovlIndex);
	reserve(dList, pairLen);
	reserve(pList, pairLen);
	for(;itOvl != itOvlEnd; goNext(itOvl)) {
		TSize count = (value(itOvl)).i2;
		//// Debug Code
		//std::cout << "Pair: " << count << ',' << "Priority: " << (value(itOvl)).i1 << std::endl;
		//std::cout << "First read: " << value(prePList, count).i1 << ',' <<  "Second read: " << value(prePList, count).i2 << std::endl;
		//std::cout << "Low Diag: " << value(preDList, count).i1 << ',' << "High Diag: " << value(preDList, count).i2 << std::endl;
		//std::cout << value(begEndPos, value(prePList, count).i1).i1 << ',' << value(begEndPos, value(prePList, count).i1).i2 << std::endl;
		//std::cout << value(begEndPos, value(prePList, count).i2).i1 << ',' << value(begEndPos, value(prePList, count).i2).i2 << std::endl;
		appendValue(dList, value(preDList, count));
		appendValue(pList, value(prePList, count));
	}
}



//////////////////////////////////////////////////////////////////////////////

template<typename TString, typename TSpec, typename TBegEndPos, typename TSize, typename TPairList, typename TPos, typename TSpec2>
inline void 
selectPairsIndel(StringSet<TString, TSpec> const& str,
				 TBegEndPos const& begEndPos,
				 TSize lookAround,
				 TPairList& pList,
				 String<Pair<TPos, TPos>, TSpec2>& dList)
{	
	SEQAN_CHECKPOINT	
	typedef String<Pair<TPos, TPos>, TSpec2>  TDistanceList;
	typedef StringSet<TString, TSpec> TStringSet;
	typedef Pair<TPos, TPos> TDiagPair;
	typedef typename Value<TPairList>::Type TPair;
	typedef typename Iterator<TBegEndPos>::Type TBegEndIter;
	
	TBegEndIter beIt = begin(begEndPos);
	TBegEndIter beItEnd = end(begEndPos);
	TSize index1 = 0;
	for(;beIt != beItEnd; ++beIt, ++index1) {
		TPos posIi1 = (value(beIt)).i1;
		TPos posIi2 = (value(beIt)).i2;
		TBegEndIter beIt2 = beIt;
		++beIt2;
		TSize index2 = index1 + 1;
		for(;beIt2 != beItEnd; ++beIt2, ++index2) {
			TPos posJi1 = (value(beIt2)).i1;
			TPos posJi2 = (value(beIt2)).i2;

			TPos beg1 = posJi1;
			TPos diagLow = -1 * (posJi2 - posJi1);
			if (posJi2 < posJi1) {
				diagLow = -1 * (posJi1 - posJi2);
				beg1 = posJi2;
			}
			TPos beg2 = posIi1;
			TPos diagHigh = posIi2 - posIi1;
			if (posIi2 < posIi1) {
				diagHigh = posIi1 - posIi2;
				beg2 = posIi2;
			}

			TPos diff = beg1 - beg2;
			if (diff < 0) diff *= -1;
			if (diff < (TPos) lookAround) {
				appendValue(pList, TPair(positionToId(str, index1), positionToId(str, index2)));
				appendValue(dList, TDiagPair(diagLow, diagHigh));
			}
		}
	}
}


//////////////////////////////////////////////////////////////////////////////

template<typename TString, typename TSpec, typename TId, typename TDiagList, typename TBegEndPos, typename TScore, typename TSize, typename TSegmentMatches, typename TScoreValues, typename TDistance>
inline void 
appendSegmentMatches(StringSet<TString, TSpec> const& str,
					 String<Pair<TId, TId> > const& pList,
					 TDiagList const& dList,
					 TBegEndPos const& begEndPos,
					 TScore const& score_type,
					 TSize thresholdMatchlength,
					 TSize thresholdQuality,
					 TSize maxOvl,
					 TSegmentMatches& matches,
					 TScoreValues& scores,
					 TDistance& dist,
					 Overlap_Library)
{
	SEQAN_CHECKPOINT
	typedef StringSet<TString, Dependent<> > TStringSet;
	typedef String<Pair<TId, TId> > TPairList;
	typedef typename Value<TScoreValues>::Type TScoreValue;
	typedef typename Iterator<TPairList>::Type TPairIter;
	typedef typename Iterator<TDiagList>::Type TDiagIter;

	// Initialization
	double qltThres = (double) thresholdQuality / 100.0;
	TSize nseq = length(str);
	__resizeWithRespectToDistance(dist, nseq);

	// "Front" and "Back"-overlap counter for each read
	String<TSize> frontOvl;
	String<TSize> backOvl;
	fill(frontOvl, nseq, 0);
	fill(backOvl, nseq, 0);
	
	// Pairwise alignments
	String<bool> aligned;
	fill(aligned, length(pList), true);
	typedef Iterator<String<bool> >::Type TBoolIter;
	TBoolIter itAligned = begin(aligned);
	TPairIter itPair = begin(pList);
	TDiagIter itDiag = begin(dList);
	TPairIter itPairEnd = end(pList);
	TSize dropCount = 0;
	for(;itPair != itPairEnd; goNext(itPair), goNext(itDiag), goNext(itAligned)) {
		TId id1 = (value(itPair)).i1;
		TId id2 = (value(itPair)).i2;
		TSize seq1 = idToPosition(str, id1);
		TSize seq2 = idToPosition(str, id2);
		if ((value(frontOvl, seq1) > maxOvl) &&
			(value(backOvl, seq1) > maxOvl) &&
			(value(frontOvl, seq2) > maxOvl) &&
			(value(backOvl, seq2) > maxOvl)) {
				++dropCount;
				continue;
		}


		// Make a pairwise string-set
		TStringSet pairSet;
		assignValueById(pairSet, const_cast<StringSet<TString, TSpec>&>(str), id1);
		assignValueById(pairSet, const_cast<StringSet<TString, TSpec>&>(str), id2);
		
		// Overlap alignment
		TSize from = length(matches);
		//TScoreValue myScore = globalAlignment(matches, pairSet, score_type, AlignConfig<true,true,true,true>(), Gotoh() );
		TScoreValue myScore = globalAlignment(matches, pairSet, score_type, AlignConfig<true,true,true,true>(), (value(itDiag)).i1, (value(itDiag)).i2, BandedGotoh() );
		TSize to = length(matches);

		// Determine a sequence weight
		TSize matchLen = 0;
		TSize overlapLen = 0;
		TSize alignLen = 0;
		getAlignmentStatistics(matches, pairSet, from, matchLen, overlapLen, alignLen);
		double quality = (double) matchLen / (double) overlapLen;

		// Get only the good overlap alignments
		if ((quality >= qltThres) && (matchLen >= thresholdMatchlength)) {

			//// Debug Code
			//Graph<Alignment<TStringSet, TSize> > tmp(pairSet);
			//globalAlignment(tmp, pairSet, score_type, AlignConfig<true,true,true,true>(), (value(itDiag)).i1, (value(itDiag)).i2, BandedGotoh() );
			////globalAlignment(tmp, pairSet, score_type, Gotoh() );
			//std::cout << "Match length: " << matchLen << std::endl;
			//std::cout << "Overlap length: " << overlapLen << std::endl;
			//std::cout << "Align length: " << alignLen << std::endl;
			//std::cout << "Quality: " << quality << std::endl;
			//std::cout << tmp << std::endl;

			// Create a corresponding edge
			TSize i = idToPosition(str, id1);
			TSize j = idToPosition(str, id2);
			if (i<j) __getAlignmentStatistics(dist, i, j, nseq, matchLen, quality);
			else __getAlignmentStatistics(dist, j, i, nseq, matchLen, quality);
			
			// Record the scores
			resize(scores, to);
			typedef typename Iterator<TScoreValues>::Type TScoreIter;
			TScoreIter itScore = begin(scores);
			TScoreIter itScoreEnd = end(scores);
			goFurther(itScore, from);
			for(;itScore != itScoreEnd; ++itScore) value(itScore) = myScore;

			// Update the overlap counter
			TSize lenLast = (value(matches, from)).len; 
			if ((value(matches, to - 1)).begin1 == 0) ++value(frontOvl, seq1);
			if ((value(matches, to - 1)).begin2 == 0) ++value(frontOvl, seq2);
			if ((value(matches, from)).begin1 + lenLast == length(value(pairSet, 0))) ++value(backOvl, seq1);
			if ((value(matches, from)).begin2 + lenLast == length(value(pairSet, 1))) ++value(backOvl, seq2);
		} else {
			resize(matches, from);
			value(itAligned) = false;
		}
	}
	std::cout << "Filtration ration: " << (double) dropCount / (double) length(pList) << std::endl;

	// Find sequences that have no overlap in the front or back
	String<TSize> noFront;
	String<TSize> noBack;
	for(TSize seqI = 0; seqI < nseq; ++seqI) {
		if (value(frontOvl, seqI) == 0) appendValue(noFront, seqI);
		else if (value(backOvl, seqI) == 0) appendValue(noBack, seqI);
	}
	// Drop the first and the last sequence
	typedef typename Iterator<TBegEndPos>::Type TBegEndIter;
	TBegEndIter begEndIt = begin(begEndPos);
	TBegEndIter begEndItEnd = end(begEndPos);
	TSize minVal = supremumValue<TSize>();
	TSize maxVal = 0;
	for(;begEndIt != begEndItEnd; goNext(begEndIt)) {
		TSize pos1 = (value(begEndIt)).i1;
		TSize pos2 = (value(begEndIt)).i2;
		if (pos1 > pos2) { TSize tmp = pos1; pos1 = pos2; pos2 = tmp;}
		if (pos1 < minVal) minVal = pos1;
		if (pos2 > maxVal) maxVal = pos2;
	}
	// Insert all remaining sequences into a set
	std::set<TSize> unalignedReads;
	for(TSize i = 0; i < length(noFront); ++i) {
		TSize p1 = value(begEndPos, value(noFront, i)).i1;
		TSize p2 = value(begEndPos, value(noFront, i)).i2;
		if (p1 > p2) {TSize tmp = p1; p1 = p2; p2 = tmp; }
		if (p1 != minVal) unalignedReads.insert(value(noFront, i));
	}
	for(TSize i = 0; i < length(noBack); ++i) {
		TSize p1 = value(begEndPos, value(noBack, i)).i1;
		TSize p2 = value(begEndPos, value(noBack, i)).i2;
		if (p1 > p2) {TSize tmp = p1; p1 = p2; p2 = tmp; }
		if (p2 != maxVal) unalignedReads.insert(value(noBack, i));
	}
	TSize countUnalignedReads = unalignedReads.size();
	std::cout << "Unaligned reads: " << countUnalignedReads << std::endl;
	if (countUnalignedReads > 0) {
		// Realign all unaligned sequences
		itPair = begin(pList);
		itDiag = begin(dList);
		itAligned = begin(aligned);
		for(;itPair != itPairEnd; goNext(itPair), goNext(itDiag), goNext(itAligned)) {
			if (value(itAligned) == true) continue;
			TId id1 = (value(itPair)).i1;
			TId id2 = (value(itPair)).i2;
			TSize seq1 = idToPosition(str, id1);
			TSize seq2 = idToPosition(str, id2);
			if ((unalignedReads.find(seq1) == unalignedReads.end()) &&
				(unalignedReads.find(seq2) == unalignedReads.end())) continue;
			if ((value(frontOvl, seq1) > maxOvl) &&
				(value(backOvl, seq1) > maxOvl) &&
				(value(frontOvl, seq2) > maxOvl) &&
				(value(backOvl, seq2) > maxOvl)) {
					continue;
			}

			// Make a pairwise string-set
			TStringSet pairSet;
			assignValueById(pairSet, const_cast<StringSet<TString, TSpec>&>(str), id1);
			assignValueById(pairSet, const_cast<StringSet<TString, TSpec>&>(str), id2);
		
			// Overlap alignment
			TSize from = length(matches);
#ifdef CELERA_OFFSET
			TScoreValue myScore = globalAlignment(matches, pairSet, score_type, AlignConfig<true,true,true,true>(), Gotoh() );
#else
			TScoreValue myScore = globalAlignment(matches, pairSet, score_type, AlignConfig<true,true,true,true>(), (value(itDiag)).i1, (value(itDiag)).i2, BandedGotoh() );
#endif
			TSize to = length(matches);

			// Determine a sequence weight
			TSize matchLen = 0;
			TSize overlapLen = 0;
			TSize alignLen = 0;
			getAlignmentStatistics(matches, pairSet, from, matchLen, overlapLen, alignLen);
			double quality = (double) matchLen / (double) overlapLen;

			if ((quality >= 0.8) && (matchLen >= 5)) {
				// Create a corresponding edge
				TSize i = idToPosition(str, id1);
				TSize j = idToPosition(str, id2);
				if (i<j) __getAlignmentStatistics(dist, i, j, nseq, matchLen, quality);
				else __getAlignmentStatistics(dist, j, i, nseq, matchLen, quality);

				// Record the scores
				resize(scores, to);
				typedef typename Iterator<TScoreValues>::Type TScoreIter;
				TScoreIter itScore = begin(scores);
				TScoreIter itScoreEnd = end(scores);
				goFurther(itScore, from);
				for(;itScore != itScoreEnd; ++itScore) value(itScore) = myScore;

				// Update the overlap counter
				TSize lenLast = (value(matches, from)).len; 
				if ((value(matches, to - 1)).begin1 == 0) ++value(frontOvl, seq1);
				if ((value(matches, to - 1)).begin2 == 0) ++value(frontOvl, seq2);
				if ((value(matches, from)).begin1 + lenLast == length(value(pairSet, 0))) ++value(backOvl, seq1);
				if ((value(matches, from)).begin2 + lenLast == length(value(pairSet, 1))) ++value(backOvl, seq2);
			} else {
				resize(matches, from);
			}
		}
	}
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...

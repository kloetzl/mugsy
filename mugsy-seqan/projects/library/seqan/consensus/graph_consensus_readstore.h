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
  $Id: graph_consensus_readstore.h 1811 2008-03-31 15:38:54Z rausch@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_GRAPH_CONSENSUS_READSTORE_H
#define SEQAN_HEADER_GRAPH_CONSENSUS_READSTORE_H

namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////
// Read Store
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet = Dna, typename TSpec = Default>
class ReadStore;

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TSpec>
class ReadStore
{
	public:
		typedef typename Size<ReadStore>::Type TSize;
		String<TAlphabet, External<> > data_reads;
		//String<TAlphabet, Alloc<> > data_reads;
		String<char, External<> > data_qualities;
		//String<char, Alloc<> > data_qualities;
		StringSet<String<char, External<> >, Owner<ConcatDirect<> > > data_names;
		//StringSet<String<char, Alloc<> > > data_names;
		String<Pair<TSize, TSize> > data_begin_end;
		String<Pair<TSize, TSize> > data_clr;
		String<TSize> data_frg_id;
		TSize data_pos_count;
		
	public:
		ReadStore() : data_pos_count(0)
		{
			SEQAN_CHECKPOINT
			clear(data_reads);
			clear(data_qualities);
			clear(data_names);
			clear(data_begin_end);
			clear(data_clr);
			clear(data_frg_id);
		}

		~ReadStore() 
		{
			SEQAN_CHECKPOINT
		}


	private:
		ReadStore(ReadStore const & _other)
		{
			SEQAN_CHECKPOINT
			data_reads = _other.data_reads;
			data_qualities = _other.data_qualities;
			data_names = _other.data_names;
			data_begin_end = _other.data_begin_end;
			data_frg_id = _other.data_frg_id;
			data_clr = _other.data_clr;
			data_pos_count = _other.data_pos_count;
		}

		ReadStore const& 
		operator = (ReadStore const& _other) 
		{
			SEQAN_CHECKPOINT
			if (this == &_other) return *this;
			data_reads = _other.data_reads;
			data_qualities = _other.data_qualities;
			data_names = _other.data_names;
			data_begin_end = _other.data_begin_end;
			data_frg_id = _other.data_frg_id;
			data_clr = _other.data_clr;
			data_pos_count = _other.data_pos_count;
			return *this;
		}
};
	


//////////////////////////////////////////////////////////////////////////////
// FUNCTIONS
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TSpec, typename TSize, typename TRead>
inline void 
loadRead(ReadStore<TAlphabet, TSpec>& readSt, 
		 TSize index,
		 TRead& seq) 
{
	SEQAN_CHECKPOINT
	clear(seq);
	seq = infix(readSt.data_reads, (value(readSt.data_begin_end, index)).i1, (value(readSt.data_begin_end, index)).i2);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TSpec, typename TSpec2, typename TSize, typename TRead>
inline void 
loadRead(ReadStore<TAlphabet, TSpec>& readSt,
		 String<TAlphabet, TSpec2>& readString, 
		 TSize index,
		 TRead& seq) 
{
	SEQAN_CHECKPOINT
	clear(seq);
	seq = infix(readString, (value(readSt.data_begin_end, index)).i1, (value(readSt.data_begin_end, index)).i2);
}




//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TSpec, typename TSize, typename TRead>
inline void 
loadQuality(ReadStore<TAlphabet, TSpec>& readSt,
			TSize index,
			TRead& seq) 
{
	SEQAN_CHECKPOINT
	clear(seq);
	seq = infix(readSt.data_qualities, (value(readSt.data_begin_end, index)).i1, (value(readSt.data_begin_end, index)).i2);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TSpec, typename TSpec2, typename TSize, typename TRead>
inline void 
loadQuality(ReadStore<TAlphabet, TSpec>& readSt,
			String<char, TSpec2>& qualityString, 
			TSize index,
			TRead& seq) 
{
	SEQAN_CHECKPOINT
	clear(seq);
	seq = infix(qualityString, (value(readSt.data_begin_end, index)).i1, (value(readSt.data_begin_end, index)).i2);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TSpec, typename TSize, typename TName>
inline void 
loadExtId(ReadStore<TAlphabet, TSpec>& readSt,
		  TSize index,
		  TName& name) 
{
	SEQAN_CHECKPOINT
	name = value(readSt.data_names, index);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TSpec, typename TSize>
inline typename Id<ReadStore<TAlphabet, TSpec> >::Type
loadFrgId(ReadStore<TAlphabet, TSpec>& readSt, 
		  TSize index) 
{
	SEQAN_CHECKPOINT
	return value(readSt.data_frg_id, index);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TSpec>
inline typename Size<ReadStore<TAlphabet, TSpec> >::Type
length(ReadStore<TAlphabet, TSpec>& readSt) 
{
	SEQAN_CHECKPOINT
	return readSt.data_pos_count;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TSpec, typename TAlph2, typename TSpec2, typename TSize, typename TStringSet, typename TLayoutPos>
inline void
loadReadsClr(ReadStore<TAlphabet, TSpec>& readSt,
			 CtgStore<TAlph2, TSpec2>& ctgSt,
			 TSize index,
			 TStringSet& strSet,
			 TLayoutPos& startEndPos) 
{
	SEQAN_CHECKPOINT
	typedef typename Value<TLayoutPos>::Type TPair;
	String<GappedRead<> >& gapReads = value(ctgSt.data_reads, index);
	TSize numReads = length(gapReads);
	resize(strSet, numReads);
	String<TAlphabet> all = readSt.data_reads;
	for(TSize i = 0; i<numReads; ++i) {
		GappedRead<>& gRead = value(gapReads, i);
		String<TAlphabet> seq;
		loadRead(readSt, all, gRead.data_source , seq);
		if (gRead.data_clr.i1 < gRead.data_clr.i2) {
			value(strSet, i) = infix(seq, gRead.data_clr.i1, gRead.data_clr.i2);
		} else {
			value(strSet, i) = infix(seq, gRead.data_clr.i2, gRead.data_clr.i1);
			reverseComplementInPlace(value(strSet, i));
		}
		appendValue(startEndPos, TPair(gRead.data_clr.i1 + gRead.data_offset, gRead.data_clr.i2 + gRead.data_offset));
	}
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TSpec, typename TAlph2, typename TSpec2, typename TSize, typename TReadClrSet, typename TQualityClrSet, typename TLayoutPos>
inline void
loadReadsClr(ReadStore<TAlphabet, TSpec>& readSt,
             CtgStore<TAlph2, TSpec2>& ctgSt,
             TSize index,
             TReadClrSet& readSet,
             TQualityClrSet& qualitySet,
             TLayoutPos& startEndPos) 
{
    SEQAN_CHECKPOINT

    typedef typename Value<TLayoutPos>::Type TPair;
    String<GappedRead<> >& gapReads = value(ctgSt.data_reads, index);
    TSize numReads = length(gapReads);
	
    resize(readSet, numReads);
    resize(qualitySet, numReads);

	String<TAlphabet> all = readSt.data_reads;
    for(TSize i = 0; i<numReads; ++i) {
        GappedRead<>& gRead = value(gapReads, i);
        String<TAlphabet> seq;
        String<char> quality;
        loadRead(readSt, all, gRead.data_source , seq);
        loadQuality(readSt,gRead.data_source, quality);
        if (gRead.data_clr.i1 < gRead.data_clr.i2) {
            value(readSet, i) = infix(seq, gRead.data_clr.i1, gRead.data_clr.i2);
            value(qualitySet, i) = infix(quality, gRead.data_clr.i1, gRead.data_clr.i2);
        } else {
            value(readSet, i) = infix(seq, gRead.data_clr.i2, gRead.data_clr.i1);
            value(qualitySet, i) = infix(quality, gRead.data_clr.i2, gRead.data_clr.i1);
            reverseComplementInPlace(value(readSet, i));
            reverseInPlace(value(qualitySet,i) );
        }
        appendValue(startEndPos, TPair(gRead.data_clr.i1 + gRead.data_offset, gRead.data_clr.i2 + gRead.data_offset));
    }
}

//////////////////////////////////////////////////////////////////////////////


template<typename TFile, typename TAlphabet, typename TSpec, typename TFragmentStore, typename TLibraryStore, typename TContigStore>
inline void 
read(TFile & file,
	 ReadStore<TAlphabet, TSpec>& readSt,
	 TFragmentStore& frgSt,
	 TLibraryStore& libSt,
	 TContigStore& ctgSt,
	 Amos) 
{
	SEQAN_CHECKPOINT

	typedef ReadStore<TAlphabet, TSpec> TReadStore;
	typedef typename Id<TReadStore>::Type TId;
	typedef typename Size<TReadStore>::Type TSize;
	typedef typename Value<TFile>::Type TValue;

	// All maps to mirror file ids to our ids
	typedef std::map<TSize, TSize> TIdMap;
	TIdMap libIdMap;
	TIdMap frgIdMap;
	TIdMap readIdMap;
	TIdMap ctgIdMap;

	// Parse the file and convert the internal ids
	TValue c;
	if (_streamEOF(file)) return;
	else c = _streamGet(file);
	while (!_streamEOF(file)) {
		if (_streamEOF(file)) break;

		// New block?
		if (c == '{') {
			c = _streamGet(file);
			String<char> blockIdentifier;
			_parse_readIdentifier(file, blockIdentifier, c);
			_parse_skipLine(file, c);

			// Library block
			if (blockIdentifier == "LIB") {
				TSize id = 0;
				TSize mean = 0;
				TSize std = 0;
				String<char> fieldIdentifier;
				String<char> eid;
				while (c != '}') {
					clear(fieldIdentifier);
					_parse_readIdentifier(file, fieldIdentifier, c);
					if (fieldIdentifier == "iid") {
						c = _streamGet(file);
						id = _parse_readNumber(file, c);
						_parse_skipLine(file, c);
					} else if (fieldIdentifier == "eid") {
						c = _streamGet(file);
						while ((c != '\n') && (c != '\r')) {
							appendValue(eid, c);
							c = _streamGet(file);
						}
						_parse_skipLine(file, c);
					} else if (fieldIdentifier == "mea") {
						c = _streamGet(file);
						mean = _parse_readNumber(file, c);
						_parse_skipLine(file, c);
					} else if (fieldIdentifier == "std") {
						c = _streamGet(file);
						std = _parse_readNumber(file, c);
						_parse_skipLine(file, c);
					} else {
						_parse_skipLine(file, c);
					}
				}
				appendValue(libSt.data_mean, mean);
				appendValue(libSt.data_std, std);
				appendValue(libSt.data_names, eid);
				libIdMap.insert(std::make_pair(id, libSt.data_pos_count));
				++libSt.data_pos_count;
			} else if (blockIdentifier == "FRG") {  // Fragment block
				TSize libId = 0;
				TSize id = 0;
				TSize rds1Id = 0;
				TSize rds2Id = 0;
				String<char> fieldIdentifier;
				String<char> eid;
				while (c != '}') {
					clear(fieldIdentifier);
					_parse_readIdentifier(file, fieldIdentifier, c);
					if (fieldIdentifier == "iid") {
						c = _streamGet(file);
						id = _parse_readNumber(file, c);
						_parse_skipLine(file, c);
					} else if (fieldIdentifier == "eid") {
						c = _streamGet(file);
						while ((c != '\n') && (c != '\r')) {
							appendValue(eid, c);
							c = _streamGet(file);
						}
						_parse_skipLine(file, c);
					} else if (fieldIdentifier == "lib") {
						c = _streamGet(file);
						libId = _parse_readNumber(file, c);
						_parse_skipLine(file, c);
					} else if (fieldIdentifier == "rds") {
						c = _streamGet(file);
						rds1Id = _parse_readNumber(file, c);
						c = _streamGet(file);
						rds2Id = _parse_readNumber(file, c);
						_parse_skipLine(file, c);
					} else {
						_parse_skipLine(file, c);
					}
				}
				appendValue(frgSt.data_lib_id, libId);
				appendValue(frgSt.data_rds, Pair<TSize, TSize>(rds1Id, rds2Id));
				appendValue(frgSt.data_names, eid);
				frgIdMap.insert(std::make_pair(id, frgSt.data_pos_count));
				++frgSt.data_pos_count;
			} else if (blockIdentifier == "RED") {   // Read block
				TId id = 0;
				TId frgId = 0;
				TSize clr1 = 0;
				TSize clr2 = 0;
				String<char> fieldIdentifier;
				String<char> eid;
				TSize begRead = 0;
				while (c != '}') {
					clear(fieldIdentifier);
					_parse_readIdentifier(file, fieldIdentifier, c);
					if (fieldIdentifier == "iid") {
						c = _streamGet(file);
						id = _parse_readNumber(file, c);
						_parse_skipLine(file, c);
					} else if (fieldIdentifier == "eid") {
						c = _streamGet(file);
						while ((c != '\n') && (c != '\r')) {
							appendValue(eid, c);
							c = _streamGet(file);
						}
						_parse_skipLine(file, c);
					} else if (fieldIdentifier == "frg") {
						c = _streamGet(file);
						frgId = _parse_readNumber(file, c);
						_parse_skipLine(file, c);
					} else if (fieldIdentifier == "seq") {
						if (readSt.data_pos_count == 0) begRead = 0;
						else begRead = (value(readSt.data_begin_end, readSt.data_pos_count - 1)).i2;
						c = _streamGet(file);
						_parse_skipWhitespace(file, c);
						while (c != '.') {
							_parse_readSequenceData(file,c,readSt.data_reads);
							_parse_skipLine(file, c);
						}
					} else if (fieldIdentifier == "qlt") {
						c = _streamGet(file);
						_parse_skipWhitespace(file, c);
						while (c != '.') {
							if ((c!=' ') && (c != '\t') && (c != '\n') && (c != '\r')) {
								appendValue(readSt.data_qualities, c);
							}
							c = _streamGet(file);
						}
					} else if (fieldIdentifier == "clr") {
						c = _streamGet(file);
						clr1 = _parse_readNumber(file, c);
						c = _streamGet(file);
						clr2 = _parse_readNumber(file, c);
						_parse_skipLine(file, c);
					} else {
						_parse_skipLine(file, c);
					}
				}
				TSize endRead = length(readSt.data_reads);
				appendValue(readSt.data_begin_end, Pair<TSize,TSize>(begRead, endRead));
				appendValue(readSt.data_frg_id, frgId);
				appendValue(readSt.data_names, eid);
				appendValue(readSt.data_clr, Pair<TSize, TSize>(clr1, clr2));
				readIdMap.insert(std::make_pair(id, readSt.data_pos_count));
				++readSt.data_pos_count;
			} else if (blockIdentifier == "CTG") {   // Contig block
				TId id = 0;
				String<char> fieldIdentifier;
				String<char> eid;
				appendValue(ctgSt.data_reads, String<GappedRead<> >());
				TSize begContig = 0;
				while (c != '}') {
					// Are we entering a TLE block
					if (c == '{') {
						GappedRead<> gapRead;
						gapRead.data_gap = String<TSize>();
						String<char> fdIdentifier;
						while (c != '}') {
							clear(fdIdentifier);
							_parse_readIdentifier(file, fdIdentifier, c);
							if (fdIdentifier == "src") {
								c = _streamGet(file);
								gapRead.data_source = _parse_readNumber(file, c);
								_parse_skipLine(file, c);
							} else if (fdIdentifier == "off") {
								c = _streamGet(file);
								gapRead.data_offset = _parse_readNumber(file, c);
								_parse_skipLine(file, c);
							} else if (fdIdentifier == "clr") {
								c = _streamGet(file);
								TSize clr1 = _parse_readNumber(file, c);
								c = _streamGet(file);
								TSize clr2 = _parse_readNumber(file, c);
								_parse_skipLine(file, c);
								gapRead.data_clr = Pair<TSize, TSize>(clr1, clr2);
							} else if (fdIdentifier == "gap") {
								c = _streamGet(file);
								_parse_skipWhitespace(file, c);
								while (c != '.') {
									if ((c!=' ') && (c != '\t') && (c != '\n') && (c != '\r')) {
										TSize nextGap = _parse_readNumber(file, c);
										appendValue(gapRead.data_gap, nextGap);
									}
									c = _streamGet(file);
								}
							} else {
								_parse_skipLine(file, c);
							}
						}
						_parse_skipLine(file, c);
						appendValue(value(ctgSt.data_reads, ctgSt.data_pos_count), gapRead);
					} else {
						clear(fieldIdentifier);
						_parse_readIdentifier(file, fieldIdentifier, c);
						if (fieldIdentifier == "iid") {
							c = _streamGet(file);
							id = _parse_readNumber(file, c);
							_parse_skipLine(file, c);
						} else if (fieldIdentifier == "eid") {
							c = _streamGet(file);
							while ((c != '\n') && (c != '\r')) {
								appendValue(eid, c);
								c = _streamGet(file);
							}
							_parse_skipLine(file, c);
						} else if (fieldIdentifier == "seq") {
							if (ctgSt.data_pos_count == 0) begContig = 0;
							else begContig = (value(ctgSt.data_begin_end, ctgSt.data_pos_count - 1)).i2;
							c = _streamGet(file);
							_parse_skipWhitespace(file, c);
							while (c != '.') {
								// Also include the gaps
								do {
									_parse_readSequenceData(file,c,ctgSt.data_contig);
								} while (c == '-');
								_parse_skipLine(file, c);
							}
						} else if (fieldIdentifier == "qlt") {
							c = _streamGet(file);
							_parse_skipWhitespace(file, c);
							while (c != '.') {
								if ((c!=' ') && (c != '\t') && (c != '\n') && (c != '\r')) {
									appendValue(ctgSt.data_quality, c);
								}
								c = _streamGet(file);
							}
						} else {
							_parse_skipLine(file, c);
						}
					}
				}			
				TSize endContig = length(ctgSt.data_contig);
				appendValue(ctgSt.data_begin_end, Pair<TSize,TSize>(begContig, endContig));
				appendValue(ctgSt.data_names, eid);
				ctgIdMap.insert(std::make_pair(id, ctgSt.data_pos_count));
				++ctgSt.data_pos_count;
			}
		} else {
			_parse_skipLine(file, c);
		}
	}

	// Renumber all ids
	for(TSize i = 0; i<length(frgSt); ++i) {
		value(frgSt.data_lib_id, i) = (libIdMap.find(value(frgSt.data_lib_id, i)))->second;
		TId id1 = (value(frgSt.data_rds, i)).i1;
		TId id2 = (value(frgSt.data_rds, i)).i2;
		if (id1 != id2) {
			(value(frgSt.data_rds, i)).i1 = (readIdMap.find(id1))->second;
			(value(frgSt.data_rds, i)).i2 = (readIdMap.find(id2))->second;
		}
	}
	for(TSize i = 0; i<length(readSt); ++i) {
		value(readSt.data_frg_id, i) = (frgIdMap.find(value(readSt.data_frg_id, i)))->second;
	}
	for(TSize i = 0; i<length(ctgSt); ++i) {
		for(TSize k = 0; k<length(value(ctgSt.data_reads, i)); ++k) {
			(value(value(ctgSt.data_reads, i), k)).data_source = (readIdMap.find((value(value(ctgSt.data_reads, i), k)).data_source))->second;
		}
	}
}


//////////////////////////////////////////////////////////////////////////////


template<typename TFile, typename TAlphabet, typename TSpec, typename TFragmentStore, typename TLibraryStore, typename TContigStore>
inline void 
write(TFile& target,
	  ReadStore<TAlphabet, TSpec>& readSt,
	  TFragmentStore& frgSt,
	  TLibraryStore& libSt,
	  TContigStore& ctgSt,
	  Amos) 
{
	SEQAN_CHECKPOINT
	typedef ReadStore<TAlphabet, TSpec> TReadStore;
	typedef typename Id<TReadStore>::Type TId;
	typedef typename Size<TReadStore>::Type TSize;
	typedef typename Value<TFile>::Type TValue;
	typedef String<char> TCharString;

	// Write Header
	_streamWrite(target,"{UNV\niid:1\ncom:\nafg file created with SeqAn\n.\n}\n");
	
	// Write Libraries
	for(TSize i = 0; i<length(libSt); ++i) {
		_streamWrite(target,"{LIB\n");
		_streamWrite(target,"iid:");
		_streamPutInt(target, i + 1);
		_streamPut(target, '\n');
		_streamWrite(target,"eid:");
		_streamWrite(target, value(libSt.data_names, i));
		_streamPut(target, '\n');
		_streamWrite(target,"{DST\n");
		_streamWrite(target,"mea:");
		_streamPutInt(target, value(libSt.data_mean, i));
		_streamPut(target, '\n');
		_streamWrite(target,"std:");
		_streamPutInt(target, value(libSt.data_std, i));
		_streamPut(target, '\n');
		_streamWrite(target,"}\n");	
		_streamWrite(target,"}\n");
	}

	// Write Fragments
	for(TSize i = 0; i<length(frgSt); ++i) {
		_streamWrite(target,"{FRG\n");
		_streamWrite(target,"iid:");
		_streamPutInt(target, i + 1);
		_streamPut(target, '\n');
		_streamWrite(target,"eid:");
		_streamWrite(target, value(frgSt.data_names, i));
		_streamPut(target, '\n');
		_streamWrite(target,"lib:");
		_streamPutInt(target, (value(frgSt.data_lib_id, i) + 1));
		_streamPut(target, '\n');
		TId id1 = (value(frgSt.data_rds, i)).i1;
		TId id2 = (value(frgSt.data_rds, i)).i2;
		if (id1 != id2) {	// Just print if this is a proper mate pair
			_streamWrite(target,"rds:");
			_streamPutInt(target, id1 + 1);
			_streamPut(target, ',');
			_streamPutInt(target, id2 + 1);
			_streamPut(target, '\n');
		}
		_streamWrite(target,"}\n");
	}

	// Write Reads
	String<TAlphabet> all = readSt.data_reads;
	String<char> allQuality = readSt.data_qualities;
	for(TSize i = 0; i<length(readSt); ++i) {
		_streamWrite(target,"{RED\n");
		_streamWrite(target,"iid:");
		_streamPutInt(target, i + 1);
		_streamPut(target, '\n');
		_streamWrite(target,"eid:");
		_streamWrite(target, value(readSt.data_names, i));
		_streamPut(target, '\n');
		_streamWrite(target,"seq:\n");
		String<TAlphabet> read;
		loadRead(readSt, all, i, read);
		for(TSize k = 0;k<length(read); k+=60) {
			TSize endK = k + 60;
			if (endK > length(read)) endK = length(read);
			_streamWrite(target, infix(read, k, endK));
			_streamPut(target, '\n');
		}
		_streamWrite(target, ".\n");
		_streamWrite(target,"qlt:\n");
		String<char> qlt;
		loadQuality(readSt, allQuality, i, qlt);
		for(TSize k = 0;k<length(qlt); k+=60) {
			TSize endK = k + 60;
			if (endK > length(qlt)) endK = length(qlt);
			_streamWrite(target, infix(qlt, k, endK));
			_streamPut(target, '\n');
		}
		_streamWrite(target, ".\n");
		_streamWrite(target,"frg:");
		_streamPutInt(target, value(readSt.data_frg_id, i) + 1);
		_streamPut(target, '\n');
		_streamWrite(target,"clr:");
		_streamPutInt(target, value(readSt.data_clr, i).i1);
		_streamPut(target, ',');
		_streamPutInt(target, value(readSt.data_clr, i).i2);
		_streamPut(target, '\n');
		_streamWrite(target,"}\n");
	}
	clear(all);
	clear(allQuality);

	// Write Contigs
	for(TSize i = 0; i<length(ctgSt); ++i) {
		_streamWrite(target,"{CTG\n");
		_streamWrite(target,"iid:");
		_streamPutInt(target, i + 1);
		_streamPut(target, '\n');
		_streamWrite(target,"eid:");
		_streamWrite(target, value(ctgSt.data_names, i));
		_streamPut(target, '\n');
		_streamWrite(target,"seq:\n");
		String<char> contig;
		loadContig(ctgSt, i, contig);
		for(TSize k = 0;k<length(contig); k+=70) {
			TSize endK = k + 70;
			if (endK > length(contig)) endK = length(contig);
			_streamWrite(target, infix(contig, k, endK));
			_streamPut(target, '\n');
		}
		_streamWrite(target, ".\n");
		_streamWrite(target,"qlt:\n");
		String<char> qlt;
		loadQuality(ctgSt, i, qlt);
		for(TSize k = 0;k<length(qlt); k+=70) {
			TSize endK = k + 70;
			if (endK > length(qlt)) endK = length(qlt);
			_streamWrite(target, infix(qlt, k, endK));
			_streamPut(target, '\n');
		}
		_streamWrite(target, ".\n");
		String<GappedRead<> >& gappedReads = value(ctgSt.data_reads, i);
		for(TSize k = 0; k<length(gappedReads); ++k) {
			_streamWrite(target,"{TLE\n");
			_streamWrite(target,"src:");
			_streamPutInt(target, (value(gappedReads, k)).data_source + 1);
			_streamPut(target, '\n');
			_streamWrite(target,"off:");
			_streamPutInt(target, (value(gappedReads, k)).data_offset);
			_streamPut(target, '\n');
			TSize clr1 = (value(gappedReads, k)).data_clr.i1;
			TSize clr2 = (value(gappedReads, k)).data_clr.i2;
			_streamWrite(target,"clr:");
			_streamPutInt(target, clr1);
			_streamPut(target, ',');
			_streamPutInt(target, clr2);
			_streamPut(target, '\n');
			String<TSize>& gaps = (value(gappedReads, k)).data_gap;
			if (!empty(gaps)) {
				_streamWrite(target,"gap:\n");
				for(TSize z = 0;z<length(gaps); ++z) {
					_streamPutInt(target, value(gaps, z));
					_streamPut(target, ' ');
				}
				_streamPut(target, '\n');
				_streamWrite(target, ".\n");
			}
			_streamWrite(target,"}\n");
		}
		_streamWrite(target,"}\n");
	}
}



//////////////////////////////////////////////////////////////////////////////


template<typename TFile, typename TAlphabet, typename TSpec, typename TFragmentStore, typename TLibraryStore, typename TContigStore>
inline void 
write(TFile& target,
	  ReadStore<TAlphabet, TSpec>& readSt,
	  TFragmentStore&,
	  TLibraryStore&,
	  TContigStore&,
	  CeleraFrg) 
{
	SEQAN_CHECKPOINT
	typedef ReadStore<TAlphabet, TSpec> TReadStore;
	typedef typename Id<TReadStore>::Type TId;
	typedef typename Size<TReadStore>::Type TSize;
	typedef typename Value<TFile>::Type TValue;
	typedef String<char> TCharString;

	// Write Reads
	for(TSize i = 0; i<length(readSt); ++i) {
		_streamWrite(target,"{FRG\n");
		_streamWrite(target,"act:");
		_streamPut(target, 'A');
		_streamPut(target, '\n');
		_streamWrite(target,"acc:");
		_streamPutInt(target, i + 1);
		_streamPut(target, '\n');
		_streamWrite(target,"typ:");
		_streamPut(target, 'R');
		_streamPut(target, '\n');
		_streamWrite(target,"src:\n");
		_streamWrite(target, value(readSt.data_names, i));
		_streamWrite(target, "\n.\n");
		_streamWrite(target,"etm:");
		_streamPut(target, '0');
		_streamPut(target, '\n');
		_streamWrite(target,"seq:\n");
		String<TAlphabet> read;
		loadRead(readSt, i, read);
		for(TSize k = 0;k<length(read); k+=70) {
			TSize endK = k + 70;
			if (endK > length(read)) endK = length(read);
			_streamWrite(target, infix(read, k, endK));
			_streamPut(target, '\n');
		}
		_streamWrite(target, ".\n");
		_streamWrite(target,"qlt:\n");
		String<char> qlt;
		loadQuality(readSt, i, qlt);
		for(TSize k = 0;k<length(qlt); k+=70) {
			TSize endK = k + 70;
			if (endK > length(qlt)) endK = length(qlt);
			_streamWrite(target, infix(qlt, k, endK));
			_streamPut(target, '\n');
		}
		_streamWrite(target, ".\n");
		_streamWrite(target,"clr:");
		_streamPutInt(target, value(readSt.data_clr, i).i1);
		_streamPut(target, ',');
		_streamPutInt(target, value(readSt.data_clr, i).i2);
		_streamPut(target, '\n');
		_streamWrite(target,"}\n");
	}
}


//////////////////////////////////////////////////////////////////////////////


template<typename TFile, typename TAlphabet, typename TSpec, typename TFragmentStore, typename TLibraryStore, typename TContigStore>
inline void 
write(TFile& target,
	  ReadStore<TAlphabet, TSpec>& readSt,
	  TFragmentStore&,
	  TLibraryStore&,
	  TContigStore& ctgSt,
	  CeleraCgb) 
{
	SEQAN_CHECKPOINT
	typedef ReadStore<TAlphabet, TSpec> TReadStore;
	typedef typename Id<TReadStore>::Type TId;
	typedef typename Size<TReadStore>::Type TSize;
	typedef typename Value<TFile>::Type TValue;
	typedef String<char> TCharString;
	String<GappedRead<> >& origGappedReads = value(ctgSt.data_reads, 0);
	String<char> contig;
	loadContig(ctgSt, 0, contig);

	// Find smallest offset
	TSize offsetLeft = length(contig);
	for(TSize k = 0; k<length(origGappedReads); ++k) if ((value(origGappedReads, k)).data_offset < offsetLeft) offsetLeft = (value(origGappedReads, k)).data_offset;

	// Sort the reads
	typedef std::set<std::pair<TSize, TSize> > TOffsetIndexPair;
	TOffsetIndexPair offsetIndexSet;
	for(TSize k = 0; k<length(origGappedReads); ++k) {
		TSize clr1 = (value(origGappedReads, k)).data_clr.i1;
		TSize clr2 = (value(origGappedReads, k)).data_clr.i2;
		clr1 += ((value(origGappedReads, k)).data_offset - offsetLeft);
		clr2 += ((value(origGappedReads, k)).data_offset - offsetLeft);
		if (clr1 > clr2) { TSize tmp = clr1; clr1 = clr2; clr2 = tmp; }
		offsetIndexSet.insert(std::make_pair(clr1, k));
	}
	String<GappedRead<> > gappedReads;
	for(typename TOffsetIndexPair::const_iterator itPos = offsetIndexSet.begin(); itPos != offsetIndexSet.end(); ++itPos) {
		appendValue(gappedReads, value(origGappedReads, (*itPos).second));
	}

	//// Write IAF record for all reads
	//for(TSize k = 0; k<length(gappedReads); ++k) {
	//	_streamWrite(target,"{IAF\n");
	//	_streamWrite(target,"acc:");
	//	_streamPutInt(target, (value(gappedReads, k)).data_source + 1);
	//	_streamPut(target, '\n');
	//	_streamWrite(target,"typ:");
	//	_streamPut(target, 'R');
	//	_streamPut(target, '\n');
	//	_streamWrite(target,"chi:0\ncha:0\nclr:-1,-1\nmst:U\n}\n");
	//}

	// Write Header
	_streamWrite(target,"{IUM\nacc:0\nsrc:\ngen> @@ [0,0]\n.\ncov:0.000\nsta:X\nfur:X\nabp:0\nbbp:0\n");
	_streamWrite(target,"len:");
	_streamPutInt(target, length(contig));
	_streamPut(target, '\n');
	_streamWrite(target,"cns:\n.\nqlt:\n.\nfor:0\n");
	_streamWrite(target,"nfr:");
	_streamPutInt(target, length(readSt));
	_streamPut(target, '\n');
	
	// Write gapped reads
	for(TSize k = 0; k<length(gappedReads); ++k) {
		_streamWrite(target,"{IMP\n");
		_streamWrite(target,"typ:");
		_streamPut(target, 'R');
		_streamPut(target, '\n');
		_streamWrite(target,"mid:");
		_streamPutInt(target, (value(gappedReads, k)).data_source + 1);
		_streamPut(target, '\n');
		TSize clr1 = (value(gappedReads, k)).data_clr.i1;
		TSize clr2 = (value(gappedReads, k)).data_clr.i2;
		clr1 += ((value(gappedReads, k)).data_offset - offsetLeft);
		clr2 += ((value(gappedReads, k)).data_offset - offsetLeft);
		_streamWrite(target,"con:");
		//TSize orig1 = clr1;
		//TSize orig2 = clr2;
		//if (orig1 > orig2) {TSize tmp = orig1; orig1 = orig2; orig2 = tmp; }
		//TSize best = 0;
		//TSize bestDist = 0;
		//for(TSize other = 0; other<length(gappedReads); ++other) {
		//	if (other == k) continue;
		//	TSize thisRead1 = (value(gappedReads, other)).data_clr.i1;
		//	TSize thisRead2 = (value(gappedReads, other)).data_clr.i2;
		//	thisRead1 += ((value(gappedReads, other)).data_offset - offsetLeft);
		//	thisRead2 += ((value(gappedReads, other)).data_offset - offsetLeft);
		//	if (thisRead1 > thisRead2) {TSize tmp = thisRead1; thisRead1 = thisRead2; thisRead2 = tmp; }
		//	if ((orig1 > thisRead1) && (orig2 < thisRead2)) {
		//		if ((best == 0) ||
		//			(bestDist > (thisRead2 - thisRead1))) {
		//			bestDist = (thisRead2 - thisRead1);
		//			best = other + 1;
		//		}
		//	}
		//}

		//_streamPutInt(target, best);
		_streamPut(target, '0');
		_streamPut(target, '\n');
		_streamWrite(target,"pos:");
		_streamPutInt(target, clr1);
		_streamPut(target, ',');
		_streamPutInt(target, clr2);
		_streamPut(target, '\n');
		_streamWrite(target,"dln:0\n");
		_streamWrite(target,"del:\n");
		_streamWrite(target,"}\n");
	}
	_streamWrite(target,"}\n");
}



}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...

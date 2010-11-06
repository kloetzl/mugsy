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
  $Id: graph_consensus_ctgstore.h 1811 2008-03-31 15:38:54Z rausch@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_GRAPH_CONSENSUS_CTGSTORE_H
#define SEQAN_HEADER_GRAPH_CONSENSUS_CTGSTORE_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Contig Store
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec = Default>
class GappedRead;

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec>
class GappedRead
{
	public:
		typedef typename Size<GappedRead>::Type TSize;
		typedef typename Id<GappedRead>::Type TId;
		TSize data_offset;
		TId data_source;
		Pair<TSize, TSize> data_clr;
		String<TSize> data_gap;
		
	public:
		GappedRead() : data_offset(0), data_source(0)
		{
			SEQAN_CHECKPOINT
			clear(data_gap);
		}

		GappedRead(GappedRead const & _other)
		{
			SEQAN_CHECKPOINT
			data_clr = _other.data_clr;
			data_gap = _other.data_gap;
			data_offset = _other.data_offset;
			data_source = _other.data_source;
		}

		GappedRead const& 
		operator = (GappedRead const& _other) 
		{
			SEQAN_CHECKPOINT
			if (this == &_other) return *this;
			data_clr = _other.data_clr;
			data_gap = _other.data_gap;
			data_offset = _other.data_offset;
			data_source = _other.data_source;
			return *this;
		}

};


//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet = char, typename TSpec = Default>
class CtgStore;

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TSpec>
class CtgStore
{
	public:
		typedef typename Size<CtgStore>::Type TSize;
		String<TAlphabet, External<> > data_contig;
		String<char, External<> > data_quality;
		StringSet<String<char, External<> >, Owner<ConcatDirect<> > > data_names;
		String<Pair<TSize, TSize> > data_begin_end;
		String<String<GappedRead<> > > data_reads;
		TSize data_pos_count;
		
	public:
		CtgStore() : data_pos_count(0)
		{
			SEQAN_CHECKPOINT
			clear(data_contig);
			clear(data_quality);
			clear(data_names);
			clear(data_begin_end);
			clear(data_reads);
		}

		~CtgStore() 
		{
			SEQAN_CHECKPOINT
		}


	private:
		CtgStore(CtgStore const & _other)
		{
			SEQAN_CHECKPOINT
			data_contig = _other.data_contig;
			data_quality = _other.data_quality;
			data_names = _other.data_names;
			data_begin_end = _other.data_begin_end;
			data_reads = _other.data_reads;
			data_pos_count = _other.data_pos_count;
		}

		CtgStore const& 
		operator = (CtgStore const& _other) 
		{
			SEQAN_CHECKPOINT
			if (this == &_other) return *this;
			data_contig = _other.data_contig;
			data_quality = _other.data_quality;
			data_names = _other.data_names;
			data_begin_end = _other.data_begin_end;
			data_reads = _other.data_reads;
			data_pos_count = _other.data_pos_count;
			return *this;
		}
};
	

//////////////////////////////////////////////////////////////////////////////
// FUNCTIONS
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TSpec, typename TSize, typename TContig>
inline void 
loadContig(CtgStore<TAlphabet, TSpec>& ctgSt, 
		   TSize index,
		   TContig& seq) 
{
	SEQAN_CHECKPOINT
	clear(seq);
	seq = infix(ctgSt.data_contig, (value(ctgSt.data_begin_end, index)).i1, (value(ctgSt.data_begin_end, index)).i2);
}




//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TSpec, typename TSize, typename TContig>
inline void 
loadQuality(CtgStore<TAlphabet, TSpec>& ctgSt,
			TSize index,
			TContig& seq) 
{
	SEQAN_CHECKPOINT
	clear(seq);
	seq = infix(ctgSt.data_quality, (value(ctgSt.data_begin_end, index)).i1, (value(ctgSt.data_begin_end, index)).i2);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TSpec, typename TSize, typename TName>
inline void 
loadExtId(CtgStore<TAlphabet, TSpec>& ctgSt,
		  TSize index,
		  TName& name) 
{
	SEQAN_CHECKPOINT
	name = value(ctgSt.data_names, index);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TSpec>
inline typename Size<CtgStore<TAlphabet, TSpec> >::Type
length(CtgStore<TAlphabet, TSpec>& ctgSt) 
{
	SEQAN_CHECKPOINT
	return ctgSt.data_pos_count;
}



}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...

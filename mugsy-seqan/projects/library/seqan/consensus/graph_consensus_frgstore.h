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
  $Id: graph_consensus_frgstore.h 1811 2008-03-31 15:38:54Z rausch@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_GRAPH_CONSENSUS_FRGSTORE_H
#define SEQAN_HEADER_GRAPH_CONSENSUS_FRGSTORE_H

namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////
// Fragment Store
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec = Default>
class FrgStore;

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec>
class FrgStore
{
	public:
		typedef typename Size<FrgStore>::Type TSize;
		StringSet<String<char, External<> >, Owner<ConcatDirect<> > > data_names;
		String<TSize> data_lib_id;
		String<Pair<TSize, TSize> > data_rds;
		TSize data_pos_count;
		
	public:
		FrgStore() : data_pos_count(0)
		{
			SEQAN_CHECKPOINT
			clear(data_names);
			clear(data_lib_id);
			clear(data_rds);
		}

		~FrgStore() 
		{
			SEQAN_CHECKPOINT
		}

	private:
		FrgStore(FrgStore const & _other)
		{
			data_names = _other.data_names;
			data_lib_id = _other.data_lib_id;
			data_rds = _other.data_rds;
			data_pos_count = _other.data_pos_count;
		}

		FrgStore const& 
		operator = (FrgStore const& _other) 
		{
			if (this == &_other) return *this;
			data_names = _other.data_names;
			data_lib_id = _other.data_lib_id;
			data_rds = _other.data_rds;
			data_pos_count = _other.data_pos_count;
			return *this;
		}
};
	


//////////////////////////////////////////////////////////////////////////////
// FUNCTIONS
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////


template<typename TSpec, typename TSize, typename TName>
inline void 
loadExtId(FrgStore<TSpec>& frgSt,
		  TSize index,
		  TName& name) 
{
	SEQAN_CHECKPOINT
	name = value(frgSt.data_names, index);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TSize>
inline typename Id<FrgStore<TSpec> >::Type
loadLibraryId(FrgStore<TSpec>& frgSt,
			  TSize index) 
{
	SEQAN_CHECKPOINT
	return value(frgSt.data_lib_id, index);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TSize, typename TId1, typename TId2>
inline void
loadRdsId(FrgStore<TSpec>& frgSt,
		  TSize index,
		  TId1& id1,
		  TId2& id2) 
{
	SEQAN_CHECKPOINT
	id1 = (value(frgSt.data_rds,index)).i1;
	id2 = (value(frgSt.data_rds,index)).i2;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec>
inline typename Size<FrgStore<TSpec> >::Type
length(FrgStore<TSpec>& frgSt) 
{
	SEQAN_CHECKPOINT
	return frgSt.data_pos_count;
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...

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
  $Id: graph_consensus_libstore.h 1811 2008-03-31 15:38:54Z rausch@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

#ifndef SEQAN_HEADER_GRAPH_CONSENSUS_LIBSTORE_H
#define SEQAN_HEADER_GRAPH_CONSENSUS_LIBSTORE_H

namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////
// Library Store
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec = Default>
class LibStore;

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec>
class LibStore
{
	public:
		typedef typename Size<LibStore>::Type TSize;
		StringSet<String<char, External<> >, Owner<ConcatDirect<> > > data_names;
		String<TSize> data_mean;
		String<TSize> data_std;
		TSize data_pos_count;
		
	public:
		LibStore() : data_pos_count(0)
		{
			SEQAN_CHECKPOINT
			clear(data_names);
			clear(data_mean);
			clear(data_std);
		}

		~LibStore() 
		{
			SEQAN_CHECKPOINT
		}

	private:
		LibStore(LibStore const & _other)
		{
			data_names = _other.data_names;
			data_mean = _other.data_mean;
			data_std = _other.data_std;
			data_pos_count = _other.data_pos_count;
		}

		LibStore const& 
		operator = (LibStore const& _other) 
		{
			if (this == &_other) return *this;
			data_names = _other.data_names;
			data_mean = _other.data_mean;
			data_std = _other.data_std;
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
loadExtId(LibStore<TSpec>& libSt,
		  TSize index,
		  TName& name) 
{
	SEQAN_CHECKPOINT
	name = value(libSt.data_names, index);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TSize>
inline typename Id<LibStore<TSpec> >::Type
loadLibraryMean(LibStore<TSpec>& libSt,
				TSize index) 
{
	SEQAN_CHECKPOINT
	return value(libSt.data_mean, index);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TSize>
inline typename Id<LibStore<TSpec> >::Type
loadLibraryStd(LibStore<TSpec>& libSt,
			   TSize index) 
{
	SEQAN_CHECKPOINT
	return value(libSt.data_std, index);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec>
inline typename Size<LibStore<TSpec> >::Type
length(LibStore<TSpec>& libSt) 
{
	SEQAN_CHECKPOINT
	return libSt.data_pos_count;
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...

//SEQAN_NO_GENERATED_FORWARDS: no forwards are generated for this file

#include <seqan/basic.h>
#include "msufsort.d/MSufSort/MSufSort.h"

namespace SEQAN_NAMESPACE_MAIN {

	namespace shawarma {

		void msufsort(unsigned char *buffer, unsigned int nBytes) 
		{
			::MSufSort * msufsort = new MSufSort;
			msufsort->Sort(buffer, nBytes);
			delete msufsort;
		}

		unsigned int msufbwt(unsigned char *buffer, unsigned int nBytes) 
		{
			::MSufSort * msufsort = new MSufSort;
			unsigned int index;
			msufsort->BWT(buffer, nBytes, index);
			delete msufsort;
			return index;
		}

	}

} // namespace SEQAN_NAMESPACE_MAIN::msufsort

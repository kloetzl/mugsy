//SEQAN_NO_GENERATED_FORWARDS: no forwards are generated for this file

#include <seqan/basic.h>

extern "C" {
	#include "ds.d/ds_ssort.h"
}

namespace SEQAN_NAMESPACE_MAIN {

	namespace shawarma {

		void ds_ssort(unsigned char *t, int *sa, int n) {
			::ds_ssort(t, sa, n);
		}

		int init_ds_ssort(int adist, int bs_ratio) {
			return ::init_ds_ssort(adist, bs_ratio);
		}

	}

} // namespace SEQAN_NAMESPACE_MAIN::ds

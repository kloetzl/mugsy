//SEQAN_NO_GENERATED_FORWARDS: no forwards are generated for this file

#include <seqan/basic.h>

extern "C" {
	void suffixsort(int *x, int *p, int n, int k, int l);
}

namespace SEQAN_NAMESPACE_MAIN {

	namespace shawarma {

		void qsufsort(int *x, int *p, int n, int k, int l) {
			::suffixsort(x, p, n, k, l);
		}

	}

} // namespace SEQAN_NAMESPACE_MAIN::ds

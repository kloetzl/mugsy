//SEQAN_NO_GENERATED_FORWARDS: no forwards are generated for this file

#include <seqan/basic.h>
#include "divsufsort.d/include/divsufsort.h"

namespace SEQAN_NAMESPACE_MAIN {

	namespace shawarma {

		int divsufsort(unsigned char *t, int *sa, int n) {
			return ::divsufsort(t, sa, n);
		}

		int divbwt(unsigned char *t, unsigned char *u, int *a, int n) {
			return ::divbwt(t, u, a, n);
		}

		const char *divsufsort_version(void) {
			return ::divsufsort_version();
		}

	}

} // namespace SEQAN_NAMESPACE_MAIN::ds

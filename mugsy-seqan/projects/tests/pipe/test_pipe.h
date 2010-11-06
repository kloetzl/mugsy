/*
 *  test_pipe.h
 *  genindex
 *
 *  Created by David Weese on 17.07.05.
 *
 */

#ifndef SEQAN_HEADER_TEST_PIPE_H
#define SEQAN_HEADER_TEST_PIPE_H

namespace SEQAN_NAMESPACE_MAIN
{

	template < typename TBuffer >
	void permute(TBuffer buf) {
		typename Size<TBuffer>::Type i, j, s = length(buf);
//        srand( (unsigned)time( NULL ) );
		for(i = 0; i < s; i++)
			buf[i] = s-i-1;
		for(i = 0; i < s; i++) {
            if (i > 0) {
			    j = i - (rand() % i) - 1;
                assert(/*0 <= j && */j < s);
            } else
                j = 0;
			unsigned tmp = buf[i];
			buf[i] = buf[j];
			buf[j] = tmp;
		}
	}

	template < typename TBuffer >
	void randomize(TBuffer buf) {
		typename Size<TBuffer>::Type i, s = length(buf);
		for(i = 0; i < s; i++)
            buf[i] = rand() % s;
	}

	template < typename TValue >
	struct IdentityMap : public ::std::unary_function< TValue, TValue > {
		inline TValue operator() (TValue const i) { return i; }
	};

	template < typename TValue >
	struct SimpleCompare : public ::std::binary_function< TValue const, TValue const, int > {
		inline int operator() (TValue const a, TValue const b) const {
            if (a < b) return -1;
            if (a > b) return 1;
            return 0;
        }
	};
	
}

#endif


#include <iostream>
#include <fstream>
#include <typeinfo>
#include <time.h>
#include <cstdio>
#include <vector>

#define SEQAN_DEBUG
//#define SEQAN_TEST

#include <seqan/pipe.h>
#include "test_pipe.h"

using namespace std;
using namespace seqan;


//////////////////////////////////////////////////////////////////////////////


	bool testExternalString(unsigned maxSize = 16*1024*1024) {
        typedef Iterator< String<unsigned,External<> > const >::Type TIter;
		SimpleBuffer<unsigned> buf;
		allocPage(buf, maxSize, buf);

		String<unsigned,External<> > vector;
		for(unsigned i = 1; i <= maxSize; i = i << 1) {
			::std::cout << i << " "; ::std::cout.flush();
			resize(buf, i);
			randomize(buf);

			Pipe<SimpleBuffer<unsigned>, Source<> > src(buf);
			vector << src;

			TIter I = begin(vector);
            for(unsigned *cur = buf.begin; cur != buf.end; ++cur) {
                if (*cur != *I) {
                    ::std::cout << ::std::endl << "testExternalString failed at position " << (cur - buf.begin) << " ";
                    return false;
                }
				++I;
			}
		}
		freePage(buf, buf);
		return true;
	}



	bool testPool(unsigned maxSize = 16*1024*1024) {
		SimpleBuffer<unsigned> buf;
		allocPage(buf, maxSize, buf);

		Pool<unsigned,PoolSpec<> > pool;
		for(unsigned i = 1; i <= maxSize; i = i << 1) {
			::std::cout << i;
			::std::cout.flush();

			resize(buf, i);
			randomize(buf);

			Pipe<SimpleBuffer<unsigned>, Source<> > src(buf);
			pool << src;

			if (pool.memBuffer.begin)
				::std::cout << "* ";
			else 
				::std::cout << " ";
			::std::cout.flush();

			beginRead(pool);
            for(unsigned *cur = buf.begin; cur != buf.end; cur++) {
                if (*cur != *pool) {
			        endRead(pool);
                    ::std::cout << ::std::endl << "testPool failed at position " << (cur - buf.begin) << " ";
                    return false;
                }
				++pool;
			}
			endRead(pool);
		}
		freePage(buf, buf);
		return true;
	}



	bool testMapper(unsigned maxSize = 16*1024*1024) {
		SimpleBuffer<unsigned> buf;
		allocPage(buf, maxSize, buf);

		Pool<unsigned,MapperSpec<MapperConfig<IdentityMap<unsigned> > > > mapper;
		for(unsigned i = 1; i <= maxSize; i = i << 1) {
			::std::cout << i;
			::std::cout.flush();

			resize(buf, i);
			permute(buf);

			Pipe<SimpleBuffer<unsigned>, Source<> > src(buf);
			mapper << src;

			if (mapper.memBuffer.begin)
				::std::cout << "* ";
			else 
				::std::cout << " ";
			::std::cout.flush();

			beginRead(mapper);
			for(unsigned j = 0; j < i; ++j) {
				if (*mapper != j) {
					freePage(buf, buf);
                    ::std::cout << ::std::endl << "testMapper failed at position " << j << " ";
					return false;
				}	
				++mapper;
			}
			endRead(mapper);
		}
		freePage(buf, buf);
		return true;
	}



	bool testPartiallyFilledMapper(unsigned maxSize = 16*1024*1024) {
		SimpleBuffer<unsigned> buf;
		allocPage(buf, maxSize, buf);

		Pool<unsigned,MapperSpec<MapperConfig<IdentityMap<unsigned> > > > mapper;
		for(unsigned i = 1; i <= maxSize; i = i << 1) {
			::std::cout << i;
			::std::cout.flush();

			resize(buf, i);
			permute(buf);

			// partially fill the mapper 
			mapper.undefinedValue = i;	// select i as an undefined value (all defined values are less than i)
			resize(mapper, i);
			resize(buf, i - i/3);
			Pipe<SimpleBuffer<unsigned>, Source<> > src(buf);
			beginWrite(mapper) && append(mapper, src) && endWrite(mapper);

			if (mapper.memBuffer.begin)
				::std::cout << "* ";
			else 
				::std::cout << " ";
			::std::cout.flush();

			unsigned undefCounter = 0, missCounter = 0;
            beginRead(mapper);
			for(unsigned j = 0; j < i; ++j) {
				if (*mapper == i) 
					++undefCounter;
				else
					if (*mapper != j) {
						++missCounter;
						if (!mapper.memBuffer.begin) { // external mapping -> no misses allowed
							freePage(buf, buf);
							::std::cout << ::std::endl << "testPartiallyFilledMapper failed at position " << j << " [ = " << *mapper << " ] ";
							return false;
						}
					}
				++mapper;
			}
			endRead(mapper);
			if (mapper.memBuffer.begin) {
				if (undefCounter + missCounter > i/3) {
					::std::cout << ::std::endl << "testPartiallyFilledMapper failed [only " << undefCounter + missCounter << " of " << i/3 << " undefined] ";
					return false;
				}
			} else
				if (undefCounter != i/3) {
					::std::cout << ::std::endl << "testPartiallyFilledMapper failed [only " << undefCounter + missCounter << " of " << i/3 << " undefined] ";
					return false;
				}
		}
		freePage(buf, buf);
		return true;
	}



	bool testSorter(unsigned maxSize = 16*1024*1024) {
		SimpleBuffer<unsigned> buf;
		allocPage(buf, maxSize, buf);

		Pool<unsigned,SorterSpec<SorterConfig<SimpleCompare<unsigned> > > > sorter;
		for(unsigned i = 1; i <= maxSize; i = i << 1) {
			::std::cout << i;
			::std::cout.flush();

			resize(buf, i);
			permute(buf);

			Pipe<SimpleBuffer<unsigned>, Source<> > src(buf);
			sorter << src;

			if (sorter.memBuffer.begin)
				::std::cout << "* ";
			else 
				::std::cout << " ";
			::std::cout.flush();

			beginRead(sorter);
            unsigned j = *sorter, pos = 0;
			while (!eof(sorter)) {
				if (*sorter < j || /* *sorter < 0 || */ *sorter >= i) {
					freePage(buf, buf);
                    ::std::cout << ::std::endl << "testSorter failed at position " << pos << " ";
					return false;
				}
                j = *sorter;
				++sorter; ++pos;
			}
			endRead(sorter);
		}
		freePage(buf, buf);
		return true;
	}



    bool testContainers(unsigned maxSize = 32*1024*1024) {
		bool result = true;
        time_t start, pool, mapperPart, mapper, sorter;
        time(&start);
//____________________________________________________________________________

		::std::cout << "Testing External Vector ... " << ::std::endl;
		if (testExternalString(maxSize) || (result = false)) {
			::std::cout << "OK.";
		} else
			::std::cout << "FAILED.";
		::std::cout << ::std::endl;
//____________________________________________________________________________

		::std::cout << "Testing Simple Pool ... " << ::std::endl;
		if (testPool(maxSize) || (result = false)) {
			::std::cout << "OK.";
		} else
			::std::cout << "FAILED.";

		time(&pool);
		::std::cout << " " << (pool - start) << " seconds" << ::std::endl;
//____________________________________________________________________________

		::std::cout << "Testing Mapper ... " << ::std::endl;
		if (testMapper(maxSize) || (result = false)) {
			::std::cout << "OK.";
		} else
			::std::cout << "FAILED.";

		time(&mapper);
		::std::cout << " " << (mapper - pool) << " seconds" << ::std::endl;
//____________________________________________________________________________

		::std::cout << "Testing partially filled Mapper ... " << ::std::endl;
		if (testPartiallyFilledMapper(maxSize) || (result = false)) {
			::std::cout << "OK.";
		} else
			::std::cout << "FAILED.";

        time(&mapperPart);
		::std::cout << " " << (mapperPart - mapper) << " seconds" << ::std::endl;
//____________________________________________________________________________

		::std::cout << "Testing Sorter ... " << ::std::endl;
		if (testSorter(maxSize) || (result = false)) {
			::std::cout << "OK.";
		} else
			::std::cout << "FAILED.";

        time(&sorter);
		::std::cout << " " << (sorter - mapperPart) << " seconds" << ::std::endl;
//____________________________________________________________________________

		return result;
	}



int main() 
{
//	SEQAN_TREPORT("TEST BEGIN")

	testContainers();

//	debug::verifyCheckpoints("projects/library/seqan/pipe/pool_base.h");

//	SEQAN_TREPORT("TEST END")

	return 0;
}

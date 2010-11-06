#include <iostream>

#define SEQAN_DEBUG
#define SEQAN_TEST

#include "seqan/basic.h"

#include <memory>
#include <vector>
#include <map>

using namespace std;
using namespace seqan;

//____________________________________________________________________________

struct TestAllocator
{
	mutable map<char *, size_t> data_allocated;
	mutable map<char *, size_t> data_deallocated;

	TestAllocator() {}
	~TestAllocator() 
	{
		map<char *, size_t>::iterator it = data_allocated.begin();
		while (it != data_allocated.end())
		{
			SEQAN_TASSERT2(data_deallocated.count(it->first), "memory block not deallocated");
            deallocate(int(), it->first, it->second);
			++it;
		}
	}
};

template <typename TValue, typename TSize, typename TUsage>
void allocate(TestAllocator & me, 
			  TValue * & data_, 
			  TSize count, 
			  Tag<TUsage> const)
{
	SEQAN_ASSERT(count)
	allocate(int(), data_, count);
	me.data_allocated[(char *) data_] = count;
}

template <typename TValue, typename TSize, typename TUsage>
void deallocate(TestAllocator & me, 
				TValue * data_, 
				TSize count, 
				Tag<TUsage> const)
{
	SEQAN_TASSERT2(me.data_allocated.count((char *) data_), "memory block was not allocated")
	SEQAN_TASSERT2(me.data_allocated[(char *) data_] == count, "memory block was allocated with different size");
	SEQAN_TASSERT2(!me.data_deallocated.count((char *) data_), "memory block already deallocated");

	me.data_deallocated[(char *) data_] = count;
}

int countAllocs(TestAllocator & me)
{
	return me.data_allocated.size();
}
int countDeallocs(TestAllocator & me)
{
	return me.data_deallocated.size();
} 
 

//____________________________________________________________________________

/*
void test()
{
	printf("start Test\n");

	int somewhat;

	char * s;
	allocate(somewhat, s, 100);
	deallocate(somewhat, s, 100);

	ToStdAllocator<int, char> std_alloc(somewhat);
	vector<char, ToStdAllocator<int, char> > vec( std_alloc );

	char c;
	for (c='A'; c <= 'Z'; ++c)
	{
		vec.push_back(c);
	}

	printf("\nend Test\n");
}
*/
//____________________________________________________________________________

void testSimpleAllocator()
{
	int * dat1;
	int * dat2;

	Allocator<SimpleAlloc<TestAllocator> > allo1;
	allocate(allo1, dat1, 100);
	allocate(allo1, dat2, 105);
	deallocate(allo1, dat1, 100);
	allocate(allo1, dat2, 201);

	SEQAN_TASSERT(countAllocs(parentAllocator(allo1)) == 3);
	SEQAN_TASSERT(countDeallocs(parentAllocator(allo1)) == 1);

	clear(allo1);

	SEQAN_TASSERT(countDeallocs(parentAllocator(allo1)) == 3);
}
//____________________________________________________________________________

void testPoolAllocator()
{
	int * dat1;
	int * dat2;

	typedef Allocator<SimpleAlloc<TestAllocator> > TParentAlloc;
	Allocator<SinglePool<20 * sizeof(int), TParentAlloc> > allo1;
	allocate(allo1, dat1, 20);
	allocate(allo1, dat2, 20);
	deallocate(allo1, dat1, 20);
	allocate(allo1, dat2, 20);

	SEQAN_TASSERT(dat1 == dat2)

	SEQAN_TASSERT(countAllocs(parentAllocator(parentAllocator(allo1))) == 1);
	SEQAN_TASSERT(countDeallocs(parentAllocator(parentAllocator(allo1))) == 0);

	allocate(allo1, dat1, 100);
	deallocate(allo1, dat1, 100);

	SEQAN_TASSERT(countAllocs(parentAllocator(parentAllocator(allo1))) == 2);
	SEQAN_TASSERT(countDeallocs(parentAllocator(parentAllocator(allo1))) == 1);

	clear(allo1);

	SEQAN_TASSERT(countDeallocs(parentAllocator(parentAllocator(allo1))) == 2);
}

//____________________________________________________________________________

void testMultiPoolAllocator()
{
	int * dat1;
	int * dat2;

	typedef Allocator<SimpleAlloc<TestAllocator> > TParentAlloc;
	Allocator<MultiPool<TParentAlloc> > allo1;
	allocate(allo1, dat1, 20);
	allocate(allo1, dat2, 20);
	deallocate(allo1, dat1, 20);
	allocate(allo1, dat2, 20);

	SEQAN_TASSERT(dat1 == dat2)

	SEQAN_TASSERT(countAllocs(parentAllocator(parentAllocator(allo1))) == 1);
	SEQAN_TASSERT(countDeallocs(parentAllocator(parentAllocator(allo1))) == 0);

	allocate(allo1, dat1, 30);
	deallocate(allo1, dat1, 30);

	SEQAN_TASSERT(countAllocs(parentAllocator(parentAllocator(allo1))) == 2);
	SEQAN_TASSERT(countDeallocs(parentAllocator(parentAllocator(allo1))) == 0);

	clear(allo1);

	SEQAN_TASSERT(countDeallocs(parentAllocator(parentAllocator(allo1))) == 2);
}

//____________________________________________________________________________

void Main_Test_Allocator() 
{ 
	SEQAN_TREPORT("TEST ALLOCATOR BEGIN")

	testSimpleAllocator();
	testPoolAllocator();
	testMultiPoolAllocator();

	debug::verifyCheckpoints("projects/library/seqan/basic/basic_allocator_interface.h");
	debug::verifyCheckpoints("projects/library/seqan/basic/basic_allocator_to_std.h");
	debug::verifyCheckpoints("projects/library/seqan/basic/basic_allocator_simple.h");
	debug::verifyCheckpoints("projects/library/seqan/basic/basic_allocator_singlepool.h");
	debug::verifyCheckpoints("projects/library/seqan/basic/basic_allocator_multipool.h");

	SEQAN_TREPORT("TEST ALLOCATOR END")
}

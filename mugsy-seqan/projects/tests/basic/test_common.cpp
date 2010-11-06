#include <iostream>
#include <cstdio>
#include <vector>

#define SEQAN_DEBUG
#define SEQAN_TEST

#include "seqan/basic.h"
#include "seqan/sequence.h"

using namespace std;
using namespace seqan;

//////////////////////////////////////////////////////////////////////////////


void Test_Definition()
{
	SEQAN_ASSERT(_ClassIdentifier<int>::getID() == _ClassIdentifier<int>::getID());
	SEQAN_ASSERT(_ClassIdentifier<char>::getID() != _ClassIdentifier<int>::getID());
}

//////////////////////////////////////////////////////////////////////////////

void Test_Type()
{
	int i;
	int const ci = 99;
	int a[10];

	_toParameter<int>(& i) = 10;
	SEQAN_ASSERT(i == 10);

	*_toParameter<int *>(& i) = 20;
	SEQAN_ASSERT(i == 20);

	_Pointer<int>::Type p1 = _toPointer(i);
	*p1 = 30;
	SEQAN_ASSERT(i == 30);

	_Pointer<int *>::Type p2 = _toPointer(p1);
	*p2 = 40;
	SEQAN_ASSERT(i == 40);

	_Pointer<int[10]>::Type p3 = _toPointer(a);
	p3[1] = 50;
	SEQAN_ASSERT(a[1] == 50);

	_Pointer<int const *>::Type p4 = _toPointer(ci);
	SEQAN_ASSERT(*p4 == 99);

}

//////////////////////////////////////////////////////////////////////////////

void Test_Iterator_Adapt_Std()
{
//test SeqAn iterator to fulfill std iterator 

	typedef ::std::iterator_traits<Iterator<char *, Rooted>::Type>::value_type T1;
	bool b1 = _isSameType<T1, char>();
	SEQAN_TASSERT(b1)

	
}

//////////////////////////////////////////////////////////////////////////////

void Main_Test_Common() 
{
	SEQAN_TREPORT("TEST COMMON BEGIN")

	Test_Definition();
	Test_Type();
	Test_Iterator_Adapt_Std();

	debug::verifyCheckpoints("projects/library/seqan/basic/basic_definition.h");
	debug::verifyCheckpoints("projects/library/seqan/basic/basic_type.h");

	debug::verifyCheckpoints("projects/library/seqan/basic/basic_iterator_adapt_std.h");

	SEQAN_TREPORT("TEST COMMON END")
}

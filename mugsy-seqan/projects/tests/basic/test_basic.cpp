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
//test iterator proxy

void Test_Proxy_Iterator()
{
	int i1[] = {10, 20, 30};
	int * pi1 = i1;
	Proxy<IteratorProxy<int *> > px(pi1);
	SEQAN_TASSERT(px == 10);

//assign
	px = 11;
	SEQAN_TASSERT(i1[0] == 11);

	int i2 = 12;
	px = i2;
	SEQAN_TASSERT(i1[0] == 12);

	int * pi2 = i1 + 1;
	Proxy<IteratorProxy<int *> > px2(pi2);
	px = px2;
	SEQAN_TASSERT(i1[0] == 20);
	SEQAN_TASSERT(px == 20);

//copy ctor
	Proxy<IteratorProxy<int *> > px3(px2);
	SEQAN_TASSERT(px3 == 20);


//assign
	char s1[100] = "";
	char * it1 = s1;
	Proxy<IteratorProxy<char *> > px4(it1);

	assign(px4, 'X');
	SEQAN_TASSERT(px4 == 'X');

	char c1 = 'a';
	assign(px4, c1);
	SEQAN_TASSERT(px4 == 'a');

}


//////////////////////////////////////////////////////////////////////////////
//helper class for reference counting
//this is needed for test of Holder

struct RefCountObj
{
	static int static_addrefs;
	static int static_releaserefs;
	static int static_ctors;
	static int static_dtors;

	mutable int data_addrefs;
	mutable int data_releaserefs;

	int data_value;

	RefCountObj():
		data_addrefs(0),
		data_releaserefs(0),
		data_value(0)
	{
		++static_ctors;
	}
	RefCountObj(RefCountObj const & other_):
		data_addrefs(0),
		data_releaserefs(0),
		data_value(other_.data_value)
	{
		++static_ctors;
	}
	~RefCountObj()
	{
SEQAN_TASSERT(data_addrefs == data_releaserefs);

		++static_dtors;
	}
	RefCountObj & operator = (RefCountObj const & other_)
	{
		data_addrefs = data_releaserefs = 0;
		data_value = other_.data_value;
		return *this;
	}
};

int RefCountObj::static_addrefs = 0;
int RefCountObj::static_releaserefs = 0;
int RefCountObj::static_ctors = 0;
int RefCountObj::static_dtors = 0;

void addRef(RefCountObj & me)
{
	++me.data_addrefs;
	++RefCountObj::static_addrefs;
}
void addRef(RefCountObj const & me)
{
	++me.data_addrefs;
	++RefCountObj::static_addrefs;
}
void releaseRef(RefCountObj & me)
{
	++me.data_releaserefs;
	++RefCountObj::static_releaserefs;
}
void releaseRef(RefCountObj const & me)
{
	++me.data_releaserefs;
	++RefCountObj::static_releaserefs;
}

//____________________________________________________________________________
// test for holder class

void Test_Holder()
{
	{
//ctors
		Holder<RefCountObj> ho1;
		SEQAN_TASSERT(empty(ho1));
		SEQAN_TASSERT(!dependent(ho1));

		create(ho1);
		SEQAN_TASSERT(!empty(ho1));
		SEQAN_TASSERT(!dependent(ho1));

		Holder<RefCountObj> ho2(ho1);

		Holder<RefCountObj> ho3(value(ho1));
		SEQAN_TASSERT(value(ho1).data_addrefs == 1);
		SEQAN_TASSERT(value(ho1).data_releaserefs == 0);
		SEQAN_TASSERT(!dependent(ho1));


//create
		RefCountObj rco1;
		create(ho3, rco1);
		SEQAN_TASSERT(value(ho3).data_addrefs == 0);
		SEQAN_TASSERT(value(ho1).data_addrefs == 1);
		SEQAN_TASSERT(value(ho1).data_releaserefs == 1);

//setValue
		setValue(ho3, rco1);
		SEQAN_TASSERT(dependent(ho3));
		SEQAN_TASSERT(value(ho3).data_addrefs == 1);

		rco1.data_value = 10;
		create(ho3);
		SEQAN_TASSERT(value(ho3).data_value == 10);

		RefCountObj rco2;
		rco2.data_value = 20;

//operator = (value) => assignValue
		ho2 = rco2;
		SEQAN_TASSERT(value(ho2).data_value == 20);
		SEQAN_TASSERT(rco2.data_addrefs == 0);
		SEQAN_TASSERT(!dependent(ho2));

		rco2.data_value = 30;
		SEQAN_TASSERT(value(ho2).data_value == 20);

//operator = (holder) => assign
		setValue(ho1, rco1);
		ho1 = ho2;
		SEQAN_TASSERT(value(ho1).data_value == 20);
		SEQAN_TASSERT(rco2.data_addrefs == 0);

//clear
		clear(ho3);
		SEQAN_TASSERT(empty(ho3));

		assign(ho2, ho3);
		SEQAN_TASSERT(empty(ho2));

//conversion operator
		rco1 = ho1;
		SEQAN_TASSERT(rco1.data_value == 20);

//moveValue
		moveValue(ho1, rco2);
		SEQAN_TASSERT(rco1.data_value == 30);

	}

	SEQAN_TASSERT(RefCountObj::static_addrefs == RefCountObj::static_releaserefs);
	SEQAN_TASSERT(RefCountObj::static_ctors == RefCountObj::static_dtors);


//test default implementations of addRef and releaseRef

	int i1;
	int const i2 = 0;

	addRef(i1);
	addRef(i2);
	releaseRef(i1);
	releaseRef(i2);

//test const object holders
/*
	typedef char Bla[100];
	Holder<Bla const> cho1 = "test";*/
}

//////////////////////////////////////////////////////////////////////////////
// Test Iterator, basic functions

namespace SEQAN_NAMESPACE_MAIN
{

struct Test_Iterator_1
{
	int data_dat_1;
	mutable int data_dat_2;

	Test_Iterator_1(int data_):
		data_dat_1(data_ + 1),
		data_dat_2(data_ + 2)
	{
	}
	Test_Iterator_1(Test_Iterator_1 const & other_):
		data_dat_1(other_.data_dat_1 + 10),
		data_dat_2(other_.data_dat_2 + 10)
	{
	}
	~Test_Iterator_1() 
	{
	}

	Test_Iterator_1 & operator = (Test_Iterator_1 const & other_)
	{
		data_dat_1 = other_.data_dat_1 + 20;
		data_dat_2 = other_.data_dat_2 + 20;
		return *this;
	}

	int & operator * ()
	{
		return data_dat_1;
	}
	int & operator * () const
	{
		return data_dat_2;
	}
};

template <>
struct Value<Test_Iterator_1>
{
	typedef int Type;
};

}
//____________________________________________________________________________


void Test_Iterator_Basic()
{
//test default iterator functions

	//value, getValue, operator *
	int i1 = 10;
	int * it1 = & i1;
	SEQAN_TASSERT(value(it1) == 10);
	SEQAN_TASSERT(getValue(it1) == 10);

	Test_Iterator_1 it2(10);
	SEQAN_TASSERT(value(it2) == 11);
	SEQAN_TASSERT(getValue(it2) == 11);
/*
	Test_Iterator_1 const it3(10);
	SEQAN_TASSERT(value(it3) == 12);
	SEQAN_TASSERT(getValue(it3) == 12);

	//assign to value reference
	value(it2) = 15;
	SEQAN_TASSERT(getValue(it2) == 15);

	//moveValue
	moveValue(it2, 50);
	SEQAN_TASSERT(value(it2) == 50);

	//defaults of some advanced functions
	SEQAN_TASSERT(position(it1) == 0);
	SEQAN_TASSERT(container(it1) == *it1);
*/
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec>
void Test_Iter()
{
	typedef Iter<char *, TSpec> TIterator;

	char arr1[] = "XYZ";
	char arr2[] = "abcdefg";

	TIterator it1 = begin(arr1);
	SEQAN_TASSERT(container(it1) == arr1);
	SEQAN_TASSERT(*it1 == 'X')

	setContainer(it1, arr2);
	SEQAN_TASSERT(*it1 == 'a')

	TIterator it2(it1);
	SEQAN_TASSERT(*it2 == 'a')

	TIterator it3;
	it3 = it2 + 1;
	SEQAN_TASSERT(*it3 == 'b')

	++it3;
	SEQAN_TASSERT(*it3 == 'c')

	it3++;
	SEQAN_TASSERT(*it3 == 'd')

	SEQAN_TASSERT(position(it3) == 3)

	setPosition(it3, 2);
	SEQAN_TASSERT(*it3 == 'c')

	--it3;
	SEQAN_TASSERT(*it3 == 'b')

	it3--;
	SEQAN_TASSERT(*it3 == 'a')
	SEQAN_TASSERT(it3 == it2)
	SEQAN_TASSERT(!(it3 != it2))

	++it3;
	assignValue(it3, 'z');
	SEQAN_TASSERT(*it3 == 'z')

	moveValue(it3, 'y');
	SEQAN_TASSERT(*it3 == 'y')

	TIterator const it4 = it3;
	assignValue(it4, 'x');
	SEQAN_TASSERT(*it4 == 'x')
	SEQAN_TASSERT(*it3 == 'x')

	moveValue(it4, 'w');
	SEQAN_TASSERT(*it4 == 'w')
	SEQAN_TASSERT(*it3 == 'w')


	it3 = 1 + it4;
	SEQAN_TASSERT(*it3 == 'c')

	it3 = it4 - 1;
	SEQAN_TASSERT(*it3 == 'a')

	it3 += 4;
	SEQAN_TASSERT(*it3 == 'e')

	it3 -= 2;
	SEQAN_TASSERT(*it3 == 'c')

	SEQAN_TASSERT(it3 - it4 == 1)

}


void Test_Iterator_Adaptor()
{
	typedef AdaptorIterator<char *> TSpec;
	typedef Iter<char *, TSpec> TIterator;

	Test_Iter<TSpec>();

	char arr1[] = "abc";
	TIterator it1 = begin(arr1);
	char * ptr1 = it1;
	SEQAN_TASSERT(*ptr1 == 'a');
}

void Test_Iterator_Position()
{
	Test_Iter<PositionIterator>();
}

//////////////////////////////////////////////////////////////////////////////
//class for testing move operations

struct MoveObj
{
	mutable int data_dat;

	MoveObj(int dat = 0): data_dat(dat) {}
	MoveObj(MoveObj const & other_): data_dat(other_.data_dat) 
	{ 
		other_.data_dat = 0; 
	}
	MoveObj const & operator = (MoveObj const & other_) const
	{ 
		data_dat = other_.data_dat; 
		other_.data_dat = 0; 
		return *this;
	}
	~MoveObj() {}
};

//____________________________________________________________________________

void Test_Transport()
{
	MoveObj o1(10);
	MoveObj o2;
	move(o2, o1);
	SEQAN_TASSERT(o1.data_dat == 0);
	SEQAN_TASSERT(o2.data_dat == 10);

	MoveObj const o3;
	move(o3, o2);
	SEQAN_TASSERT(o2.data_dat == 0);
	SEQAN_TASSERT(o3.data_dat == 10);

	MoveObj const o4;
	move(o4, o3);
	SEQAN_TASSERT(o3.data_dat == 0);
	SEQAN_TASSERT(o4.data_dat == 10);

	move(o1, o4);
	SEQAN_TASSERT(o4.data_dat == 0);
	SEQAN_TASSERT(o1.data_dat == 10);
}

//////////////////////////////////////////////////////////////////////////////

void Main_Test_Common();
void Main_Test_Allocator();
void Main_Test_Alphabet();

//////////////////////////////////////////////////////////////////////////////

int main() 
{
	SEQAN_TREPORT("TEST BEGIN")

	Test_Proxy_Iterator();
	Test_Holder();
	Test_Iterator_Basic();

	Test_Iterator_Adaptor();
	Test_Iterator_Position();

	Test_Transport();

	Main_Test_Common();
	Main_Test_Allocator();
	Main_Test_Alphabet();

//	Test_Alterator_Iterator_Converter();

//	debug::verifyCheckpoints("projects/library/seqan/basic/basic_operator.h");
	debug::verifyCheckpoints("projects/library/seqan/basic/basic_transport.h");
	debug::verifyCheckpoints("projects/library/seqan/basic/basic_proxy.h");
	debug::verifyCheckpoints("projects/library/seqan/basic/basic_holder.h");
	debug::verifyCheckpoints("projects/library/seqan/basic/basic_iterator.h");
	debug::verifyCheckpoints("projects/library/seqan/basic/basic_iterator_base.h");
	debug::verifyCheckpoints("projects/library/seqan/basic/basic_iterator_adaptor.h");
	debug::verifyCheckpoints("projects/library/seqan/basic/basic_iterator_position.h");

	SEQAN_TREPORT("TEST END")

	return 0;
}

#include <iostream>
#include <cstdio>

#define SEQAN_TEST
#define SEQAN_NOSRAN

#include <seqan/sequence.h>
#include <seqan/file.h>
#include <seqan/map.h>



using namespace std;
using namespace seqan;

//////////////////////////////////////////////////////////////////////////////
// Simple Sum List Reference Implementation 
//////////////////////////////////////////////////////////////////////////////

namespace seqan
{

struct _Dummy;
struct _DummySumListIterator;

template <unsigned int DIM, typename TValue>
class SumList<DIM, TValue, _Dummy>
{
public:
	typedef SumListValues<DIM, TValue> TValues;
	typedef String<TValues> TString;

	TString data;
};

//////////////////////////////////////////////////////////////////////////////

template <unsigned int DIM, typename TValue, typename TIteratorSpec>
struct Iterator<SumList<DIM, TValue, _Dummy>, TIteratorSpec>
{
	typedef SumList<DIM, TValue, _Dummy> TSumList;
	typedef Iter<TSumList, _DummySumListIterator> Type;
};

//////////////////////////////////////////////////////////////////////////////

template <typename T>
struct Values;

template <unsigned int DIM, typename TValue>
struct Values<SumList<DIM, TValue, _Dummy> >
{
	typedef SumListValues<DIM, TValue> Type;
};

//////////////////////////////////////////////////////////////////////////////


template <unsigned int DIM, typename TValue>
inline typename Size< SumList<DIM, TValue, _Dummy> >::Type
length(SumList<DIM, TValue, _Dummy> & me)
{
	return length(me.data);
}

//////////////////////////////////////////////////////////////////////////////

template <unsigned int DIM, typename TValue>
inline typename Value< SumList<DIM, TValue, _Dummy> >::Type 
getSum(SumList<DIM, TValue, _Dummy> & me,
	   unsigned int dim)
{
	typedef SumListValues<DIM, TValue> TValues;
	typedef String<TValues> TString;

	typedef typename Iterator<TString>::Type TIterator;
	TValue sum = 0;
	
	for (TIterator it = begin(me.data); it != end(me.data); ++it)
	{
		sum += value(it)[dim];
	}

	return sum;
}

//////////////////////////////////////////////////////////////////////////////

template <unsigned int DIM, typename TValue>
inline void
clear(SumList<DIM, TValue, _Dummy> & me)
{
	clear(me.data);
}

//////////////////////////////////////////////////////////////////////////////

template <unsigned int DIM, typename TValue>
inline typename Iterator< SumList<DIM, TValue, _Dummy> >::Type
begin(SumList<DIM, TValue, _Dummy> & me)
{
	typedef SumList<DIM, TValue, _Dummy>  TMe;
	typedef typename Iterator<TMe>::Type TIterator;
	return TIterator(me);
}

//////////////////////////////////////////////////////////////////////////////

template <unsigned int DIM, typename TValue>
inline typename Iterator< SumList<DIM, TValue, _Dummy> >::Type
end(SumList<DIM, TValue, _Dummy> & me)
{
	typedef SumList<DIM, TValue, _Dummy>  TMe;
	typedef typename Iterator<TMe>::Type TIterator;
	return TIterator(me, GoEnd());
}

//////////////////////////////////////////////////////////////////////////////

template <unsigned int DIM, typename TValue, typename TValues>
inline void 
appendValues(SumList<DIM, TValue, _Dummy> & me,
			 TValues const & new_values)
{
	appendValue(me.data, new_values);
}
template <unsigned int DIM, typename TValue, typename TValue2>
inline void 
appendValues(SumList<DIM, TValue, _Dummy> & me,
			 TValue2 const * new_values)
{
	SumListValues<DIM, TValue> vals(new_values);
	appendValues(me, vals);
}

//////////////////////////////////////////////////////////////////////////////

//template <typename TSumList, typename TValue>
//inline void
//searchSumList(Iter<TSumList, _DummySumListIterator> & it,
//			  TValue const & val,
//			  int dim)
//{
//	goBegin(it);
//	while (!atEnd(it) && ((getSum(it, dim) + getValue(it, dim))< val))
//	{
//		goNext(it);
//	}
//}

//////////////////////////////////////////////////////////////////////////////

template <unsigned int DIM, typename TValue>
class Iter<SumList<DIM, TValue, _Dummy>, _DummySumListIterator>
{
public:
	typedef SumList<DIM, TValue, _Dummy> TSumList;
	typedef SumListValues<DIM, TValue> TValues;
	typedef String<TValues> TString;
	typedef typename Iterator<TString>::Type TStringIterator;

	TSumList * container;
	TStringIterator iter;

	Iter()
	{
	}
	Iter(TSumList & sl)
		: container(& sl)
		, iter(begin(sl.data))
	{
	}
	Iter(TSumList & sl, GoEnd)
		: container(& sl)
		, iter(end(sl.data))
	{
	}
	Iter(Iter const & other)
		: container(other.container)
		, iter(other.iter)
	{
	}
	~Iter()
	{
	}
	Iter const & operator = (Iter const & other)
	{
		container = other.container;
		iter = other.iter;
		return *this;
	}
};


//////////////////////////////////////////////////////////////////////////////

template <typename TSumList>
inline void
goNext(Iter<TSumList, _DummySumListIterator> & it)
{
	goNext(it.iter);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSumList>
inline void
goBegin(Iter<TSumList, _DummySumListIterator> & it)
{
	it.iter = begin(it.container->data);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSumList>
inline bool
atEnd(Iter<TSumList, _DummySumListIterator> & it)
{
	return it.iter == end(it.container->data);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSumList>
inline typename Value<TSumList>::Type
getValue(Iter<TSumList, _DummySumListIterator > & it,
		 int dim)
{
	return value(it.iter)[dim];
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSumList>
inline typename Values<TSumList>::Type
getValues(Iter<TSumList, _DummySumListIterator > & it)
{
	return value(it.iter);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSumList>
inline typename Value<TSumList>::Type
getSum(Iter<TSumList, _DummySumListIterator > & it,
	   int dim)
{
	typedef typename Value<TSumList>::Type TValue;
	typedef SumListValues<DIMENSION<TSumList>::VALUE, TValue> TValues;
	typedef String<TValues> TString;

	typedef typename Iterator<TString>::Type TIterator;
	TValue sum = 0;
	
	for (TIterator it2 = begin(it.container->data); it2 != it.iter; ++it2)
	{
		sum += value(it2)[dim];
	}

	return sum;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSumList, typename TValue>
inline void
searchSumList(Iter<TSumList, _DummySumListIterator> & it,
			  TValue const & val,
			  int dim)
{
	goBegin(it);
	TValue sum = 0;
	while (!atEnd(it) && ((sum + getValue(it, dim))<= val)) // SEARCH SEMANTICS
	{
		sum += getValue(it, dim);
		goNext(it);
	}
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSumList, typename TValue2>
inline void 
assignValue(Iter<TSumList, _DummySumListIterator > & it,
			int dim,
			TValue2 val)
{
	value(it.iter)[dim] = val;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSumList, typename TValues>
inline void 
insertValues(Iter<TSumList, _DummySumListIterator > & it,
			 TValues const & vals)
{
	size_t position = it.iter - begin(it.container->data);
	insertValue(it.container->data, it.iter - begin(it.container->data), vals);
	it.iter = begin(it.container->data) + position;
}
template <typename TSumList, typename TValue>
inline void 
insertValues(Iter<TSumList, _DummySumListIterator > & it,
			 TValue const * p_vals)
{
	SumListValues<DIMENSION<TSumList>::VALUE, typename Value<TSumList>::Type > vals(p_vals);
	return insertValues(it, vals);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSumList>
inline void 
removeValues(Iter<TSumList, _DummySumListIterator > & it)
{
	arrayCopyForward(it.iter + 1, end(it.container->data, Standard()), it.iter);
	resize(it.container->data, length(it.container->data) - 1);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSumList>
inline bool
operator == (Iter<TSumList, _DummySumListIterator> const & left,
			 Iter<TSumList, _DummySumListIterator> const & right)
{
	return left.iter == right.iter;
}

template <typename TSumList>
inline bool
operator != (Iter<TSumList, _DummySumListIterator> const & left,
			 Iter<TSumList, _DummySumListIterator> const & right)
{
	return left.iter != right.iter;
}


//////////////////////////////////////////////////////////////////////////////

} //namespace seqan

//////////////////////////////////////////////////////////////////////////////
// Tests for Mini Sum List
//////////////////////////////////////////////////////////////////////////////

void Test_MiniSumList_Entry()
{
	typedef _MiniListEntry<size_t> TEntry;

	unsigned char buf [1 + sizeof(size_t)];

	unsigned char & buf2char = * buf;
	unsigned short & buf2short = * reinterpret_cast<unsigned short *>(buf);
	unsigned int & buf2int = * reinterpret_cast<unsigned int *>(buf);
	size_t & buf2full = * reinterpret_cast<size_t *>(buf + 1);

	TEntry & entr = * (reinterpret_cast<TEntry *>(buf));

	entr.assignValue(0x12);
	SEQAN_TASSERT(entr.getValue() == 0x12);
	SEQAN_TASSERT(entr.size() == 1);
	SEQAN_TASSERT(buf2char == 0x12 << 2);

	entr.assignValue(0x1234);
	SEQAN_TASSERT(entr.getValue() == 0x1234);
	SEQAN_TASSERT(entr.size() == 2);
	SEQAN_TASSERT((buf2short & 0xfffc) == (0x1234 << 2));

	entr.assignValue(0x1234567);
	SEQAN_TASSERT(entr.getValue() == 0x1234567);
	SEQAN_TASSERT(entr.size() == 4);
	SEQAN_TASSERT((buf2int & 0xfffffffc) == (0x1234567 << 2));

	entr.assignValue(0xfedcba98);
	SEQAN_TASSERT(entr.getValue() == 0xfedcba98);
	SEQAN_TASSERT(entr.size() == 1 + sizeof(size_t));
	SEQAN_TASSERT(buf2full == 0xfedcba98);

}

//////////////////////////////////////////////////////////////////////////////

template <int DIM>
void Test_MiniSumList()
{
	typedef SumList<DIM, size_t, MiniSumList< > > TSumList;
	typedef	SumListValues<DIM, size_t> TValues;
	typedef typename Iterator<TSumList>::Type TIterator;

	TSumList sumlist1;

	TValues vals_sum;

	while (true)
	{
		TValues vals;
		for (int i = 0; i < DIM; ++i)
		{
			vals[i] = rand() % 0x1000;
		}
		if (!appendValues(sumlist1, vals)) break;
		vals_sum += vals;
	}

	SEQAN_TASSERT(sumlist1.data_sum == vals_sum)

	TSumList sumlist2;
	splitSumList(sumlist1, sumlist2);
	TValues vals_sum2 = sumlist1.data_sum;
	vals_sum2 += sumlist2.data_sum;
	SEQAN_TASSERT(vals_sum2 == vals_sum)

}

//////////////////////////////////////////////////////////////////////////////

void Test_MiniSumList2()
{
	typedef SumList<3, int, MiniSumList< > > TSumList;
	typedef Iterator<TSumList>::Type TIterator;
	typedef	SumListValues<3, int> TValues;

	int const VALUES [][3] = 
	{
		{1, 2, 3},
		{2, 3, 1},
		{0, 1, 2},
		{3, 1, 1},
		{1, 2, 1},
	};

	TSumList sl;

	for (unsigned int i = 0; i < sizeof(VALUES) / sizeof(VALUES[0]); ++i)
	{
		appendValues(sl, VALUES[i]);
	}
	SEQAN_TASSERT(length(sl) == 5)

	SEQAN_TASSERT(getSum(sl, 0) == 7)
	SEQAN_TASSERT(getSum(sl, 1) == 9)
	SEQAN_TASSERT(getSum(sl, 2) == 8)

	TValues vals;

	TIterator it(sl);

	searchSumList(it, 2, 0); // SEARCH SEMANTICS
	SEQAN_TASSERT(getValue(it, 0) == 2)
	SEQAN_TASSERT(getValue(it, 1) == 3)
	SEQAN_TASSERT(getValue(it, 2) == 1)
	SEQAN_TASSERT(getSum(it, 0) == 1)
	SEQAN_TASSERT(getSum(it, 1) == 2)
	SEQAN_TASSERT(getSum(it, 2) == 3)

	searchSumList(it, 10, 0);
	SEQAN_TASSERT(getSum(it, 0) == 7)
	SEQAN_TASSERT(getSum(it, 1) == 9)
	SEQAN_TASSERT(getSum(it, 2) == 8)

	SEQAN_TASSERT(it == end(sl))
	SEQAN_TASSERT(it != begin(sl))


	searchSumList(it, 3, 0);// SEARCH SEMANTICS
	SEQAN_TASSERT(getValue(it, 0) == 3)
	SEQAN_TASSERT(getValue(it, 1) == 1)
	SEQAN_TASSERT(getValue(it, 2) == 1)
	SEQAN_TASSERT(getSum(it, 0) == 3)
	SEQAN_TASSERT(getSum(it, 1) == 6)
	SEQAN_TASSERT(getSum(it, 2) == 6)

	assignValue(it, 0, 80000);
	SEQAN_TASSERT(getValue(it, 0) == 80000)
	SEQAN_TASSERT(getSum(sl, 0) == 80004)

	//removeValues
	removeValues(it);
	SEQAN_TASSERT(length(sl) == 4)
	SEQAN_TASSERT(!atEnd(it))
	SEQAN_TASSERT(getValue(it, 0) == 1)
	SEQAN_TASSERT(getValue(it, 1) == 2)
	SEQAN_TASSERT(getValue(it, 2) == 1)
	SEQAN_TASSERT(getSum(sl, 0) == 4)
	SEQAN_TASSERT(getSum(sl, 1) == 8)
	SEQAN_TASSERT(getSum(sl, 2) == 7)

	removeValues(it);
	SEQAN_TASSERT(length(sl) == 3)
	SEQAN_TASSERT(atEnd(it))

	//insertValues
	searchSumList(it, 2, 0);
	SEQAN_TASSERT(getValue(it, 0) == 2)
	SEQAN_TASSERT(getValue(it, 1) == 3)
	SEQAN_TASSERT(getValue(it, 2) == 1)

	insertValues(it, VALUES[3]);
	SEQAN_TASSERT(length(sl) == 4)
	SEQAN_TASSERT(getValue(it, 0) == 3)
	SEQAN_TASSERT(getValue(it, 1) == 1)
	SEQAN_TASSERT(getValue(it, 2) == 1)
	SEQAN_TASSERT(getSum(sl, 0) == 6)
	SEQAN_TASSERT(getSum(sl, 1) == 7)
	SEQAN_TASSERT(getSum(sl, 2) == 7)

	goNext(it);
	SEQAN_TASSERT(getValue(it, 0) == 2)
	SEQAN_TASSERT(getValue(it, 1) == 3)
	SEQAN_TASSERT(getValue(it, 2) == 1)
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//helper functions

template <unsigned int DIM, typename TValue>
inline void
createRandomSumListValues(SumListValues<DIM, TValue> & values)
{
	for (unsigned int i = 0; i < DIM; ++i)
	{
		values[i] = rand() % 1000;
	}
}


//____________________________________________________________________________

template <typename TSumList1, typename TSumList2>
void compareSumLists(TSumList1 & s1,
					 TSumList2 & s2)
{
	typedef typename Iterator<TSumList1>::Type TIterator1;
	typedef typename Iterator<TSumList2>::Type TIterator2;
	typedef typename Values<TSumList2>::Type TValues;

	TIterator1 it1(s1);
	TIterator2 it2(s2);

	int DIM = DIMENSION<TSumList1>::VALUE;

	SEQAN_TASSERT(length(s1) == length(s2))

	for (int i = 0; ; ++i)
	{
		if (atEnd(it1)) break;
		SEQAN_TASSERT(!atEnd(it2))
		SEQAN_TASSERT(getValues(it1) == getValues(it2))

		goNext(it1);
		goNext(it2);
	}
	SEQAN_TASSERT(atEnd(it2))

	for (int i = 0; i < DIM; ++i)
	{
		SEQAN_TASSERT(getSum(s1, i) == getSum(s2, i))
	}
}

//____________________________________________________________________________

template <int DIM, typename TIter1, typename TIter2>
void compareIterators(TIter1 & it1,
					  TIter2 & it2)
{
	for (int i = 0; i < DIM; ++i)
	{
		SEQAN_TASSERT(getSum(it1, i) == getSum(it2, i))
		SEQAN_TASSERT(getValue(it1, i) == getValue(it2, i))
	}
}


//____________________________________________________________________________


//testet, ob die Summen an den Kanten stimmen
template <unsigned int DIM, typename TValue>
void testSkipSumListIntegrity(SumList<DIM, TValue, SkipSumList< > > & ssl)
{
	typedef SumListValues<DIM, TValue> TValues;
	typedef SkiplistPath<TValue, _SkipSumList<DIM> > TPath;
	typedef SkiplistElement<TValue, _SkipSumList<DIM> > TElement;

	TPath path(ssl.map);
	for (int counter = 0; !atEnd(path); ++counter)
	{
		//is the bottom line edges correct?
		SEQAN_TASSERT(path.data_elements[0]->data_next[0].values == path.data_elements[0]->minilist.data_sum)

		//are the other edges correct?
		for (int i = 1; i <= ssl.map.data_height; ++i)
		{
			if (path.data_elements[i] != path.data_elements[0]) break;
			SEQAN_TASSERT(path.sums[i] == path.sums[0])
		}

		goNext(path, ssl.map);
	}
}

//____________________________________________________________________________

template <unsigned int DIM, typename TValue>
void Test_SkipSumListStress()
{
	typedef SumList<DIM, TValue, SkipSumList< > > TSkipSumList;
	typedef SumList<DIM, TValue, _Dummy> TDummySumList;
	typedef SumListValues<DIM, TValue> TValues;

	TSkipSumList ssl;
	TSkipSumList ssl2;
	TDummySumList dsl;

	//appendValues
	for (unsigned int i = 0 ; i < 1000; ++i)
	{
		TValues vals;
		createRandomSumListValues(vals);
		appendValues(ssl, vals);
		appendValues(dsl, vals);
	}
	testSkipSumListIntegrity(ssl);

	compareSumLists(ssl, dsl);


	//iteration
	typedef typename Iterator<TSkipSumList>::Type TSkipSumListIterator;
	typedef typename Iterator<TDummySumList>::Type TDummySumListIterator;

	TSkipSumListIterator sit = begin(ssl);
	TDummySumListIterator dit = begin(dsl);

	for (int i = 0; !atEnd(sit); ++i)
	{
		SEQAN_TASSERT(!atEnd(dit))

		compareIterators<DIM>(sit, dit);

		goNext(sit);
		goNext(dit);
	}
	SEQAN_TASSERT(atEnd(dit))


/*
	//searchSumList => is implicitly done by the assignValue test
	for (int i = 0; i < 1000; ++i)
	{
		int dim = rand() % DIM;
		int val = (rand() + (rand() << 16)) % getSum(ssl, dim);

		searchSumList(sit, val, dim);
		searchSumList(dit, val, dim);

		compareIterators<DIM>(sit, dit);
	}

*/
/*
	//_splitMiniList => is implicitly done by the assignValue test

	typedef SkiplistPath<TValue, _SkipSumList<DIM> > TPath;
	typedef SkiplistElement<TValue, _SkipSumList<DIM> > TElement;
	typedef Pair<TValue, TElement *> TPair;

	for (int i = 0; i < 100; ++i)
	{
		cout << i << "\n";
		int dim = rand() % DIM;
		int find_sum = (rand() + (rand() << 16)) % getSum(ssl, dim);

		searchSumList(sit, find_sum, dim);
		TPath path;
		_skipsumlistFind(ssl.map, TPair(getSum(sit, 0), sit.element), 0, path);
		_splitMiniList(sit, path);
		testSkipSumListIntegrity(ssl);
	}
*/

	cout << "a thousand dots: ";
	//assignValue, insertValues, removeValues
	for (int i = 0; i < 1000; ++i)
	{
		cout << "."; //progression

		//assignValue
		int dim = rand() % DIM;
		int find_sum = (rand() + (rand() << 16)) % getSum(ssl, dim);
		int new_value = rand() | (rand() << 10);

		searchSumList(sit, find_sum, dim);
		searchSumList(dit, find_sum, dim);

		compareIterators<DIM>(sit, dit);

		assignValue(sit, dim, new_value);
		testSkipSumListIntegrity(ssl);

		assignValue(dit, dim, new_value);

		compareIterators<DIM>(sit, dit);
		compareSumLists(ssl, dsl);


		//insertValues
		find_sum = (rand() + (rand() << 16)) % getSum(ssl, dim);
		TValues new_values;
		createRandomSumListValues(new_values);

		searchSumList(sit, find_sum, dim);
		searchSumList(dit, find_sum, dim);

		insertValues(sit, new_values);
		testSkipSumListIntegrity(ssl);

		insertValues(dit, new_values);

		compareIterators<DIM>(sit, dit);
		compareSumLists(ssl, dsl);


		//removeValues
		if (i % 2) 
		{
			dim = rand() % DIM;
			find_sum = (rand() + (rand() << 16)) % getSum(ssl, dim);

			searchSumList(sit, find_sum, dim);
			searchSumList(dit, find_sum, dim);

			removeValues(sit);
			testSkipSumListIntegrity(ssl);

			removeValues(dit);

			compareIterators<DIM>(sit, dit);
			compareSumLists(ssl, dsl);
		}

		ssl2 = ssl;
		testSkipSumListIntegrity(ssl);
		testSkipSumListIntegrity(ssl2);
		compareSumLists(ssl, ssl2);
	}
	cout << "\n";

}

//////////////////////////////////////////////////////////////////////////////

void Main_TestSumlist() 
{
	SEQAN_TREPORT("TEST SUMLIST BEGIN")

	Test_MiniSumList_Entry();
	Test_MiniSumList<5>();
	Test_MiniSumList<1>();
	Test_MiniSumList2();

	Test_SkipSumListStress<3, size_t>();

	SEQAN_TREPORT("TEST SUMLIST END")
}

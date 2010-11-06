#include <iostream>

#define SEQAN_DEBUG
#define SEQAN_TEST

#include "seqan/basic.h"
#include "seqan/sequence.h"


using namespace std;
using namespace seqan;

//Helper class that counts ctors, dtors and copys
struct Test1
{
	static int m_ctor_count;
	static int m_dtor_count;
	static int m_copy_count;
	static int m_move_count;

	Test1(): x(0xfade)
	{
		++m_ctor_count;
	}
	Test1(Test1 const & obj): x(obj.x)
	{
		++m_ctor_count;
		++m_copy_count;
	}
	Test1(Test1 const & obj, Move): x(obj.x)
	{
		++m_move_count;
		obj.x = 0x105e;
	}
	~Test1()
	{
		x = 0xdead;
		++m_dtor_count;
	}
	Test1 & operator = (Test1 const & obj)
	{
		x = obj.x;
		++m_copy_count;
		return *this;
	}
	mutable int x;
};

int Test1::m_ctor_count = 0;
int Test1::m_dtor_count = 0;
int Test1::m_copy_count = 0;
int Test1::m_move_count = 0;


inline void 
move(Test1 & target, Test1 & source)
{
	++Test1::m_move_count;
	target.x = source.x;
	source.x = 0x105e;
}
inline void 
move(Test1 const & target, Test1 & source)
{
	++Test1::m_move_count;
	target.x = source.x;
	source.x = 0x105e;
}
inline void 
move(Test1 & target, Test1 const & source)
{
	++Test1::m_move_count;
	target.x = source.x;
	source.x = 0x105e;
}
inline void 
move(Test1 const & target, Test1 const & source)
{
	++Test1::m_move_count;
	target.x = source.x;
	source.x = 0x105e;
}

//////////////////////////////////////////////////////////////////////////////
//Test value array function for a class type

void TestAlphabetInterface()
{
	{
		Test1 a; //1 ctor
		a.x = 0xbeef;

		char c_buf1[200 * sizeof(Test1)];
		Test1 * a_buf1 = (Test1 *) c_buf1;

		char c_buf2[200 * sizeof(Test1)];
		Test1 * a_buf2 = (Test1 *) c_buf2;

		arrayConstruct(a_buf1, a_buf1 + 100); //100 ctor
		for (int i=0; i < 100; ++i)
		{
			SEQAN_TASSERT(a_buf1[i].x == 0xfade)
		}

		arrayConstruct(a_buf2, a_buf2 + 100, a); //100 ctor, 100 copy
		for (int i=0; i < 100; ++i)
		{
			SEQAN_TASSERT(a_buf2[i].x == 0xbeef)
		}

		arrayConstruct(a_buf1, a_buf1+ 100); //100 ctor
		arrayDestruct(a_buf1, a_buf1 + 100); //100 dtor
		for (int i=0; i < 100; ++i)
		{
			SEQAN_TASSERT(a_buf1[i].x == 0xdead)
		}

		arrayConstructCopy(a_buf1, a_buf1 + 100, a_buf2); //100 ctor, 100 copy
		arrayDestruct(a_buf1, a_buf1 + 100); //100 dtor
		for (int i=0; i < 100; ++i)
		{
			SEQAN_TASSERT(a_buf2[i].x == 0xdead)
		}

		arrayFill(a_buf1, a_buf1 + 100, a); //100 copy
		for (int i=0; i < 100; ++i)
		{
			SEQAN_TASSERT(a_buf1[i].x == 0xbeef)
		}

		arrayCopyForward(a_buf1, a_buf1 + 100, a_buf2); //100 copy
		
		for (int i=0; i < 100; ++i) a_buf1[i].x = i;

		arrayCopy(a_buf1, a_buf1 + 50, a_buf1 + 20); //50 copy
		for (int i=0; i < 50; ++i)
		{
			SEQAN_TASSERT(a_buf1[i+20].x == i)
		}

		arrayCopy(a_buf1 + 80, a_buf1 + 100, a_buf1 + 75); //20 copy
		for (int i=80; i < 100; ++i)
		{
			SEQAN_TASSERT(a_buf1[i-5].x == i)
		}

		for (int i=0; i < 100; ++i) a_buf1[i].x = i;

		arrayClearSpace(a_buf1, 100, 50, 70); //20 ctor, 70 dtor, 50 copy
		for (int i=50; i < 100; ++i)
		{
			SEQAN_TASSERT(a_buf1[i+20].x == i)
		}
		arrayConstruct(a_buf1, a_buf1 + 70, a); //70 ctor, 70 copy

		arrayClearSpace(a_buf1, 120, 70, 50); //70 dtor, 50 copy
		for (int i=50; i < 100; ++i)
		{
			SEQAN_TASSERT(a_buf1[i].x == i)
		}

		arrayConstruct(a_buf1, a_buf1 + 50, a); //70 ctor, 70 copy

		for (int i=0; i < 100; ++i) a_buf1[i].x = i;

		arrayClearSpace(a_buf1, 100, 90, 110); //10 ctor, 90 dtor, 10 copy
		for (int i=110; i < 120; ++i)
		{
			SEQAN_TASSERT(a_buf1[i].x == i-20)
		}

		arrayDestruct(a_buf1 + 110, a_buf1 + 120);	//10 dtor
		arrayDestruct(a_buf2, a_buf2 + 100); //100 dtor


		arrayConstruct(a_buf2, a_buf2 + 23); //23 ctor
		arrayConstructMove(a_buf2, a_buf2 + 23, a_buf1); // 23 move 
		for (int i = 0; i < 23; ++i)
		{
			SEQAN_TASSERT(a_buf1[i].x == 0xfade);
			SEQAN_TASSERT(a_buf2[i].x == 0x105e);
		}

		arrayMove(a_buf1, a_buf1 + 23, a_buf1 + 5); // 23 move 
		for (int i = 0; i < 23; ++i)
		{
			SEQAN_TASSERT(a_buf1[i + 5].x == 0xfade);
		}
		for (int i = 0; i < 5; ++i)
		{
			SEQAN_TASSERT(a_buf1[i].x == 0x105e);
		}

		arrayMove(a_buf1 + 5, a_buf1 + 28, a_buf1); // 23 move 

		arrayDestruct(a_buf1, a_buf1 + 23); //23 dtor

		//1 dtor for a
	}

	SEQAN_TASSERT(Test1::m_ctor_count == 574)
	SEQAN_TASSERT(Test1::m_dtor_count == 574)
	SEQAN_TASSERT(Test1::m_copy_count == 700)
	SEQAN_TASSERT(Test1::m_move_count == 69)


	SEQAN_TASSERT(gapValue<char>() == '-');
	SEQAN_TASSERT(gapValue<int>() == int());
}

//////////////////////////////////////////////////////////////////////////////
//Test value array functions for some types

template <typename T, typename _T>
void TestArrayFunctions(_T const _val1, _T const _val2)
{
	T val1 = (T)_val1;
	T val2 = (T)_val2;

	T a = val1;

	T a_buf1[200];

	T a_buf2[200];

	arrayConstruct(a_buf1, a_buf1 + 100); //nothing happens

	arrayConstruct(a_buf2, a_buf2 + 100, a);
	for (int i=0; i < 100; ++i)
	{
		SEQAN_TASSERT(a_buf2[i] == val1)
	}

	arrayDestruct(a_buf1, a_buf1 + 100); //nothing happens

	arrayConstructCopy(a_buf2, a_buf2 + 100, a_buf1); 
	for (int i=0; i < 100; ++i)
	{
		SEQAN_TASSERT(a_buf1[i] == val1)
	}

	a = val2;

	arrayFill(a_buf1, a_buf1 + 100, a); 
	for (int i=0; i < 100; ++i)
	{
		SEQAN_TASSERT(a_buf1[i] == val2)
	}

	arrayCopyForward(a_buf1, a_buf1 + 100, a_buf2);
	for (int i=0; i < 100; ++i)
	{
		SEQAN_TASSERT(a_buf2[i] == val2)
	}
	
	for (int i=0; i < 100; ++i) a_buf1[i] = (T) i;

	arrayCopy(a_buf1, a_buf1 + 50, a_buf1 + 20); 
	for (int i=0; i < 50; ++i)
	{
		SEQAN_TASSERT(a_buf1[i+20] == (T)i)
	}

	arrayCopy(a_buf1 + 80, a_buf1 + 100, a_buf1 + 75); 
	for (int i=80; i < 100; ++i)
	{
		SEQAN_TASSERT(a_buf1[i-5] == (T)i)
	}

	for (int i=0; i < 100; ++i) a_buf1[i] = (T) i;

	arrayClearSpace(a_buf1, 100, 50, 70); 
	for (int i=50; i < 100; ++i)
	{
		SEQAN_TASSERT(a_buf1[i+20] == (T)i)
	}

	arrayClearSpace(a_buf1, 120, 70, 50); 
	for (int i=50; i < 100; ++i)
	{
		SEQAN_TASSERT(a_buf1[i] == (T)i)
	}

}

//////////////////////////////////////////////////////////////////////////////
//Test SimpleType instances

template <typename T>
void TestSimpleType()
{
	T a;
	T b(a);
	b = a;

	typename Value<T>::Type val;
	val = a;
	a = val;

	T a_buf1[200];
	T a_buf2[200];
	T a_buf3[200];

	arrayConstruct(a_buf1, a_buf1 + 100);
	arrayConstruct(a_buf1, a_buf1 + 100, a);
	arrayDestruct(a_buf1, a_buf1 + 100);
	arrayConstructCopy(a_buf2, a_buf2 + 100, a_buf1); 
	arrayFill(a_buf1, a_buf1 + 100, a); 
	arrayCopyForward(a_buf1, a_buf1 + 100, a_buf2);
	arrayCopy(a_buf1, a_buf1 + 50, a_buf1 + 20); 
	arrayCopy(a_buf1 + 80, a_buf1 + 100, a_buf1 + 75); 
	arrayMoveForward(a_buf1, a_buf1 + 10, a_buf3); 
	arrayMoveBackward(a_buf3, a_buf3 + 10, a_buf1); 
	arrayClearSpace(a_buf1, 100, 50, 70); 
	arrayClearSpace(a_buf1, 120, 70, 50); 
}

//////////////////////////////////////////////////////////////////////////////

template <typename TTarget, typename TSource>
void TestConversion()
{
	static TSource source;
	TTarget target = source;

	SEQAN_ASSERT(target == source)
	SEQAN_ASSERT(source == target)
	SEQAN_ASSERT(!(target != source))
	SEQAN_ASSERT(!(source != target))
	SEQAN_ASSERT(!(target < source))
	SEQAN_ASSERT(!(source < target))
	SEQAN_ASSERT(target <= source)
	SEQAN_ASSERT(source <= target)
	SEQAN_ASSERT(!(target > source))
	SEQAN_ASSERT(!(source > target))
	SEQAN_ASSERT(target >= source)
	SEQAN_ASSERT(source >= target)

	TSource const c_source = TSource();
	target = c_source;

	target = TSource();
	assign(target, source);

	TSource a_source_1[200];
	TTarget a_target[200];

	arrayCopyForward(a_source_1, a_source_1 + 100, a_target);
	arrayCopy(a_source_1, a_source_1 + 50, a_target + 20); 
	arrayCopyBackward(a_source_1, a_source_1 + 100, a_target);
	arrayMoveForward(a_source_1, a_source_1 + 100, a_target);
	arrayMove(a_source_1, a_source_1 + 50, a_target + 20); 
	arrayMoveBackward(a_source_1, a_source_1 + 100, a_target);
}

void TestSimpleTypeConversions()
{
	TestConversion<Ascii, Dna>();
	TestConversion<Ascii, Dna5>();
	TestConversion<Ascii, Iupac>();
	TestConversion<Ascii, AminoAcid>();

	TestConversion<Dna, Ascii>();
	TestConversion<Dna, Byte>();
	TestConversion<Dna, Unicode>();
	TestConversion<Dna, int>();
	TestConversion<Dna, Dna5>();
	TestConversion<Dna, Iupac>();

	TestConversion<Dna5, Ascii>();
	TestConversion<Dna5, Byte>();
	TestConversion<Dna5, Unicode>();
	TestConversion<Dna5, Dna>();
	TestConversion<Dna5, Iupac>();

	TestConversion<Iupac, Ascii>();
	TestConversion<Iupac, Byte>();
	TestConversion<Iupac, Unicode>();
	TestConversion<Iupac, Dna>();
	TestConversion<Iupac, Dna5>();

	TestConversion<AminoAcid, Ascii>();
	TestConversion<AminoAcid, Byte>();
	TestConversion<AminoAcid, Unicode>();

	typedef SimpleType<int, void> ST;

	TestConversion<int, ST>();
	TestConversion<unsigned int, ST>();
	TestConversion<short, ST>();
	TestConversion<unsigned short, ST>();
	TestConversion<char, ST>();
	TestConversion<signed char, ST>();
	TestConversion<unsigned char, ST>();
	TestConversion<ST, ST>();

}

//////////////////////////////////////////////////////////////////////////////
//Test infimum / supremum values

template <typename T>
void TestExtremeValuesSigned()
{
	long double minVal = -1;
	for(int e = 1; e < BitsPerValue<T>::VALUE; ++e)
		minVal = 2*minVal;

	long double maxVal = -minVal - 1;

/*
	::std::cout << ::std::endl << "Max/Min of " << typeid(T).name() << ::std::endl;
	::std::cout << maxVal << " == " << SupremumValue<T>::VALUE << "(" << (double)SupremumValue<T>::VALUE << ")  " << supremumValue<T>() << ::std::endl;
	::std::cout << minVal << " == " << InfimumValue<T>::VALUE << "(" << (double)InfimumValue<T>::VALUE << ")  " << infimumValue<T>() << ::std::endl;
*/

	bool isSigned = TYPECMP< typename _MakeSigned<T>::Type, T >::VALUE;
	SEQAN_ASSERT(isSigned);

	SEQAN_ASSERT(supremumValue<T>() == SupremumValue<T>::VALUE);
	SEQAN_ASSERT(infimumValue<T>()  == InfimumValue<T>::VALUE);

	long double maxDelta = maxVal - SupremumValue<T>::VALUE;
	long double minDelta = minVal - (long double)InfimumValue<T>::VALUE;
	SEQAN_ASSERT(maxDelta <= maxVal/1000);
	SEQAN_ASSERT(-maxVal/1000 <= maxDelta);
	SEQAN_ASSERT(minDelta <= maxVal/1000);
	SEQAN_ASSERT(-maxVal/1000 <= minDelta);
}

template <typename T>
void TestExtremeValuesUnsigned()
{
	long double maxVal = 1;
	for(int e = 0; e < BitsPerValue<T>::VALUE; ++e)
		maxVal = 2*maxVal;
	maxVal = maxVal - 1;

/*
	::std::cout << ::std::endl << "Max/Min of " << typeid(T).name() << ::std::endl;
	::std::cout << maxVal << " == " << SupremumValue<T>::VALUE << "(" << (double)SupremumValue<T>::VALUE << ")  " << supremumValue<T>() << ::std::endl;
	::std::cout << 0 << " == " << InfimumValue<T>::VALUE << "(" << (double)InfimumValue<T>::VALUE << ")  " << infimumValue<T>() << ::std::endl;
*/

	bool isUnsigned = TYPECMP< typename _MakeUnsigned<T>::Type, T >::VALUE;
	SEQAN_ASSERT(isUnsigned);

	SEQAN_ASSERT(supremumValue<T>() == SupremumValue<T>::VALUE);
	SEQAN_ASSERT(infimumValue<T>()  == InfimumValue<T>::VALUE);

	long double maxDelta = maxVal - SupremumValue<T>::VALUE;
	SEQAN_ASSERT(maxDelta <= maxVal/1000);
	SEQAN_ASSERT(-maxVal/1000 <= maxDelta);
	SEQAN_ASSERT(0 == InfimumValue<T>::VALUE);
}

void TestExtremeValues()
{
	TestExtremeValuesSigned<signed char>();
	TestExtremeValuesSigned<signed short>();
	TestExtremeValuesSigned<signed int>();
	TestExtremeValuesSigned<signed long>();
	TestExtremeValuesUnsigned<unsigned char>();
	TestExtremeValuesUnsigned<unsigned short>();
	TestExtremeValuesUnsigned<unsigned int>();
	TestExtremeValuesUnsigned<unsigned long>();
	TestExtremeValuesSigned<__int64>();
/*	TestExtremeValues<float>();
	TestExtremeValues<double>();
	TestExtremeValues<long double>();*/
}

//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

void Main_Test_Alphabet() 
{   
	SEQAN_TREPORT("TEST ALPHABET BEGIN")

	TestAlphabetInterface();

	TestSimpleType<Dna>();
	TestSimpleType<Dna5>();
	TestSimpleType<Iupac>();
	TestSimpleType<AminoAcid>();
	TestSimpleTypeConversions();
	TestExtremeValues();
  
	TestSimpleType<bool>();
	TestArrayFunctions<char>(0xde, 0xad);
	TestArrayFunctions<signed char>(0xde, 0xad);
	TestArrayFunctions<unsigned char>(0xde, 0xad);
	TestArrayFunctions<short>(0xdead, 0xbeef);
	TestArrayFunctions<unsigned short>(0xdead, 0xbeef);
	TestArrayFunctions<int>(0xdead, 0xbeef);
	TestArrayFunctions<unsigned int>(0xdead, 0xbeef);
	TestArrayFunctions<float>(3.1, 1.2);
	TestArrayFunctions<double>(3.1, 1.2);
	TestArrayFunctions<long double>(3.1, 1.2);

	debug::verifyCheckpoints("projects/library/seqan/basic/basic_alphabet_interface.h");
	debug::verifyCheckpoints("projects/library/seqan/basic/basic_alphabet_simple.h");
	debug::verifyCheckpoints("projects/library/seqan/basic/basic_alphabet_trait_basic.h");

	SEQAN_TREPORT("TEST ALPHABET END")
}

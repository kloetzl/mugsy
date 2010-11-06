#include <iostream>
#include <fstream>
#include <typeinfo>
#include <time.h>

#define SEQAN_DEBUG
#define SEQAN_TEST

#include "seqan/sequence.h"

using namespace std;
using namespace seqan;

//////////////////////////////////////////////////////////////////////////////
//test sequence default interface:
//non-container objects are treated like containers of length 1

struct DummyClass
{
};

void Test_Sequence_Interface()
{
//* ???Anti Default Sequences
	DummyClass const c = DummyClass();
	DummyClass d;
	DummyClass e;

	SEQAN_TASSERT(id(c) != id(e));
	SEQAN_TASSERT(id(d) != id(e));
	SEQAN_TASSERT(id(c) != 0);
	SEQAN_TASSERT(id(d) != 0);

	SEQAN_TASSERT(begin(c, Standard()) == & c);		//begin
	SEQAN_TASSERT(begin(d, Standard()) == & d);

	SEQAN_TASSERT(end(c, Standard()) == & c + 1);	//end
	SEQAN_TASSERT(end(d, Standard()) == & d + 1);

	SEQAN_TASSERT(beginPosition(c) == 0);			//beginPosition
	SEQAN_TASSERT(beginPosition(d) == 0);			
	SEQAN_TASSERT(endPosition(c) == (size_t)(end(c) - begin(c))); //endPosition
	SEQAN_TASSERT(endPosition(d) == (size_t)(end(d) - begin(d)));

	SEQAN_TASSERT(iter(c, endPosition(c)) == end(c));	//iter
	SEQAN_TASSERT(iter(d, endPosition(d)) == end(d));

	SEQAN_TASSERT(& getValue(c, 10) == & c);
	SEQAN_TASSERT(& getValue(d, 10) == & d);

	SEQAN_TASSERT(length(c) == 1);
	SEQAN_TASSERT(length(d) == 1);

	SEQAN_TASSERT(capacity(c) == 1);
	SEQAN_TASSERT(capacity(d) == 1);

	SEQAN_TASSERT(empty(c) == false);
	SEQAN_TASSERT(empty(d) == false);

	Iterator<DummyClass, Standard>::Type it_1;
	it_1 = begin(d, Standard());

	goNext(it_1);
	goPrevious(it_1);
	SEQAN_TASSERT(it_1 == begin(d, Standard()))
//*/

//____________________________________________________________________________
//test interfaces for assignment functions:
//semantics are tested somewhere else

	char str[100] = "hallo";

	assign(str, str);
	assign(str, infix(str, 0, 5));
	assign(infix(str, 0, 5), str);
	assign(infix(str, 0, 5), infix(str, 0, 5));
	SEQAN_TASSERT(isEqual(str, "hallo"));
	
	assign(str, str, 5);
	assign(str, infix(str, 0, 5), 5);
	assign(infix(str, 0, 5), str, 5);
	assign(infix(str, 0, 5), infix(str, 0, 5), 5);

	append(str, str);
	append(str, infix(str, 0, 5));
	append(infix(str, 0, 5), str);
	append(infix(str, 0, 5), infix(str, 0, 5));
	
	append(str, str, 5);
	append(str, infix(str, 0, 5), 5);
	append(infix(str, 0, 5), str, 5);
	append(infix(str, 0, 5), infix(str, 0, 5), 5);

	replace(str, 2, 3, str);
	replace(str, 2, 3, infix(str, 0, 5));
	replace(infix(str, 0, 5), 2, 3, str); 
	replace(infix(str, 0, 5), 2, 3, infix(str, 0, 5)); 
	
	replace(str, 2, 3, str, 5);
	replace(str, 2, 3, infix(str, 0, 5), 5);
	replace(infix(str, 0, 5), 2, 3, str, 5); 
	replace(infix(str, 0, 5), 2, 3, infix(str, 0, 5), 5); 

//____________________________________________________________________________

	String<char, Alloc<> > str2;

	Size<String<char, Alloc<> > >::Type cap = reserve(str2, 200);
	SEQAN_TASSERT(cap <= capacity(str2));
	SEQAN_TASSERT(cap == 200);

	cap = reserve(str2, 400, Insist());
	SEQAN_TASSERT(cap == 400);

	cap = reserve(str2, 1000, Limit());
	SEQAN_TASSERT(cap == capacity(str2));


	Size<String<char, Alloc<> > >::Type len = resize(str2, 100);
	SEQAN_TASSERT(len <= length(str2));
	SEQAN_TASSERT(len <= 100);

	len = fill(str2, 150, 'C');
	SEQAN_TASSERT(len <= length(str2));
	SEQAN_TASSERT(len <= 150);


	len = resizeSpace(str2, 100, 50, 100);
	SEQAN_TASSERT(len == 100);
	SEQAN_TASSERT(length(str2) == 200);

	len = resizeSpace(str2, 100, 150, 200, 220);
	SEQAN_TASSERT(len == 70);
	SEQAN_TASSERT(length(str2) == 220);

//____________________________________________________________________________
//test interfaces for appending different string types

	assign(str, "test mixed append");
	String<char> str3;
	clear(str3);
	append(str3, str);
	SEQAN_TASSERT(isEqual(str3, str));
	SEQAN_TASSERT(isEqual(str, str3));

	const char *str4 = "test const char*";
	clear(str3);
	append(str3, str4);
	SEQAN_TASSERT(isEqual(str4, str3));
	SEQAN_TASSERT(isEqual(str3, str4));

	clear(str3);
	append(str3, "test");
	SEQAN_TASSERT(isEqual("test", str3));
	SEQAN_TASSERT(isEqual(str3, "test"));
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//Test semantics of assignment functions like assign, append, ...

template <typename TExpand, typename TMe>
void Test_String_Base_Assignments(TMe & str)
{
	assign(str, "goldfishgoldfishgoldfishgoldfishgoldfish", TExpand());
	SEQAN_TASSERT(isEqual(str, "goldfishgoldfishgoldfishgoldfishgoldfish"));

	assign(str, "x", TExpand());
	SEQAN_TASSERT(isEqual(str,"x"));

	assign(str, "goldfishgoldfishgoldfishgoldfishgoldfish", 8, TExpand());
	SEQAN_TASSERT(isEqual(str,"goldfish"));

	assign(str, "x", 8, TExpand());
	SEQAN_TASSERT(isEqual(str,"x"));

	//test append
	assign(str, "goldfish");

	append(str, "goldfishgoldfishgoldfishgoldfish", TExpand());
	SEQAN_TASSERT(isEqual(str,"goldfishgoldfishgoldfishgoldfishgoldfish"));

	assign(str, "goldfish");
	append(str, "", TExpand());
	SEQAN_TASSERT(isEqual(str,"goldfish"));

	append(str, "goldfish", 4, TExpand());
	SEQAN_TASSERT(isEqual(str,"gold"));

	append(str, "goldfish", 8, TExpand());
	SEQAN_TASSERT(isEqual(str,"goldgold"));
	
	append(str, "goldfish", 16, TExpand());
	SEQAN_TASSERT(isEqual(str,"goldgoldgoldfish"));

	//test replace

	str = "goldfish";
	replace(str, 2, 4, "ndel chips and ", TExpand());
	SEQAN_TASSERT(str == "gondel chips and fish");

	replace(str, 7, 16, "is", TExpand());
	SEQAN_TASSERT(str == "gondel is fish");

	replace(str, 6, 10, "", TExpand());
	SEQAN_TASSERT(str == "gondelfish");

	replace(str, 2, 2, "ld in my ma", TExpand());
	SEQAN_TASSERT(str == "gold in my mandelfish");

	replace(str, (size_t) 0, 0, "here is ", TExpand());
	SEQAN_TASSERT(str == "here is gold in my mandelfish");

	replace(str, length(str), length(str), " and silver", TExpand());
	SEQAN_TASSERT(str == "here is gold in my mandelfish and silver");

	replace(str, 8, length(str), "nothing", TExpand());
	SEQAN_TASSERT(str == "here is nothing");


	assignValue(str, 1, 'a');		//assignValue
	SEQAN_TASSERT(str == "hare is nothing");

	moveValue(str, 1, 'e');			//moveValue
	SEQAN_TASSERT(str == "here is nothing");

	appendValue(str, '!');			//appendValue
	SEQAN_TASSERT(str == "here is nothing!");

	insertValue(str, 7, 't');		//insertValue
	SEQAN_TASSERT(str == "here ist nothing!");

	erase(str, 7);					//erase
	SEQAN_TASSERT(str == "here is nothing!");
	erase(str, 1, 3);
	SEQAN_TASSERT(str == "he is nothing!");

	TMe str2 = "another string";
	move(str, str2);				//move
	SEQAN_TASSERT(str == "another string");
}

void Test_String_Base()
{
	String<char> str1("hello");
	String<char> const str2("HELLO");

	SEQAN_TASSERT(getValue(str1, 0) == 'h');
	SEQAN_TASSERT(getValue(str1, 1) == 'e');
	SEQAN_TASSERT(getValue(str2, 1) == 'E');

	clear(str1);
	SEQAN_TASSERT(length(str1) == 0);

//____________________________________________________________________________
//note: the following tests destroy the same items several times
//note: most of the tests only check syntax not semantic.
//semantic should be tested via the calling assignment functions
/*
	str1 = "hello";
	_clearSpace(str1, 200, Exact());
	SEQAN_TASSERT(length(str1) == 200);
	SEQAN_TASSERT(capacity(str1) == 200);

	_clearSpace(str1, 200, Limit());
	_clearSpace(str1, 300, Generous());
	_clearSpace(str1, 200, Insist());

	_clearSpace(str1, 300, 100, Exact());
	SEQAN_TASSERT(length(str1) == 100);
	_clearSpace(str1, 300, 10000, Limit());
	_clearSpace(str1, 400, 100, Generous());
	_clearSpace(str1, 300, 100, Insist());

	_clearSpace(str1, 300, 10, 20, Exact());
	_clearSpace(str1, 300, 10, 20, Limit());
	_clearSpace(str1, 300, 10, 20, Generous());
	_clearSpace(str1, 100, 10, 200, Insist());

	_clearSpace(str1, 300, begin(str1), end(str1), Exact());
	_clearSpace(str1, 300, begin(str1), end(str1), Limit());
	_clearSpace(str1, 300, begin(str1), end(str1), Generous());
	_clearSpace(str1, 100, begin(str1), end(str1), Insist());

	_clearSpace(str1, 300, 10, 20, 100, Exact());
	_clearSpace(str1, 300, 10, 20, 100, Limit());
	_clearSpace(str1, 300, 10, 20, 100, Generous());
	_clearSpace(str1, 100, 10, 50, 100, Insist());

	_clearSpace(str1, 300, 10, 20, 20000, Exact());
	_clearSpace(str1, 300, 10, 20, 20000, Limit());
	_clearSpace(str1, 300, 10, 20, 20000, Generous());
	_clearSpace(str1, 100, 10, 50, 20000, Insist());
*/

/*
	_clearSpace(str1, 300, begin(str1), end(str1), 100, Exact());
	_clearSpace(str1, 300, begin(str1), end(str1), 100, Limit());
	_clearSpace(str1, 300, begin(str1), end(str1), 100, Generous());
	_clearSpace(str1, 100, begin(str1), end(str1), 100, Insist());
*/

//____________________________________________________________________________
// test assignment functions

	Test_String_Base_Assignments<Exact>(str1);
	Test_String_Base_Assignments<Generous>(str1);
	reserve(str1, 10000);
	Test_String_Base_Assignments<Insist>(str1);

	str1 += str2;
	str1 += str1;

//* ???Anti Default Sequences
	str1 += 'x';
//*/
//____________________________________________________________________________

	String<char> str3;
	resize(str3, 0);
	SEQAN_TASSERT(length(str3) == 0);

	resize(str3, 20);
	SEQAN_TASSERT(length(str3) == 20);
	resize(str3, 10);
	SEQAN_TASSERT(length(str3) == 10);
	resize(str3, 20);
	SEQAN_TASSERT(length(str3) == 20);
	resize(str3, 200);
	SEQAN_TASSERT(length(str3) == 200);

	fill(str3, 100, 'x');
	SEQAN_TASSERT(length(str3) == 100);
	fill(str3, 200, 'x');
	SEQAN_TASSERT(length(str3) == 200);
	fill(str3, 400, 'y');
	SEQAN_TASSERT(length(str3) == 400);

//____________________________________________________________________________

	str3 = "abc";
	SEQAN_TASSERT(str3 < "abcd");
}

//////////////////////////////////////////////////////////////////////////////
//test some basic string features

template <typename TMe>
void TestStringBasics()
{
	//test default ctor
	TMe str1;
	SEQAN_TASSERT(str1 == "");
	SEQAN_TASSERT(length(str1) == 0);
	SEQAN_TASSERT(length(str1) <= capacity(str1));

	//test assignment ctor (with length == 0)
	TMe str2 = str1;
	SEQAN_TASSERT(str1 == str2);
	SEQAN_TASSERT(length(str1) == 0);
	SEQAN_TASSERT(length(str1) <= capacity(str1));

	//test assignment with char const []
	str1 = "hamster";
	SEQAN_TASSERT(str1 == "hamster");
	SEQAN_TASSERT(length(str1) == 7);
	SEQAN_TASSERT(length(str1) <= capacity(str1));

	//test assignment with char *
	char * s1 = (char *) "goldfish";
	str1 = s1;
	SEQAN_TASSERT(str1 == "goldfish");
	SEQAN_TASSERT(length(str1) == 8);
	SEQAN_TASSERT(length(str1) <= capacity(str1));

	//test assignment ctor (with length > 0)
	TMe str3 = str1;
	SEQAN_TASSERT(str3 == str1);
	SEQAN_TASSERT(length(str3) == 8);
	SEQAN_TASSERT(length(str3) <= capacity(str3));

	//test independency
	{
		TMe str4 = "hamster";
		str3 = str4;
		str4 = "...";
	}
	SEQAN_TASSERT(str3 == "hamster");
	SEQAN_TASSERT(length(str3) == 7);
	SEQAN_TASSERT(length(str3) <= capacity(str3));

	typename Size<TMe>::Type len = length(str3);
	SEQAN_TASSERT(len == 7);

	//test begin and end
	SEQAN_TASSERT(end(str3) == begin(str3) + 7);
	typename Iterator<TMe, Rooted>::Type str3_begin = begin(str3, Rooted());
	*str3_begin = 'X';
	SEQAN_TASSERT(str3 == "Xamster");

	//test at
	value(str3, 1) = 'Y';
	str3[2] = 'Z';
	int i1 = 3;
	value(str3, i1) = 'A';
	i1 = 4;
	str3[i1] = 'B';
	SEQAN_TASSERT(str3 == "XYZABer");
	SEQAN_TASSERT(getValue(str3, 5) == 'e');
	SEQAN_TASSERT(str3[6] == 'r');

	//test clear and empty
	SEQAN_TASSERT(!empty(str3));
	clear(str3);
	SEQAN_TASSERT(empty(str3));
	SEQAN_TASSERT(begin(str3) == end(str3));
	SEQAN_TASSERT(length(str3) == 0);
}

//////////////////////////////////////////////////////////////////////////////
//test some basic string features for strings that can change length
//note: capacity of TMe strings must be >= 200

template <typename TMe>
void TestStringResize()
{
	TMe str1;
	resize(str1, 50);
	SEQAN_TASSERT(length(str1) == 50);

	fill(str1, 100, 3);
	SEQAN_TASSERT(length(str1) == 100);
	SEQAN_TASSERT(getValue(str1, 51) == 3);
	SEQAN_TASSERT(getValue(str1, 99) == 3);

}

//////////////////////////////////////////////////////////////////////////////

void Test_String_Alloc()
{
	TestStringBasics<String<char> >();
	TestStringResize<String<char> >();

	String<char, Alloc<> > str1 = "hello";
	SEQAN_TASSERT(str1[1] == 'e');

	String<char, Alloc<> > const str2 = "hello";
	SEQAN_TASSERT(str2[4] == 'o');
	SEQAN_TASSERT(capacity(str2) >= length(str2));

	String<char, Alloc<> > str3;
	move(str3, str1);
	SEQAN_TASSERT(str3 == "hello");
	SEQAN_TASSERT(length(str1) == 0);


}

//////////////////////////////////////////////////////////////////////////////

void Test_String_Array()
{
	TestStringBasics<String<char, Array<100> > >();

	String<char, Array<100> > str1 = "hello";
	SEQAN_TASSERT(str1[1] == 'e');

	String<char, Array<100> > const str2 = "hello";
	SEQAN_TASSERT(str2[4] == 'o');
	SEQAN_TASSERT(capacity(str2) >= length(str2));
}

//////////////////////////////////////////////////////////////////////////////

void Test_String_Stack()
{
	TestStringBasics<String<char, Block<3> > >();

	String<char, Block<3> > str1 = "hello";
	SEQAN_TASSERT(str1[1] == 'e');

	String<char, Block<3> > const str2 = "hello";
	SEQAN_TASSERT(str2[4] == 'o');
	SEQAN_TASSERT(capacity(str2) >= length(str2));

//	resize(str1, 30);
}

//////////////////////////////////////////////////////////////////////////////

void Test_String_Packed()
{
	TestStringBasics<String<char, Packed<> > >();
	TestStringBasics<String<Dna, Packed<> > >();

	TestStringResize<String<char, Packed<> > >();

}

//////////////////////////////////////////////////////////////////////////////

void Test_String_Pointer()
{
	char str1[200] = "hello";

	SEQAN_TASSERT(getValue(str1, 0) == 'h')
	SEQAN_TASSERT(getValue(str1, 1) == 'e')
	SEQAN_TASSERT(getValue("hello", 1) == 'e')

	int str2[100] = { 10, 20, 30, 0 };
	int const str3[100] = { 10, 20, 30, 0 };

	SEQAN_TASSERT(length(str2) == 3)
	SEQAN_TASSERT(length(str3) == 3)

	SEQAN_TASSERT(length(str1) == 5)
	SEQAN_TASSERT(length("hello") == 5)

	clear(str1);
	SEQAN_TASSERT(length(str1) == 0)	
	SEQAN_TASSERT(empty(str1))	

	SEQAN_TASSERT(reserve(str1, 100, Insist()) == 100);
	SEQAN_TASSERT(reserve(str1, 100, Limit()) == capacity(str1));

	fill(str1, 20, 'A');
	SEQAN_TASSERT(isEqual(str1, "AAAAAAAAAAAAAAAAAAAA"))

	resize(str1, 10);
	SEQAN_TASSERT(isEqual(str1, "AAAAAAAAAA"))

//____________________________________________________________________________
// compare operators 

//operators disabled due to ambiguity pointer/iterator vs. c-style string
/*
	assign(str1, "hello");
	String<char> str4 = str1;
	SEQAN_TASSERT(str1 == str4);

	str1 += str4;
	SEQAN_TASSERT(isEqual(str1, "hellohello"));

	SEQAN_TASSERT(str1 != str4);
	SEQAN_TASSERT(!isNotEqual(str1, "hellohello"));

	SEQAN_TASSERT(!(str1 < str4));
	SEQAN_TASSERT(!isLess(str1, str4));

	SEQAN_TASSERT(!(str1 <= str4));
	SEQAN_TASSERT(!isLessOrEqual(str1, str4));

	SEQAN_TASSERT(str1 > str4);
	SEQAN_TASSERT(isGreater(str1, str4));

	SEQAN_TASSERT(str1 >= str4);
	SEQAN_TASSERT(isGreaterOrEqual(str1, str4));
*/
}

//////////////////////////////////////////////////////////////////////////////

void Test_String_CStyle()
{
	String<char, CStyle> str1;
	char strq [200] = "hello seqan";
	
	String<char, CStyle> str2(strq);
	SEQAN_TASSERT(str2 == strq);
	SEQAN_TASSERT(strlen(str2) == strlen(strq));
	SEQAN_TASSERT(begin(str2) == begin(strq));
	SEQAN_TASSERT(end(str2) == end(strq));

	SEQAN_TASSERT(capacity(str2) == capacity(strq));

	String<char, CStyle> str3(str2);

	String<char const, CStyle> str4("a const string");

	String<char> stra("alloc string");
	String<char, CStyle>str5(stra);
	SEQAN_TASSERT(str5 == stra);

	String<char> const strac("const alloc string");
	String<char, CStyle> str7(strac);
	SEQAN_TASSERT(str7 == strac);

	str2 = stra;
	SEQAN_TASSERT(str2 == stra)

	str2 = strac;
	SEQAN_TASSERT(str2 == strac)

	str1 = str5;
	SEQAN_TASSERT(str1 == str5)

	char * cp1 = str1;
	SEQAN_TASSERT(cp1 != NULL);

	String<char, CStyle> const str6(strq);
	char const * cp2 = str6;
	SEQAN_TASSERT(cp2 != NULL);

	str1 = str6;
	SEQAN_TASSERT(str1 == str6)

	String<char, CStyle> str8(str6);
	SEQAN_TASSERT(str8 == str6)

	str2 = strac;
	SEQAN_TASSERT(str2 == strac)
	SEQAN_TASSERT(id(str2) != id(strac))

	clear(str2);
	SEQAN_TASSERT(length(str2) == 0);

	assign(str2, stra);
	SEQAN_TASSERT(str2 == stra)
	SEQAN_TASSERT(id(str2) == id(stra))
	SEQAN_TASSERT(str2 == toCString(stra));

	str2 = strac;
	SEQAN_TASSERT(str2 == strac)
	SEQAN_TASSERT(id(str2) != id(strac))

	String<Dna> str_dna("acgt");
	str2 = str_dna;
	SEQAN_TASSERT(str2 == str_dna)
	SEQAN_TASSERT(id(str2) != id(strac))

	String<Dna, CStyle> str9(str_dna);
	SEQAN_TASSERT(str2 == str_dna)

	char * strp = (char *) "this is a long array of chars";
	create(str2, strp);
	SEQAN_TASSERT(str2 == strp)
	SEQAN_TASSERT(id(str2) != id(strp))

	assign(str2, strp);
	SEQAN_TASSERT(str2 == strp)
	SEQAN_TASSERT(id(str2) == id(strp))

	str2 = "hello";
	String<char, CStyle > str10;
	move(str10, str2);
	SEQAN_TASSERT(str10 == "hello");
	SEQAN_TASSERT(length(str2) == 0);
}

//////////////////////////////////////////////////////////////////////////////

void Test_Segment()
{
//____________________________________________________________________________
// infix

	Infix<String<char> >::Type infix_1;
	String<char> str_1 = "this is a string";
	setHost(infix_1, str_1);
	SEQAN_TASSERT(length(infix_1) == 0)
	SEQAN_TASSERT(id(infix_1) == id(str_1))

	setEnd(infix_1, end(str_1));
	SEQAN_TASSERT(infix_1 == str_1)
	SEQAN_TASSERT(length(infix_1) == length(infix_1))

	setEndPosition(infix_1, 9);
	SEQAN_TASSERT(infix_1 == infix(str_1, 0, 9))

	Infix<String<char> >::Type infix_2(infix_1);
	SEQAN_TASSERT(infix_2 == infix(str_1, 0, 9))
	SEQAN_TASSERT(infix_2 == infix_1)
	SEQAN_TASSERT(id(infix_1) == id(infix_2))

	setBeginPosition(infix_2, 5);
	SEQAN_TASSERT(infix_2 == "is a")

	setBegin(infix_2, begin(str_1));
	SEQAN_TASSERT(infix_2 == "this is a")

	Infix<String<char> >::Type infix_3(str_1);
	SEQAN_TASSERT(infix_3 == getValue(str_1, 0))
	SEQAN_TASSERT(id(infix_3) == id(str_1))

	Infix<String<char> >::Type infix_4(str_1, 5, 9);
	SEQAN_TASSERT(infix_4 == "is a")

	SEQAN_TASSERT(capacity(infix_4) == capacity(str_1) - length(str_1) + length(infix_4))

	Infix<String<char> >::Type infix_5(str_1, begin(str_1), end(str_1));
	SEQAN_TASSERT(infix_5 == str_1)

	SEQAN_TASSERT(begin(infix_5) == begin(str_1))
	SEQAN_TASSERT(beginPosition(infix_5) == 0)
	SEQAN_TASSERT(end(infix_5) == end(str_1))
	SEQAN_TASSERT(endPosition(infix_5) == length(str_1))

	SEQAN_TASSERT(begin(infix(str_1, 0, length(str_1))) == begin(str_1))
	SEQAN_TASSERT(end(infix(str_1, 0, length(str_1))) == end(str_1))
	SEQAN_TASSERT(length(infix(str_1, 0, length(str_1))) == length(str_1))

	str_1 = "begin middle end";
	assign(infix(str_1, 6, 12),  "to");
	SEQAN_TASSERT(str_1 == "begin to end")

	assign(infix(str_1, 6, 8), "the test", 14);
	SEQAN_TASSERT(str_1 == "begin the test")

//	setEnd(infix_1);
//	SEQAN_TASSERT(infix_1 == "")

//____________________________________________________________________________
// test infix iteration

	str_1 = "begin middle end";
	goBegin(infix_1);
	SEQAN_TASSERT(infix_1 == "b")

	goBegin(infix_1, str_1);
	SEQAN_TASSERT(infix_1 == "b")

	goPrevious(infix_1);
	SEQAN_TASSERT(atBegin(infix_1))

	goEnd(infix_1);
	SEQAN_TASSERT(infix_1 == str_1)

	goEnd(infix_1, str_1);
	SEQAN_TASSERT(infix_1 == str_1)

	goNext(infix_1);
	SEQAN_TASSERT(atEnd(infix_1))

//____________________________________________________________________________
// suffix

	str_1 = "this is a string";

	Suffix<String<char> >::Type suffix_1;
	setHost(suffix_1, str_1);
	SEQAN_TASSERT(length(suffix_1) == length(str_1))
	SEQAN_TASSERT(id(suffix_1) == id(str_1))

	Suffix<String<char> >::Type suffix_2(suffix_1);
	SEQAN_TASSERT(suffix_2 == suffix_1)
	SEQAN_TASSERT(id(suffix_1) == id(suffix_2))

	setBeginPosition(suffix_2, 5);
	SEQAN_TASSERT(suffix_2 == "is a string")

	setBegin(suffix_2, begin(str_1));
	SEQAN_TASSERT(suffix_2 == "this is a string")

	Suffix<String<char> >::Type suffix_3(str_1);
	SEQAN_TASSERT(suffix_3 == str_1)
	SEQAN_TASSERT(id(suffix_3) == id(str_1))

	Suffix<String<char> >::Type suffix_4(str_1, 5);
	SEQAN_TASSERT(suffix_4 == "is a string")

	SEQAN_TASSERT(capacity(suffix_4) == capacity(str_1) - length(str_1) + length(suffix_4))

	Suffix<String<char> >::Type suffix_5(str_1, begin(str_1));
	SEQAN_TASSERT(suffix_5 == str_1)

	SEQAN_TASSERT(begin(suffix_5) == begin(str_1))
	SEQAN_TASSERT(beginPosition(suffix_5) == 0)
	SEQAN_TASSERT(end(suffix_5) == end(str_1))
	SEQAN_TASSERT(endPosition(suffix_5) == length(str_1))

	SEQAN_TASSERT(begin(suffix(str_1, 0)) == begin(str_1))
	SEQAN_TASSERT(end(suffix(str_1, 3)) == end(str_1))
	SEQAN_TASSERT(length(suffix(str_1, 0)) == length(str_1))

	str_1 = "begin middle end";
	assign(suffix(str_1, 6), "to panic");
	SEQAN_TASSERT(str_1 == "begin to panic")

	assign(suffix(str_1, 6), "the test", 9);
	SEQAN_TASSERT(str_1 == "begin the")

	char str_2[200] = "begin middle end";
	assign(suffix(str_2, 6), "to panic");
	SEQAN_TASSERT(isEqual(str_2, "begin to panic"))

	assign(suffix(str_2, 6), "the test", 9);
	SEQAN_TASSERT(isEqual(str_2, "begin the"))

//____________________________________________________________________________
// test suffix iteration
/*
	str_1 = "begin middle end";
	goBegin(suffix_1);
	SEQAN_TASSERT(suffix_1 == str_1)

	goBegin(suffix_1, str_1);
	SEQAN_TASSERT(suffix_1 == str_1)
	SEQAN_TASSERT(atBegin(suffix_1))

	goEnd(suffix_1);
	SEQAN_TASSERT(suffix_1 == "d")

	goEnd(suffix_1, str_1);
	SEQAN_TASSERT(suffix_1 == "d")

	goNext(suffix_1);
	SEQAN_TASSERT(atEnd(suffix_1))

	goPrevious(suffix_1);
	SEQAN_TASSERT(atEnd(suffix_1))
*/

//____________________________________________________________________________
// compare operators 

	str_1 = "hello";
	Infix<String<char> >::Type infix_6(str_1, 0, 5);

	SEQAN_TASSERT(infix_6 == str_1);

	infix_6 += str_1;
	SEQAN_TASSERT(isEqual(infix_6, "hellohello"));

	SEQAN_TASSERT(infix_6 != "bla");
	SEQAN_TASSERT(!isNotEqual(infix_6, "hellohello"));

	SEQAN_TASSERT(!(infix_6 < "hello"));
	SEQAN_TASSERT(!isLess(infix_6, "hello"));

	SEQAN_TASSERT(!(infix_6 <= "hello"));
	SEQAN_TASSERT(!isLessOrEqual(infix_6, "hello"));

	SEQAN_TASSERT(infix_6 > "hello");
	SEQAN_TASSERT(isGreater(infix_6, "hello"));

	SEQAN_TASSERT(infix_6 >= "hello");
	SEQAN_TASSERT(isGreaterOrEqual(infix_6, "hello"));
//____________________________________________________________________________

	clear(infix_6);
	SEQAN_TASSERT(infix_6 == "");
//____________________________________________________________________________
}

//////////////////////////////////////////////////////////////////////////////

void Test_Std_String()
{
//____________________________________________________________________________

	::std::string str_1("hamster");
	SEQAN_TASSERT((size_t)(end(str_1) - begin(str_1)) == length(str_1))

	::std::string const str_2("goldfish");
	SEQAN_TASSERT((size_t)(end(str_2) - begin(str_2)) == length(str_2))

	SEQAN_TASSERT(getValue(str_1, 0) == 'h')
	SEQAN_TASSERT(getValue(str_1, 1) == 'a')
	SEQAN_TASSERT(getValue(str_2, 1) == 'o')

	SEQAN_TASSERT(length(str_1) <= capacity(str_1))

	clear(str_1);
	SEQAN_TASSERT(empty(str_1))

	reserve(str_1, 200);
	SEQAN_TASSERT(capacity(str_1) >= 200)

	resize(str_1, 100);
	SEQAN_TASSERT(length(str_1) == 100)
	fill(str_1, 150, 'x');
	SEQAN_TASSERT(length(str_1) == 150)

//____________________________________________________________________________

}


//////////////////////////////////////////////////////////////////////////////

void Test_Lexical()
{
	Lexical<> lex1;
	compare(lex1, "abc", "abcd");

	SEQAN_TASSERT(isLess(lex1))
	SEQAN_TASSERT(isLess(lex1, TagPrefixLess()))
	SEQAN_TASSERT(!isLess(lex1, TagPrefixGreater()))

	SEQAN_TASSERT(isLessOrEqual(lex1))
	SEQAN_TASSERT(isLessOrEqual(lex1, TagPrefixLess()))
	SEQAN_TASSERT(!isLessOrEqual(lex1, TagPrefixGreater()))

	Lexical<> lex2(lex1);

	SEQAN_TASSERT(isPrefix(lex2))
	SEQAN_TASSERT(!hasPrefix(lex2))

	Lexical<> lex3("abcd", "abc");

	SEQAN_TASSERT(isGreater(lex3))
	SEQAN_TASSERT(isGreater(lex3, TagPrefixLess()))
	SEQAN_TASSERT(!isGreater(lex3, TagPrefixGreater()))

	SEQAN_TASSERT(isGreaterOrEqual(lex3))
	SEQAN_TASSERT(isGreaterOrEqual(lex3, TagPrefixLess()))
	SEQAN_TASSERT(!isGreaterOrEqual(lex3, TagPrefixGreater()))

	lex2 = lex3;
	SEQAN_TASSERT(!isPrefix(lex2))
	SEQAN_TASSERT(hasPrefix(lex2))

	String<char> str1 = "alpha";

	SEQAN_TASSERT(isEqual(str1, "alpha"))
	SEQAN_TASSERT(str1 == "alpha")

	SEQAN_TASSERT(isNotEqual(str1, "beta"))
	SEQAN_TASSERT(str1 != "beta")

	SEQAN_TASSERT(isLess(str1, "beta"))
	SEQAN_TASSERT(str1 < "beta")
	SEQAN_TASSERT(!isLess(str1, "aaa"))

	SEQAN_TASSERT(isLessOrEqual(str1, "beta"))
	SEQAN_TASSERT(str1 <= "beta")
	SEQAN_TASSERT(!isLessOrEqual(str1, "aaa"))
	SEQAN_TASSERT(isLessOrEqual(str1, "alpha"))

	SEQAN_TASSERT(isGreater(str1, "aaa"))
	SEQAN_TASSERT(str1 > "aaa")
	SEQAN_TASSERT(!isGreater(str1, "beta"))

	SEQAN_TASSERT(isGreaterOrEqual(str1, "aaa"))
	SEQAN_TASSERT(str1 >= "aaa")
	SEQAN_TASSERT(!isGreaterOrEqual(str1, "beta"))
	SEQAN_TASSERT(isGreaterOrEqual(str1, "alpha"))

	SEQAN_TASSERT(isPrefix(str1, "alpha romeo"))
	SEQAN_TASSERT(isPrefix(str1, "alpha"))
	SEQAN_TASSERT(!isPrefix(str1, ""))
	SEQAN_TASSERT(!isPrefix(str1, "alp"))
	SEQAN_TASSERT(!isPrefix(str1, "b"))

	SEQAN_TASSERT(hasPrefix(str1, "alp"))
	SEQAN_TASSERT(hasPrefix(str1, "alpha"))
	SEQAN_TASSERT(hasPrefix(str1, ""))
	SEQAN_TASSERT(!hasPrefix(str1, "alphas"))
	SEQAN_TASSERT(!hasPrefix(str1, "b"))

	SEQAN_TASSERT(lcpLength("hello", "hellmaker") == 4)
	SEQAN_TASSERT(lcpLength("", "not empty") == 0)
	SEQAN_TASSERT(lcpLength("hello", "hello") == 5)
	SEQAN_TASSERT(lcpLength("hello", "good evening") == 0)

	SEQAN_TASSERT(lcpLength("hello", 'h') == 1)
}

//////////////////////////////////////////////////////////////////////////////

template <typename TTarget, typename TSource, typename TExpand>
void Test_Assignments_Combinatoric(TTarget & target, TSource source, Tag<TExpand> const tag, size_t limit = ~0)
{
	assign(target, source, tag);
	SEQAN_TASSERT(infix(source, 0, length(target)) == target);

	assign(target, source, limit, tag);
	SEQAN_TASSERT(infix(source, 0, length(target)) == target);

	TSource const source_const(source);
	assign(target, source_const, tag);
	SEQAN_TASSERT(infix(source, 0, length(target)) == target);

	assign(target, source_const, limit, tag);
	SEQAN_TASSERT(infix(source, 0, length(target)) == target);

	typename Size<TTarget>::Type len = length(target);

	append(target, source, tag);
	if (len < length(target))
	{
		SEQAN_TASSERT(infix(source, 0, length(target) - len) == infix(target, len, length(target)) );
	}

	len = length(target);
	append(target, source, limit, tag);
	if (len < length(target))
	{
		SEQAN_TASSERT(infix(source, 0, length(target) - len) == infix(target, len, length(target)) );
	}

	append(target, source_const, tag);
	append(target, source_const, limit, tag); //p

	assign(target, "my miss is a hippi");
	replace(target, 9, 11, source, tag);
	SEQAN_TASSERT(infix(target, 0, 9) == "my miss i");

	len = length(target);
	if (len >= length(source) + 9)
	{
		SEQAN_TASSERT(infix(target, 9, 9 + length(source)) == source);
	}

	assign(target, "my miss is a hippi");
	replace(target, 9, 11, source, limit, tag);

	assign(target, "my miss is a hippi");
	replace(target, 9, 11, source, tag);
	SEQAN_TASSERT(infix(target, 0, 9) == "my miss i");

	len = length(target);
	if (len >= length(source) + 9)
	{
		SEQAN_TASSERT(infix(target, 9, 9 + length(source)) == source);
	}

	assign(target, "my miss is a hippi");
	replace(target, 9, 11, source, limit, tag);

	assign(target, "my miss is a hippi");
	replace(target, 9, 11, source_const, tag);

	assign(target, "my miss is a hippi");
	replace(target, 9, 11, source_const, limit, tag);
}

void Test_Combinatoric()
{
	String<char> str1("hello");
	String<char> str2("this is test");
	String<char> const str3("this is const string");

	Test_Assignments_Combinatoric(str1, str2, Exact());
	Test_Assignments_Combinatoric(str1, str2, Generous());
	Test_Assignments_Combinatoric(str1, str2, Limit());

	Test_Assignments_Combinatoric(str1, str2, Exact(), 3);
	Test_Assignments_Combinatoric(str1, str2, Generous(), 3);
	Test_Assignments_Combinatoric(str1, str2, Limit(), 3);

	Test_Assignments_Combinatoric(str1, str1, Exact());
	Test_Assignments_Combinatoric(str1, str1, Generous());
	Test_Assignments_Combinatoric(str1, str1, Limit());

	Test_Assignments_Combinatoric(str1, str1, Exact(), 3);
	Test_Assignments_Combinatoric(str1, str1, Generous(), 3);
	Test_Assignments_Combinatoric(str1, str1, Limit(), 3);

	Test_Assignments_Combinatoric(str1, str3, Exact());
	Test_Assignments_Combinatoric(str1, str3, Generous());
	Test_Assignments_Combinatoric(str1, str3, Limit());

	Test_Assignments_Combinatoric(str1, str3, Exact(), 3);
	Test_Assignments_Combinatoric(str1, str3, Generous(), 3);
	Test_Assignments_Combinatoric(str1, str3, Limit(), 3);

//____________________________________________________________________________

	String<char> str4("hello");
	Segment<String<char> > str5(str4, 1, 2);

	Test_Assignments_Combinatoric(str1, str5, Exact());
	Test_Assignments_Combinatoric(str1, str5, Generous());
	Test_Assignments_Combinatoric(str1, str5, Limit());

	Test_Assignments_Combinatoric(str1, str5, Exact(), 3);
	Test_Assignments_Combinatoric(str1, str5, Generous(), 3);
	Test_Assignments_Combinatoric(str1, str5, Limit(), 3);

	char str6 = 'x';
	Test_Assignments_Combinatoric(str1, str6, Exact());
	Test_Assignments_Combinatoric(str1, str6, Generous());
	Test_Assignments_Combinatoric(str1, str6, Limit());

	char str7[800] = "hello again";
	reserve(str1, 10000);
	Test_Assignments_Combinatoric(str1, str7, Insist());
	Test_Assignments_Combinatoric(str7, str4, Insist());
	Test_Assignments_Combinatoric(str7, str7, Insist());
	Test_Assignments_Combinatoric(str7, "sisyphos", Insist());

	Test_Assignments_Combinatoric(str7, 'c', Insist());

	Test_Assignments_Combinatoric(str1, str7, Insist(), 3);
	Test_Assignments_Combinatoric(str7, str4, Insist(), 3);
	Test_Assignments_Combinatoric(str7, str7, Insist(), 3);
//____________________________________________________________________________

	assign(str7, "begin middle end");
	Segment<char *> infix_1(str7, 6, 12);
	SEQAN_TASSERT(infix_1 == "middle")

	Test_Assignments_Combinatoric(infix_1, str1, Insist());
	SEQAN_TASSERT(lcpLength(end(infix_1, Standard()), " end") >= 4)
	SEQAN_TASSERT(beginPosition(infix_1) == 6)
	SEQAN_TASSERT(infix(str7, 0, 6) == "begin ")

	Test_Assignments_Combinatoric(infix_1, str1, Insist(), 10);
	str4 = "begin middle end";
	Infix<String<char> >::Type infix_2(str4, 6, 12);
	SEQAN_TASSERT(infix_2 == "middle");

	Test_Assignments_Combinatoric(infix_2, str1, Exact());
	Test_Assignments_Combinatoric(infix_2, str1, Generous());
	Test_Assignments_Combinatoric(infix_2, str1, Limit());

	SEQAN_TASSERT(beginPosition(infix_2) == 6)
	SEQAN_TASSERT(infix(str4, 0, 6) == "begin ")

	Test_Assignments_Combinatoric(infix_2, str1, Exact(), 10);
	Test_Assignments_Combinatoric(infix_2, str1, Generous(), 10);
	Test_Assignments_Combinatoric(infix_2, str1, Limit(), 10);

//____________________________________________________________________________

	str1 = "seqan string";
	std::string str8("i am the standard");
	Test_Assignments_Combinatoric(str8, str1, Generous());
	Test_Assignments_Combinatoric(str8, str1, Limit());

	Test_Assignments_Combinatoric(str8, str1, Generous(), 10);
	Test_Assignments_Combinatoric(str8, str1, Limit(), 10);

	str8 = "standard string";
	Test_Assignments_Combinatoric(str1, str8, Generous());
	Test_Assignments_Combinatoric(str1, str8, Generous(), 10);

	str8 = "standard string";
	Test_Assignments_Combinatoric(str7, str8, Insist());
	Test_Assignments_Combinatoric(str7, str8, Insist(), 10);

//____________________________________________________________________________

	str1 = "this is a test string";
	String<char, Array<100> > str9;
	Test_Assignments_Combinatoric(str9, str1, Limit());
	Test_Assignments_Combinatoric(str9, str1, Limit(), 10);
}

//////////////////////////////////////////////////////////////////////////////


int mainTestString()  
{
	SEQAN_TREPORT("TEST STRING BEGIN")

	Test_Sequence_Interface();
	Test_String_Base();
	Test_String_Alloc();
	Test_String_Array();
	Test_String_Stack();
	Test_String_Pointer();
	Test_String_CStyle();
	Test_String_Packed();
	Test_Std_String();

	Test_Lexical();
	Test_Combinatoric();

	Test_Segment();

	debug::verifyCheckpoints("projects/library/seqan/sequence.h");
	debug::verifyCheckpoints("projects/library/seqan/sequence/sequence_interface.h");
	debug::verifyCheckpoints("projects/library/seqan/sequence/string_base.h");
	debug::verifyCheckpoints("projects/library/seqan/sequence/string_alloc.h");
	debug::verifyCheckpoints("projects/library/seqan/sequence/string_pointer.h");
	debug::verifyCheckpoints("projects/library/seqan/sequence/string_array.h");
	debug::verifyCheckpoints("projects/library/seqan/sequence/string_cstyle.h");
	debug::verifyCheckpoints("projects/library/seqan/sequence/lexical.h");
	debug::verifyCheckpoints("projects/library/seqan/sequence/std_string.h");

	SEQAN_TREPORT("TEST STRING END")

	return 0;
}

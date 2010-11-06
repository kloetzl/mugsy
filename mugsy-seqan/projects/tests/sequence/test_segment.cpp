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

void Test_Infix()
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

void Test_Suffix()
{
//____________________________________________________________________________
// suffix

	String<char> str_1 = "this is a string";

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

}

//////////////////////////////////////////////////////////////////////////////

int mainTestSegment()  
{
	SEQAN_TREPORT("TEST SEGMENT BEGIN")

	Test_Infix();
	Test_Suffix();

	debug::verifyCheckpoints("projects/library/seqan/sequence/segment_infix.h");
	debug::verifyCheckpoints("projects/library/seqan/sequence/segment_suffix.h");
	debug::verifyCheckpoints("projects/library/seqan/sequence/segment_prefix.h");
	//debug::verifyCheckpoints("projects/library/seqan/sequence/segment_gram.h");

	SEQAN_TREPORT("TEST SEGMENT END")

	return 0;
}

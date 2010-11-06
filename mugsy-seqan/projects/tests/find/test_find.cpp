#include <iostream>
#include <fstream>
#include <typeinfo>
#include <time.h>
#include <cstdio>
#include <vector>
#include <time.h>

#define SEQAN_DEBUG
#define SEQAN_TEST

#define TEST_PATH "projects/tests/find/"
#define LIB_PATH "projects/library/seqan/find/"

#include "seqan/find.h"

using namespace std;
using namespace seqan;

//////////////////////////////////////////////////////////////////////////////

template <typename TAlgorithmSpec>
void Test_OnlineAlg()
{
	String<unsigned int> pos;

//____________________________________________________________________________
// Test1 - small needle

	String<char> haystack("Dies ist ein Haystack. Ja, das ist wirklich einer!");
	Finder<String<char> > finder(haystack);

	String<char> needle("ist");
	Pattern<String<char>, TAlgorithmSpec> pattern(needle);

	while (find(finder, pattern))
	{
		append(pos,position(finder));
		SEQAN_TASSERT(position(finder) == beginPosition(finder))
		SEQAN_TASSERT(endPosition(finder) == beginPosition(finder) + length(finder))
		SEQAN_TASSERT(length(finder) == length(needle))
		SEQAN_TASSERT(begin(finder) == begin(haystack) + beginPosition(finder))
		SEQAN_TASSERT(end(finder) == begin(haystack) + endPosition(finder))
		SEQAN_TASSERT(infix(finder) == needle)
	}

	SEQAN_TASSERT(host(pattern) == needle);
	SEQAN_TASSERT(host(reinterpret_cast<Pattern<String<char>, TAlgorithmSpec> const &>(pattern)) == needle);
	SEQAN_TASSERT(pos[0] == 5);
	SEQAN_TASSERT(pos[1] == 31);
	SEQAN_TASSERT(length(pos) == 2);

//____________________________________________________________________________
// Test2 - large needle

	haystack = "abcdefghijklmnopqrstuvwxyzabcdefghijklmnopqrstuvwxyzabcdefgaabcdef";
	setHost(finder, haystack);
	clear(finder);

	needle = "abcdefghijklmnopqrstuvwxyzabcdefg";
	setHost(pattern, needle);

	clear(pos);
	while (find(finder, pattern))
	{
		append(pos,position(finder));
		SEQAN_TASSERT(position(finder) == beginPosition(finder))
		SEQAN_TASSERT(endPosition(finder) == beginPosition(finder) + length(finder))
		SEQAN_TASSERT(length(finder) == length(needle))
		SEQAN_TASSERT(begin(finder) == begin(haystack) + beginPosition(finder))
		SEQAN_TASSERT(end(finder) == begin(haystack) + endPosition(finder))
		SEQAN_TASSERT(infix(finder) == needle)
	}

	SEQAN_TASSERT(pos[0] == 0);
	SEQAN_TASSERT(pos[1] == 26);
	SEQAN_TASSERT(length(pos) == 2);

//____________________________________________________________________________
// Test3 - different alphabet, small needle

	String<Dna> hstk = "aaaaaaacaa";
	Finder<String<Dna> > finderDna(hstk);

	String<Dna> ndl = "aa";
	setHost(pattern, ndl);

	clear(pos);
	while (find(finderDna, pattern))
	{
		append(pos,position(finderDna));
		SEQAN_TASSERT(position(finderDna) == beginPosition(finderDna))
		SEQAN_TASSERT(endPosition(finderDna) == beginPosition(finderDna) + length(finderDna))
		SEQAN_TASSERT(length(finderDna) == length(ndl))
		SEQAN_TASSERT(begin(finderDna) == begin(hstk) + beginPosition(finderDna))
		SEQAN_TASSERT(end(finderDna) == begin(hstk) + endPosition(finderDna))
		SEQAN_TASSERT(infix(finderDna) == ndl)
	}

	SEQAN_TASSERT(pos[0] == 0);
	SEQAN_TASSERT(pos[1] == 1);
	SEQAN_TASSERT(pos[2] == 2);
	SEQAN_TASSERT(pos[3] == 3);
	SEQAN_TASSERT(pos[4] == 4);
	SEQAN_TASSERT(pos[5] == 5);
	SEQAN_TASSERT(pos[6] == 8);
	SEQAN_TASSERT(length(pos) == 7);

//____________________________________________________________________________
// Test3b - different alphabet, small needle, jumping finder

	goBegin(finderDna); // That's a repositioning
	clear(finderDna);	// That's why, clear state
	clear(pos);

	bool firstHit = true;
	while (find(finderDna, pattern)) {
		if (firstHit) {
			firstHit = false;
			finderDna += 2;
			clear(finderDna);  // clear the state of the finder
		} else {
			//unsigned int p = position(finderDna);
			append(pos,position(finderDna));
		}
	}

	SEQAN_TASSERT(pos[0] == 2);
	SEQAN_TASSERT(pos[1] == 3);
	SEQAN_TASSERT(pos[2] == 4);
	SEQAN_TASSERT(pos[3] == 5);
	SEQAN_TASSERT(pos[4] == 8);
	SEQAN_TASSERT(length(pos) == 5);

//____________________________________________________________________________
// Test4 - different alphabet, large needle
	String<Dna> text = "taaaataaaataaaataaaataaaataaaataaaataaaataaaataaaataaaataaaataaaataaaat";
	Finder<String<Dna> > finderText(text);

	String<Dna> query = "taaaataaaataaaataaaataaaataaaataaaataaaataaaat";
	setHost(pattern, query);

	clear(pos);
	while (find(finderText, pattern))
	{
		append(pos,position(finderText));
		SEQAN_TASSERT(position(finderText) == beginPosition(finderText))
		SEQAN_TASSERT(endPosition(finderText) == beginPosition(finderText) + length(finderText))
		SEQAN_TASSERT(length(finderText) == length(query))
		SEQAN_TASSERT(begin(finderText) == begin(text) + beginPosition(finderText))
		SEQAN_TASSERT(end(finderText) == begin(text) + endPosition(finderText))
		SEQAN_TASSERT(infix(finderText) == query)
	}

	SEQAN_TASSERT(pos[0] == 0);
	SEQAN_TASSERT(pos[1] == 5);
	SEQAN_TASSERT(pos[2] == 10);
	SEQAN_TASSERT(pos[3] == 15);
	SEQAN_TASSERT(pos[4] == 20);
	SEQAN_TASSERT(pos[5] == 25);
	SEQAN_TASSERT(length(pos) == 6);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TAlgorithmSpec>
void Test_OnlineAlgMulti(bool order_by_begin_position)
{
	String<unsigned int> pos;

//____________________________________________________________________________
// Test1 - Single keyword
	String<char> haystack("Dies ist ein Haystack. Ja, das ist wirklich einer!");
	Finder<String<char> > finder(haystack);

	typedef String<String<char> > TNeedle;
	TNeedle keywords;
	appendValue(keywords, String<char>("ist"));
	Pattern<TNeedle, TAlgorithmSpec> pattern(keywords);

	while (find(finder, pattern))
	{
		append(pos,position(finder));
		SEQAN_TASSERT(position(finder) == beginPosition(finder))
		SEQAN_TASSERT(endPosition(finder) == beginPosition(finder) + length(finder))
		SEQAN_TASSERT(length(finder) == length(keywords[position(pattern)]))
		SEQAN_TASSERT(begin(finder) == begin(haystack) + beginPosition(finder))
		SEQAN_TASSERT(end(finder) == begin(haystack) + endPosition(finder))
		SEQAN_TASSERT(infix(finder) == keywords[position(pattern)])
	}

	SEQAN_TASSERT(host(pattern) == keywords);
	SEQAN_TASSERT(host(reinterpret_cast<Pattern<TNeedle, TAlgorithmSpec> const &>(pattern)) == keywords);
	SEQAN_TASSERT(pos[0] == 5);
	SEQAN_TASSERT(pos[1] == 31);
	SEQAN_TASSERT(length(pos) == 2);

	haystack = "abcdefghijklmnopqrstuvwxyzabcdefghijklmnopqrstuvwxyzabcdefgaabcdef";
	setHost(finder, haystack);
	clear(finder);

	clear(keywords);
	appendValue(keywords, String<char>("abcdefghijklmnopqrstuvwxyzabcdefg"));
	setHost(pattern, keywords);
	clear(pos);

	while (find(finder, pattern))
	{
		append(pos,position(finder));
		SEQAN_TASSERT(position(finder) == beginPosition(finder))
		SEQAN_TASSERT(endPosition(finder) == beginPosition(finder) + length(finder))
		SEQAN_TASSERT(length(finder) == length(keywords[position(pattern)]))
		SEQAN_TASSERT(begin(finder) == begin(haystack) + beginPosition(finder))
		SEQAN_TASSERT(end(finder) == begin(haystack) + endPosition(finder))
		SEQAN_TASSERT(infix(finder) == keywords[position(pattern)])
	}

	SEQAN_TASSERT(pos[0] == 0);
	SEQAN_TASSERT(pos[1] == 26);
	SEQAN_TASSERT(length(pos) == 2);


	String<Dna> hstk = "aaaaaaacaa";
	Finder<String<Dna> > finderDna(hstk);

	typedef String<String<Dna> > TDnaNeedle;
	Pattern<TDnaNeedle, TAlgorithmSpec> pattern_dna(keywords);

	TDnaNeedle dna_keywords;
	appendValue(dna_keywords, String<Dna>("aa"));
	setHost(pattern_dna, dna_keywords);

	clear(pos);
	while (find(finderDna, pattern_dna))
	{
		append(pos,position(finderDna));
		SEQAN_TASSERT(position(finderDna) == beginPosition(finderDna))
		SEQAN_TASSERT(endPosition(finderDna) == beginPosition(finderDna) + length(finderDna))
		SEQAN_TASSERT(length(finderDna) == length(dna_keywords[position(pattern_dna)]))
		SEQAN_TASSERT(begin(finderDna) == begin(hstk) + beginPosition(finderDna))
		SEQAN_TASSERT(end(finderDna) == begin(hstk) + endPosition(finderDna))
		SEQAN_TASSERT(infix(finderDna) == dna_keywords[position(pattern_dna)])
	}

	SEQAN_TASSERT(pos[0] == 0);
	SEQAN_TASSERT(pos[1] == 1);
	SEQAN_TASSERT(pos[2] == 2);
	SEQAN_TASSERT(pos[3] == 3);
	SEQAN_TASSERT(pos[4] == 4);
	SEQAN_TASSERT(pos[5] == 5);
	SEQAN_TASSERT(pos[6] == 8);
	SEQAN_TASSERT(length(pos) == 7);

	String<Dna> text = "taaaataaaataaaataaaataaaataaaataaaataaaataaaataaaataaaataaaataaaataaaat";
	Finder<String<Dna> > finderText(text);

	clear(dna_keywords);
	appendValue(dna_keywords, String<Dna>("taaaataaaataaaataaaataaaataaaataaaataaaataaaat"));
	setHost(pattern_dna, dna_keywords);

	clear(pos);
	while (find(finderText, pattern_dna)) 
	{
		append(pos,position(finderText));
		SEQAN_TASSERT(position(finderText) == beginPosition(finderText))
		SEQAN_TASSERT(endPosition(finderText) == beginPosition(finderText) + length(finderText))
		SEQAN_TASSERT(length(finderText) == length(dna_keywords[position(pattern_dna)]))
		SEQAN_TASSERT(begin(finderText) == begin(text) + beginPosition(finderText))
		SEQAN_TASSERT(end(finderText) == begin(text) + endPosition(finderText))
		SEQAN_TASSERT(infix(finderText) == dna_keywords[position(pattern_dna)])
	}

	SEQAN_TASSERT(pos[0] == 0);
	SEQAN_TASSERT(pos[1] == 5);
	SEQAN_TASSERT(pos[2] == 10);
	SEQAN_TASSERT(pos[3] == 15);
	SEQAN_TASSERT(pos[4] == 20);
	SEQAN_TASSERT(pos[5] == 25);
	SEQAN_TASSERT(length(pos) == 6);

//____________________________________________________________________________
// Test2 - Multiple keywords
	String<char> hst("annual_announce_any_annually");
	Finder<String<char> > fd(hst);

	typedef String<String<char> > TN;
	TN kyw;
	appendValue(kyw, String<char>("announce"));
	appendValue(kyw, String<char>("annual"));
	appendValue(kyw, String<char>("annually"));
	Pattern<TN, TAlgorithmSpec> pt(kyw);

	String<unsigned int> finderPos;
	String<unsigned int> keywordIndex;
	while (find(fd, pt)) 
	{
		append(finderPos,position(fd));
		append(keywordIndex,position(pt));
		SEQAN_TASSERT(position(fd) == beginPosition(fd))
		SEQAN_TASSERT(endPosition(fd) == beginPosition(fd) + length(fd))
		SEQAN_TASSERT(length(fd) == length(kyw[position(pt)]))
		SEQAN_TASSERT(begin(fd) == begin(hst) + beginPosition(fd))
		SEQAN_TASSERT(end(fd) == begin(hst) + endPosition(fd))
		SEQAN_TASSERT(infix(fd) == kyw[position(pt)])
	}

	SEQAN_TASSERT(length(finderPos) == 4);
	SEQAN_TASSERT(length(keywordIndex) == 4);
	SEQAN_TASSERT(finderPos[0] == 0);
	SEQAN_TASSERT(keywordIndex[0] == 1);
	SEQAN_TASSERT(finderPos[1] == 7);
	SEQAN_TASSERT(keywordIndex[1] == 0);
	SEQAN_TASSERT(finderPos[2] == 20);
	SEQAN_TASSERT(keywordIndex[2] == 1);
	SEQAN_TASSERT(finderPos[3] == 20);
	SEQAN_TASSERT(keywordIndex[3] == 2);

	String<Dna> hstDna("AGATACGATATATAC");
	Finder<String<Dna> > fdDna(hstDna);

	typedef String<String<Dna> > TNDna;
	TNDna kywDna;
	appendValue(kywDna, String<Dna>("ATATATA"));
	appendValue(kywDna, String<Dna>("TATAT"));
	appendValue(kywDna, String<Dna>("ACGATAT"));
	Pattern<TNDna, TAlgorithmSpec> ptDna(kywDna);

	clear(finderPos);
	clear(keywordIndex);
	while (find(fdDna, ptDna)) 
	{
		append(finderPos,position(fdDna));
		append(keywordIndex,position(ptDna));
		SEQAN_TASSERT(position(fdDna) == beginPosition(fdDna))
		SEQAN_TASSERT(endPosition(fdDna) == beginPosition(fdDna) + length(fdDna))
		SEQAN_TASSERT(length(fdDna) == length(kywDna[position(ptDna)]))
		SEQAN_TASSERT(begin(fdDna) == begin(hstDna) + beginPosition(fdDna))
		SEQAN_TASSERT(end(fdDna) == begin(hstDna) + endPosition(fdDna))
		SEQAN_TASSERT(infix(fdDna) == kywDna[position(ptDna)])
	}

	SEQAN_TASSERT(length(finderPos) == 3);
	SEQAN_TASSERT(length(keywordIndex) == 3);
	if (order_by_begin_position)
	{
		SEQAN_TASSERT(finderPos[0] == 4);
		SEQAN_TASSERT(keywordIndex[0] == 2);
		SEQAN_TASSERT(finderPos[1] == 7);
		SEQAN_TASSERT(keywordIndex[1] == 0);
		SEQAN_TASSERT(finderPos[2] == 8);
		SEQAN_TASSERT(keywordIndex[2] == 1);
	}
	else
	{
		SEQAN_TASSERT(finderPos[0] == 4);
		SEQAN_TASSERT(keywordIndex[0] == 2);
		SEQAN_TASSERT(finderPos[1] == 8);
		SEQAN_TASSERT(keywordIndex[1] == 1);
		SEQAN_TASSERT(finderPos[2] == 7);
		SEQAN_TASSERT(keywordIndex[2] == 0);
	}
//____________________________________________________________________________
// Test2 - Multiple keywords that do not fit into a machine word
	String<Dna> my_haystack("AGATACGATATATACAGATACGATATATACAGATACGATATATACAGATACGATATATACAGATACGATATATAC");
	Finder<String<Dna> > my_finder(my_haystack);

	typedef String<String<Dna> > TNeedle_My;
	TNeedle_My my_keywords;
	appendValue(my_keywords, String<Dna>("ATATATA"));
	appendValue(my_keywords, String<Dna>("ACCGATCCAT"));
	appendValue(my_keywords, String<Dna>("TATAT"));
	appendValue(my_keywords, String<Dna>("ACCGAT"));
	appendValue(my_keywords, String<Dna>("ACGATAT"));
	appendValue(my_keywords, String<Dna>("CCAA"));
	Pattern<TNeedle_My, TAlgorithmSpec> my_pattern(my_keywords);

	clear(finderPos);
	clear(keywordIndex);
	while (find(my_finder, my_pattern)) {
		//std::cout << position(my_finder) << "-" << position(my_pattern) << ::std::endl;
		append(finderPos,position(my_finder));
		append(keywordIndex,position(my_pattern));
		SEQAN_TASSERT(position(my_finder) == beginPosition(my_finder))
		SEQAN_TASSERT(endPosition(my_finder) == beginPosition(my_finder) + length(my_finder))
		SEQAN_TASSERT(length(my_finder) == length(my_keywords[position(my_pattern)]))
		SEQAN_TASSERT(begin(my_finder) == begin(my_haystack) + beginPosition(my_finder))
		SEQAN_TASSERT(end(my_finder) == begin(my_haystack) + endPosition(my_finder))
		SEQAN_TASSERT(infix(my_finder) == my_keywords[position(my_pattern)])
	}

	SEQAN_TASSERT(length(finderPos) == 15);
	SEQAN_TASSERT(length(keywordIndex) == 15);
	if (order_by_begin_position)
	{
		SEQAN_TASSERT(finderPos[0] == 4);
		SEQAN_TASSERT(keywordIndex[0] == 4);
		SEQAN_TASSERT(finderPos[1] == 7);
		SEQAN_TASSERT(keywordIndex[1] == 0);
		SEQAN_TASSERT(finderPos[2] == 8);
		SEQAN_TASSERT(keywordIndex[2] == 2);
		SEQAN_TASSERT(finderPos[3] == 19);
		SEQAN_TASSERT(keywordIndex[3] == 4);
		SEQAN_TASSERT(finderPos[4] == 22);
		SEQAN_TASSERT(keywordIndex[4] == 0);
		SEQAN_TASSERT(finderPos[5] == 23);
		SEQAN_TASSERT(keywordIndex[5] == 2);
		SEQAN_TASSERT(finderPos[6] == 34);
		SEQAN_TASSERT(keywordIndex[6] == 4);
		SEQAN_TASSERT(finderPos[7] == 37);
		SEQAN_TASSERT(keywordIndex[7] == 0);
		SEQAN_TASSERT(finderPos[8] == 38);
		SEQAN_TASSERT(keywordIndex[8] == 2);
		SEQAN_TASSERT(finderPos[9] == 49);
		SEQAN_TASSERT(keywordIndex[9] == 4);
		SEQAN_TASSERT(finderPos[10] == 52);
		SEQAN_TASSERT(keywordIndex[10] == 0);
		SEQAN_TASSERT(finderPos[11] == 53);
		SEQAN_TASSERT(keywordIndex[11] == 2);
		SEQAN_TASSERT(finderPos[12] == 64);
		SEQAN_TASSERT(keywordIndex[12] == 4);
		SEQAN_TASSERT(finderPos[13] == 67);
		SEQAN_TASSERT(keywordIndex[13] == 0);
		SEQAN_TASSERT(finderPos[14] == 68);
		SEQAN_TASSERT(keywordIndex[14] == 2);
	}
	else
	{
		SEQAN_TASSERT(finderPos[0] == 4);
		SEQAN_TASSERT(keywordIndex[0] == 4);
		SEQAN_TASSERT(finderPos[1] == 8);
		SEQAN_TASSERT(keywordIndex[1] == 2);
		SEQAN_TASSERT(finderPos[2] == 7);
		SEQAN_TASSERT(keywordIndex[2] == 0);
		SEQAN_TASSERT(finderPos[3] == 19);
		SEQAN_TASSERT(keywordIndex[3] == 4);
		SEQAN_TASSERT(finderPos[4] == 23);
		SEQAN_TASSERT(keywordIndex[4] == 2);
		SEQAN_TASSERT(finderPos[5] == 22);
		SEQAN_TASSERT(keywordIndex[5] == 0);
		SEQAN_TASSERT(finderPos[6] == 34);
		SEQAN_TASSERT(keywordIndex[6] == 4);
		SEQAN_TASSERT(finderPos[7] == 38);
		SEQAN_TASSERT(keywordIndex[7] == 2);
		SEQAN_TASSERT(finderPos[8] == 37);
		SEQAN_TASSERT(keywordIndex[8] == 0);
		SEQAN_TASSERT(finderPos[9] == 49);
		SEQAN_TASSERT(keywordIndex[9] == 4);
		SEQAN_TASSERT(finderPos[10] == 53);
		SEQAN_TASSERT(keywordIndex[10] == 2);
		SEQAN_TASSERT(finderPos[11] == 52);
		SEQAN_TASSERT(keywordIndex[11] == 0);
		SEQAN_TASSERT(finderPos[12] == 64);
		SEQAN_TASSERT(keywordIndex[12] == 4);
		SEQAN_TASSERT(finderPos[13] == 68);
		SEQAN_TASSERT(keywordIndex[13] == 2);
		SEQAN_TASSERT(finderPos[14] == 67);
		SEQAN_TASSERT(keywordIndex[14] == 0);
	}

//____________________________________________________________________________
// Multiple keywords with overlapping matches
	String<Dna> my2_haystack("aaaacaaa");
	Finder<String<Dna> > my2_finder(my2_haystack);

	typedef String<String<Dna> > TNeedle_My2;
	TNeedle_My2 my2_keywords;
	appendValue(my2_keywords, String<Dna>("aa"));
	appendValue(my2_keywords, String<Dna>("aaa"));
	appendValue(my2_keywords, String<Dna>("ac"));
	appendValue(my2_keywords, String<Dna>("aac"));
	appendValue(my2_keywords, String<Dna>("gctccacctgacctagcccatggggcccaaatttccggccttaattcccattt"));
	Pattern<TNeedle_My2, TAlgorithmSpec> my2_pattern(my2_keywords);

	clear(finderPos);
	clear(keywordIndex);
	while (find(my2_finder, my2_pattern)) {
		//std::cout << position(my2_finder) << ":" << position(my2_pattern) << ::std::endl;
		append(finderPos,position(my2_finder));
		append(keywordIndex,position(my2_pattern));
		SEQAN_TASSERT(position(my2_finder) == beginPosition(my2_finder))
		SEQAN_TASSERT(endPosition(my2_finder) == beginPosition(my2_finder) + length(my2_finder))
		SEQAN_TASSERT(length(my2_finder) == length(my2_keywords[position(my2_pattern)]))
		SEQAN_TASSERT(begin(my2_finder) == begin(my2_haystack) + beginPosition(my2_finder))
		SEQAN_TASSERT(end(my2_finder) == begin(my2_haystack) + endPosition(my2_finder))
		SEQAN_TASSERT(infix(my2_finder) == my2_keywords[position(my2_pattern)])
	}

	if (order_by_begin_position)
	{
		SEQAN_TASSERT(finderPos[0] == 0);
		SEQAN_TASSERT(keywordIndex[0] == 0);
		SEQAN_TASSERT(finderPos[1] == 0);
		SEQAN_TASSERT(keywordIndex[1] == 1);
		SEQAN_TASSERT(finderPos[2] == 1);
		SEQAN_TASSERT(keywordIndex[2] == 0);
		SEQAN_TASSERT(finderPos[3] == 1);
		SEQAN_TASSERT(keywordIndex[3] == 1);
		SEQAN_TASSERT(finderPos[4] == 2);
		SEQAN_TASSERT(keywordIndex[4] == 0);
		SEQAN_TASSERT(finderPos[5] == 2);
		SEQAN_TASSERT(keywordIndex[5] == 3);
		SEQAN_TASSERT(finderPos[6] == 3);
		SEQAN_TASSERT(keywordIndex[6] == 2);
		SEQAN_TASSERT(finderPos[7] == 5);
		SEQAN_TASSERT(keywordIndex[7] == 0);
		SEQAN_TASSERT(finderPos[8] == 5);
		SEQAN_TASSERT(keywordIndex[8] == 1);
		SEQAN_TASSERT(finderPos[9] == 6);
		SEQAN_TASSERT(keywordIndex[9] == 0);
	}
	else
	{
		SEQAN_TASSERT(finderPos[0] == 0);
		SEQAN_TASSERT(keywordIndex[0] == 0);
		SEQAN_TASSERT(finderPos[1] == 1);
		SEQAN_TASSERT(keywordIndex[1] == 0);
		SEQAN_TASSERT(finderPos[2] == 0);
		SEQAN_TASSERT(keywordIndex[2] == 1);
		SEQAN_TASSERT(finderPos[3] == 2);
		SEQAN_TASSERT(keywordIndex[3] == 0);
		SEQAN_TASSERT(finderPos[4] == 1);
		SEQAN_TASSERT(keywordIndex[4] == 1);
		SEQAN_TASSERT(finderPos[5] == 3);
		SEQAN_TASSERT(keywordIndex[5] == 2);
		SEQAN_TASSERT(finderPos[6] == 2);
		SEQAN_TASSERT(keywordIndex[6] == 3);
		SEQAN_TASSERT(finderPos[7] == 5);
		SEQAN_TASSERT(keywordIndex[7] == 0);
		SEQAN_TASSERT(finderPos[8] == 6);
		SEQAN_TASSERT(keywordIndex[8] == 0);
		SEQAN_TASSERT(finderPos[9] == 5);
		SEQAN_TASSERT(keywordIndex[9] == 1);
	}

//____________________________________________________________________________
// Multiple duplicated keywords with overlapping matches, jumping finder
	goBegin(my2_finder); // That's a repositioning
	clear(my2_finder);	// That's why, clear state
	clear(finderPos);
	clear(keywordIndex);

	unsigned int hits = 0;
	while (find(my2_finder, my2_pattern)) {
		if (hits == 2) break;
		//std::cout << position(my2_finder) << ":" << position(my2_pattern) << ::std::endl;
		append(finderPos,position(my2_finder));
		append(keywordIndex,position(my2_pattern));
		++hits;
	}
	goBegin(my2_finder);
	my2_finder+=5;
	clear(my2_finder);
	clear(finderPos);
	clear(keywordIndex);
	while (find(my2_finder, my2_pattern)) {
		//std::cout << position(my2_finder) << ":" << position(my2_pattern) << ::std::endl;
		append(finderPos,position(my2_finder));
		append(keywordIndex,position(my2_pattern));
	}

	if (order_by_begin_position)
	{
		SEQAN_TASSERT(finderPos[0] == 5);
		SEQAN_TASSERT(keywordIndex[0] == 0);
		SEQAN_TASSERT(finderPos[1] == 5);
		SEQAN_TASSERT(keywordIndex[1] == 1);
		SEQAN_TASSERT(finderPos[2] == 6);
		SEQAN_TASSERT(keywordIndex[2] == 0);
	}
	else
	{
		SEQAN_TASSERT(finderPos[0] == 5);
		SEQAN_TASSERT(keywordIndex[0] == 0);
		SEQAN_TASSERT(finderPos[1] == 6);
		SEQAN_TASSERT(keywordIndex[1] == 0);
		SEQAN_TASSERT(finderPos[2] == 5);
		SEQAN_TASSERT(keywordIndex[2] == 1);
	}
}



//////////////////////////////////////////////////////////////////////////////

template <typename TAlgorithmSpec>
void Test_OnlineAlgWildcards()
{
	String<unsigned int> pos;

//____________________________________________________________________________
// Test1 - simple find wo wildcards
	String<char> haystack("Dies ist ein Haystack. Ja, das ist wirklich einer!");
	Finder<String<char> > finder(haystack);

	String<char> needle("ist");
	Pattern<String<char>, TAlgorithmSpec> pattern(needle);
	clear(pos);
	while (find(finder, pattern))
		append(pos,position(finder));

	SEQAN_TASSERT(pos[0] == 7);
	SEQAN_TASSERT(pos[1] == 33);
	SEQAN_TASSERT(host(pattern) == needle);
	SEQAN_TASSERT(host(reinterpret_cast<Pattern<String<char>, TAlgorithmSpec> const &>(pattern)) == needle);

//____________________________________________________________________________
// Test - validation of patterns
	needle = "ist";
	setHost(pattern, needle);
	SEQAN_TASSERT(valid(pattern) == true);
	SEQAN_TASSERT(valid(reinterpret_cast<Pattern<String<char>, TAlgorithmSpec> const &>(pattern)) == true);

	needle = "i[a-z]s{3,4}t?a*a+c..\\a";
	setHost(pattern, needle);
	SEQAN_TASSERT(valid(pattern) == true);

	needle = "i[st";
	setHost(pattern, needle);
	SEQAN_TASSERT(valid(pattern) == false);

	needle = "ist\\";
	setHost(pattern, needle);
	SEQAN_TASSERT(valid(pattern) == false);

	needle = "ist?*";
	setHost(pattern, needle);
	SEQAN_TASSERT(valid(pattern) == false);

	needle = "ist{5,4}";
	setHost(pattern, needle);
	SEQAN_TASSERT(valid(pattern) == false);

	needle = "ist{5,a}";
	setHost(pattern, needle);
	SEQAN_TASSERT(valid(pattern) == false);

	needle = "ist{5,";
	setHost(pattern, needle);
	SEQAN_TASSERT(valid(pattern) == false);

//____________________________________________________________________________
// Test - searching with invalid needles
	haystack = "Dies i[st ein Haystack. Ja, das i[st wirklich einer!";
	setHost(finder, haystack);
	clear(finder);

	needle = "i[st";
	setHost(pattern, needle);

	clear(pos);
	while (find(finder, pattern))
		append(pos,position(finder));
	
	SEQAN_TASSERT(length(pos) == 0);

//____________________________________________________________________________
// Test - handle needles with wildcards
	// to produce two \ in the pattern you need to escape both of them
	needle = "aa+c*[a-z]xx?aa\\\\";
	SEQAN_TASSERT(_length_wo_wild(needle) == 9);
	
//____________________________________________________________________________
// Test - optional characters (?)
	haystack = "abc__ac";
	setHost(finder, haystack);
	clear(finder);

	needle = "ab?c";
	setHost(pattern, needle);

	clear(pos);
	while (find(finder, pattern))
		append(pos,position(finder));
	
	SEQAN_TASSERT(pos[0] == 2);
	SEQAN_TASSERT(pos[1] == 6);
	SEQAN_TASSERT(length(pos) == 2);

//____________________________________________________________________________
// Test - repeatable characters (+)
	haystack = "abc__abbbbbc_ac";
	setHost(finder, haystack);
	clear(finder);

	needle = "ab+c";
	setHost(pattern, needle);

	clear(pos);
	while (find(finder, pattern))
		append(pos,position(finder));
	
	SEQAN_TASSERT(pos[0] == 2);
	SEQAN_TASSERT(pos[1] == 11);
	SEQAN_TASSERT(length(pos) == 2);

//____________________________________________________________________________
// Test - repeatable characters (*)
	haystack = "abc__abbbbbc_ac";
	setHost(finder, haystack);
	clear(finder);

	needle = "ab*c";
	setHost(pattern, needle);

	clear(pos);
	while (find(finder, pattern))
		append(pos,position(finder));
	
	SEQAN_TASSERT(pos[0] == 2);
	SEQAN_TASSERT(pos[1] == 11);
	SEQAN_TASSERT(pos[2] == 14);
	
	SEQAN_TASSERT(length(pos) == 3);

//____________________________________________________________________________
// Test - wildcard matching
	haystack = "acccdfabdeeef";
	setHost(finder, haystack);
	clear(finder);

	needle = "ab?c*de+f";
	setHost(pattern, needle);

	clear(pos);
	while (find(finder, pattern))
		append(pos,position(finder));
	
	SEQAN_TASSERT(pos[0] == 12);
	SEQAN_TASSERT(length(pos) == 1);

//____________________________________________________________________________
// Test - wildcard matching (hard case)
	haystack = "aacccdfacccdeeef";
	setHost(finder, haystack);
	clear(finder);

	needle = "a*c*de+f";
	setHost(pattern, needle);

	clear(pos);
	while (find(finder, pattern))
		append(pos,position(finder));
	
	SEQAN_TASSERT(pos[0] == 15);
	SEQAN_TASSERT(length(pos) == 1);


//____________________________________________________________________________
// Test - character classes matching
	haystack = "annual_Annual_znnual_Znnual";
	setHost(finder, haystack);
	clear(finder);

	needle = "[a-zA]nnual";	
	setHost(pattern, needle);

	clear(pos);
	while (find(finder, pattern))
		append(pos,position(finder));
	
	SEQAN_TASSERT(pos[0] == 5);
	SEQAN_TASSERT(pos[1] == 12);
	SEQAN_TASSERT(pos[2] == 19);
	SEQAN_TASSERT(length(pos) == 3);

//____________________________________________________________________________
// Test - long needles
	haystack = "abcdefghijklmnopqrstuvwxyzabcdefghijklmnopqrstuvwxyzabcdefgaabcdef";
	setHost(finder, haystack);
	clear(finder);

	needle = "abcdefghijklmnopqrstuvwxyzabcdefg";	
	setHost(pattern, needle);

	clear(pos);
	while (find(finder, pattern))
		append(pos,position(finder));

	SEQAN_TASSERT(pos[0] == 32);
	SEQAN_TASSERT(pos[1] == 58);
	SEQAN_TASSERT(length(pos) == 2);

//____________________________________________________________________________
// Test - long needles with character classes
	//	        abcdefghijklmnopqrstuvwxyzabcdefghijkl
	//	                                  abcdefghijklmnopqrstuvwxyzabcdefghijkl
	haystack = "abcdefghijklmnopqrstuvwxyzabcdefghijklmnopqrstuzwxyzabcdefghzjklaabcdefhijkl";
	setHost(finder, haystack);
	clear(finder);

	needle = "abcdefghijklmnopqrstu[vz]wxyzabcdefgh[iz]jkl";	
	setHost(pattern, needle);

	clear(pos);
	while (find(finder, pattern))
		append(pos,position(finder));

	SEQAN_TASSERT(pos[0] == 37);
	SEQAN_TASSERT(pos[1] == 63);
	SEQAN_TASSERT(length(pos) == 2);


//____________________________________________________________________________
// Test - long needles with repeating characters
	//	        abcdefghijklmnopqrstuvwxyzabcdefghijkl
	//													  abcdefghijklmnopqrstuvwxyzabcdefghijkl
	haystack = "abcdefghijklmnopqrstuvwxyzabcdefghiiiiijkl____aaaaabcdefghijklmnopqrstuvwxyzabcdeghijkl__aaaaabcdefghijklmnopqrstuvwxyzabcdefghjkl";
	setHost(finder, haystack);
	clear(finder);

	needle = "aa*bcdefghijklmnopqrstuvwxyzabcdef?g?hi+jkl";	
	setHost(pattern, needle);

	clear(pos);
	while (find(finder, pattern))
		append(pos,position(finder));

	SEQAN_TASSERT(pos[0] == 41);
	SEQAN_TASSERT(pos[1] == 86);
	SEQAN_TASSERT(length(pos) == 2);

//____________________________________________________________________________
// Test - handle .
	haystack = "annual_Annual_znnual";
	setHost(finder, haystack);
	clear(finder);

	needle = ".nnual";	
	setHost(pattern, needle);
	clear(pos);
	while (find(finder, pattern)){
		append(pos,position(finder));
	}
	SEQAN_TASSERT(pos[0] == 5);
	SEQAN_TASSERT(pos[1] == 12);
	SEQAN_TASSERT(pos[2] == 19);
	SEQAN_TASSERT(length(pos) == 3);

//____________________________________________________________________________
// Test - handle backslash 
	haystack = "annual_Annual_.nnual";
	setHost(finder, haystack);
	clear(finder);

	needle = "\\.nnual";	
	setHost(pattern, needle);
	clear(pos);
	while (find(finder, pattern)){
		append(pos,position(finder));
	}
	SEQAN_TASSERT(pos[0] == 19);
	SEQAN_TASSERT(length(pos) == 1);



//____________________________________________________________________________
// Test - handle bounded length repeats {n,m}
	haystack = "aannual_aaaannual_annual";
	setHost(finder, haystack);
	clear(finder);

	needle = "a{2,5}n{2}ual";	
	
	SEQAN_TASSERT(_length_wo_wild(needle) == 10);
	
	setHost(pattern, needle);
	clear(pos);
	while (find(finder, pattern)){
		append(pos,position(finder));
	}
	SEQAN_TASSERT(pos[0] == 6);
	SEQAN_TASSERT(pos[1] == 16);
	SEQAN_TASSERT(length(pos) == 2);


//____________________________________________________________________________
// Test - handle different types of Pattern and Needle
	String<Dna> dna_haystack("AAACCTATGGGTTTAAAACCCTGAAACCCC");
	Finder<String<Dna> > dna_finder(dna_haystack);

	String<char> char_needle("a{3}c+t[ag].");
	Pattern<String<Dna>, TAlgorithmSpec> dna_pattern(char_needle);
	clear(pos);

	while (find(dna_finder, dna_pattern))
		append(pos,position(dna_finder));

	SEQAN_TASSERT(pos[0] == 7);
	SEQAN_TASSERT(pos[1] == 23);
	SEQAN_TASSERT(length(pos) == 2);

}


//////////////////////////////////////////////////////////////////////////////

template <typename TPatternSpec>
void Test_Approx_EditDist()
{
//test DPSearch
  	String<char> hstk("any_annealing");
	String<char> nl("annual");

	Finder<String<char> > fd(hstk);

	Pattern<String<char>, TPatternSpec> pt(nl, -2);
	SEQAN_TASSERT(find(fd, pt))
	SEQAN_TASSERT(position(fd) == 8)
	SEQAN_TASSERT(getScore(pt) == -2)
	SEQAN_TASSERT(find(fd, pt))
	SEQAN_TASSERT(position(fd) == 9)
	SEQAN_TASSERT(getScore(pt) == -1)
	SEQAN_TASSERT(find(fd, pt))
	SEQAN_TASSERT(position(fd) == 10)
	SEQAN_TASSERT(getScore(pt) == -2)

	SEQAN_TASSERT(!find(fd,pt))

	String<char> haystk("Dies ist der Haystack des Tests. Ja, das ist er wirklich!");
	String<char> ndl("des");

	Finder<String<char> > fnd(haystk);

	Pattern<String<char>, TPatternSpec> pat(ndl, -2);
	SEQAN_TASSERT(host(pat) == ndl);
	SEQAN_TASSERT(host(reinterpret_cast<Pattern<String<char>, TPatternSpec> const &>(pat)) == ndl);

	SEQAN_TASSERT(scoreLimit(pat) == -2)
	setScoreLimit(pat, -1);
	SEQAN_TASSERT(scoreLimit(pat) == -1)

	SEQAN_TASSERT(find(fnd, pat))
	SEQAN_TASSERT(position(fnd) == 3)
	SEQAN_TASSERT(getScore(pat) == -1)

	SEQAN_TASSERT(find(fnd, pat))
	SEQAN_TASSERT(position(fnd) == 10)
	SEQAN_TASSERT(getScore(pat) == -1)

	SEQAN_TASSERT(find(fnd, pat))
	SEQAN_TASSERT(position(fnd) == 11)
	SEQAN_TASSERT(getScore(pat) == -1)

	SEQAN_TASSERT(find(fnd, pat))
	SEQAN_TASSERT(position(fnd) == 23)
	SEQAN_TASSERT(getScore(pat) == -1)

	SEQAN_TASSERT(find(fnd, pat))
	SEQAN_TASSERT(position(fnd) == 24)
	SEQAN_TASSERT(getScore(pat) == 0)

	SEQAN_TASSERT(find(fnd, pat))
	SEQAN_TASSERT(position(fnd) == 25)
	SEQAN_TASSERT(getScore(pat) == -1)

	SEQAN_TASSERT(find(fnd, pat))
	SEQAN_TASSERT(position(fnd) == 28)
	SEQAN_TASSERT(getScore(pat) == -1)

	SEQAN_TASSERT(find(fnd, pat))
	SEQAN_TASSERT(position(fnd) == 39)
	SEQAN_TASSERT(getScore(pat) == -1)

	SEQAN_TASSERT(!find(fnd, pat))

	// Test with long needles and a Dna Alphabet
  	String<Dna> long_haystk("taaaataaaatacaataaaataaaataaaataaaataaaataaaataaaataaaataaaataaaataaaat");
	String<Dna> long_ndl("taaaataaaatacaataaaataaaatataataaaataaaataaaat");

	Finder<String<Dna> > long_fnd(long_haystk);

	Pattern<String<Dna>, TPatternSpec> long_pat(long_ndl, -2);

	SEQAN_TASSERT(find(long_fnd,long_pat))
	SEQAN_TASSERT(position(long_fnd) == 44)
	SEQAN_TASSERT(getScore(long_pat) == -2)

	SEQAN_TASSERT(find(long_fnd,long_pat))
	SEQAN_TASSERT(position(long_fnd) == 45)
	SEQAN_TASSERT(getScore(long_pat) == -1)

	SEQAN_TASSERT(find(long_fnd,long_pat))
	SEQAN_TASSERT(position(long_fnd) == 46)
	SEQAN_TASSERT(getScore(long_pat) == -2)

	SEQAN_TASSERT(find(long_fnd,long_pat))
	SEQAN_TASSERT(position(long_fnd) == 60)
	SEQAN_TASSERT(getScore(long_pat) == -2)

	SEQAN_TASSERT(find(long_fnd,long_pat))
	SEQAN_TASSERT(position(long_fnd) == 65)
	SEQAN_TASSERT(getScore(long_pat) == -2)

	SEQAN_TASSERT(find(long_fnd,long_pat))
	SEQAN_TASSERT(position(long_fnd) == 70)
	SEQAN_TASSERT(getScore(long_pat) == -2)

	SEQAN_TASSERT(!find(long_fnd,long_pat))

//____________________________________________________________________________

	String<char> haystack_1 = "123XXXabaXXX45aba123";
	String<char> needle_1 = "XXXaba";
	Finder<String<char> > finder_1(haystack_1);
	Pattern<String<char>, TPatternSpec> pattern_1(needle_1, -2);

	SEQAN_TASSERT(find(finder_1, pattern_1))
	SEQAN_TASSERT(endPosition(finder_1) == 7)
	SEQAN_TASSERT(getScore(pattern_1) == -2)
	SEQAN_TASSERT(findBegin(finder_1, pattern_1))
	SEQAN_TASSERT(infix(finder_1) == "XXXa")
	SEQAN_TASSERT(getBeginScore(pattern_1) == -2)
	SEQAN_TASSERT(!findBegin(finder_1, pattern_1))

	SEQAN_TASSERT(find(finder_1, pattern_1))
	SEQAN_TASSERT(endPosition(finder_1) == 8)
	SEQAN_TASSERT(getScore(pattern_1) == -1)
	SEQAN_TASSERT(findBegin(finder_1, pattern_1))
	SEQAN_TASSERT(infix(finder_1) == "XXab")
	SEQAN_TASSERT(getBeginScore(pattern_1) == -2)
	SEQAN_TASSERT(findBegin(finder_1, pattern_1))
	SEQAN_TASSERT(infix(finder_1) == "XXXab")
	SEQAN_TASSERT(getBeginScore(pattern_1) == -1)
	SEQAN_TASSERT(findBegin(finder_1, pattern_1))
	SEQAN_TASSERT(infix(finder_1) == "3XXXab")
	SEQAN_TASSERT(getBeginScore(pattern_1) == -2)
	SEQAN_TASSERT(!findBegin(finder_1, pattern_1))

	SEQAN_TASSERT(find(finder_1, pattern_1))
	SEQAN_TASSERT(endPosition(finder_1) == 9)
	SEQAN_TASSERT(getScore(pattern_1) == 0)
	SEQAN_TASSERT(findBegin(finder_1, pattern_1))
	SEQAN_TASSERT(infix(finder_1) == "Xaba")
	SEQAN_TASSERT(getBeginScore(pattern_1) == -2)
	SEQAN_TASSERT(findBegin(finder_1, pattern_1))
	SEQAN_TASSERT(infix(finder_1) == "XXaba")
	SEQAN_TASSERT(getBeginScore(pattern_1) == -1)
	SEQAN_TASSERT(findBegin(finder_1, pattern_1))
	SEQAN_TASSERT(infix(finder_1) == "XXXaba")
	SEQAN_TASSERT(getBeginScore(pattern_1) == 0)
	SEQAN_TASSERT(findBegin(finder_1, pattern_1))
	SEQAN_TASSERT(infix(finder_1) == "3XXXaba")
	SEQAN_TASSERT(getBeginScore(pattern_1) == -1)
	SEQAN_TASSERT(findBegin(finder_1, pattern_1))
	SEQAN_TASSERT(infix(finder_1) == "23XXXaba")
	SEQAN_TASSERT(getBeginScore(pattern_1) == -2)
	SEQAN_TASSERT(!findBegin(finder_1, pattern_1))

	SEQAN_TASSERT(find(finder_1, pattern_1))
	SEQAN_TASSERT(endPosition(finder_1) == 10)
	SEQAN_TASSERT(getScore(pattern_1) == -1)
	SEQAN_TASSERT(findBegin(finder_1, pattern_1))
	SEQAN_TASSERT(infix(finder_1) == "XXabaX")
	SEQAN_TASSERT(getBeginScore(pattern_1) == -2)
	SEQAN_TASSERT(findBegin(finder_1, pattern_1))
	SEQAN_TASSERT(infix(finder_1) == "XXXabaX")
	SEQAN_TASSERT(getBeginScore(pattern_1) == -1)
	SEQAN_TASSERT(findBegin(finder_1, pattern_1))
	SEQAN_TASSERT(infix(finder_1) == "3XXXabaX")
	SEQAN_TASSERT(getBeginScore(pattern_1) == -2)
	SEQAN_TASSERT(!findBegin(finder_1, pattern_1))

	SEQAN_TASSERT(find(finder_1, pattern_1))
	SEQAN_TASSERT(endPosition(finder_1) == 11)
	SEQAN_TASSERT(getScore(pattern_1) == -2)
	SEQAN_TASSERT(findBegin(finder_1, pattern_1))
	SEQAN_TASSERT(infix(finder_1) == "XXXabaXX")
	SEQAN_TASSERT(getBeginScore(pattern_1) == -2)
	SEQAN_TASSERT(!findBegin(finder_1, pattern_1))

	SEQAN_TASSERT(find(finder_1, pattern_1))
	SEQAN_TASSERT(endPosition(finder_1) == 15)
	SEQAN_TASSERT(getScore(pattern_1) == -2)
	SEQAN_TASSERT(findBegin(finder_1, pattern_1))
	SEQAN_TASSERT(infix(finder_1) == "XXX45a")
	SEQAN_TASSERT(getBeginScore(pattern_1) == -2)
	SEQAN_TASSERT(!findBegin(finder_1, pattern_1))

	SEQAN_TASSERT(find(finder_1, pattern_1))
	SEQAN_TASSERT(endPosition(finder_1) == 17)
	SEQAN_TASSERT(getScore(pattern_1) == -2)
	SEQAN_TASSERT(findBegin(finder_1, pattern_1))
	SEQAN_TASSERT(infix(finder_1) == "X45aba")
	SEQAN_TASSERT(getBeginScore(pattern_1) == -2)
	SEQAN_TASSERT(findBegin(finder_1, pattern_1))
	SEQAN_TASSERT(infix(finder_1) == "XX45aba")
	SEQAN_TASSERT(getBeginScore(pattern_1) == -2)
	SEQAN_TASSERT(findBegin(finder_1, pattern_1))
	SEQAN_TASSERT(infix(finder_1) == "XXX45aba")
	SEQAN_TASSERT(getBeginScore(pattern_1) == -2)
	SEQAN_TASSERT(!findBegin(finder_1, pattern_1))

	SEQAN_TASSERT(!find(finder_1, pattern_1))
//____________________________________________________________________________

}

//////////////////////////////////////////////////////////////////////////////

void Test_Approx()
{
//test DPSearch
	Test_Approx_EditDist<DPSearch<SimpleScore> >();
	Pattern<String<char>, DPSearch<SimpleScore> > pat1	;

	SimpleScore sc;
	scoreGap(sc) = -10;
	setScoringScheme(pat1, sc);
	SEQAN_TASSERT(scoreGap(scoringScheme(pat1)) == -10);

	scoreGap(sc) = -1;
	SEQAN_TASSERT(scoreGap(scoringScheme(pat1)) == -10);
	setScoringScheme(pat1, sc);
	SEQAN_TASSERT(scoreGap(scoringScheme(pat1)) == -1);


//test other edit distance algorithm
	Test_Approx_EditDist<Myers<> >();
	Test_Approx_EditDist<AbndmAlgo >();
	Test_Approx_EditDist<PexNonHierarchical>();
	Test_Approx_EditDist<PexHierarchical>();
	// test with different multifinder
	Test_Approx_EditDist< Pex<NonHierarchical,AhoCorasick > >();
	Test_Approx_EditDist< Pex<NonHierarchical,MultiBFAM<> > >();

//____________________________________________________________________________
/*
	Pattern<String<char> > patrn(ndl);
	setBeginPosition(patrn, 0);
	SEQAN_TASSERT(beginPosition(patrn) == 0);

	setEndPosition(patrn, length(ndl));
	SEQAN_TASSERT(endPosition(patrn) == length(ndl));

	SEQAN_TASSERT(segment(patrn) == ndl);
*/
}

//////////////////////////////////////////////////////////////////////////////
//test prefix search

template <typename TPatternSpec>
void Test_Approx_Prefix_EditDist()
{
	String<char> haystack_1 = "mississippi";
	String<char> needle_1 = "misssi";
	Finder<String<char> > finder_1(haystack_1);
	Pattern<String<char>, TPatternSpec> pattern_1(needle_1, -2);
	SEQAN_TASSERT(find(finder_1, pattern_1))
	SEQAN_TASSERT(position(finder_1) == 3)
	SEQAN_TASSERT(length(finder_1) == 4)
	SEQAN_TASSERT(beginPosition(finder_1) == 0)
	SEQAN_TASSERT(endPosition(finder_1) == 4)
	SEQAN_TASSERT(infix(finder_1) == "miss");
	SEQAN_TASSERT(find(finder_1, pattern_1))
	SEQAN_TASSERT(position(finder_1) == 4)
	SEQAN_TASSERT(infix(finder_1) == "missi");
	SEQAN_TASSERT(find(finder_1, pattern_1))
	SEQAN_TASSERT(position(finder_1) == 5)
	SEQAN_TASSERT(infix(finder_1) == "missis");
	SEQAN_TASSERT(find(finder_1, pattern_1))
	SEQAN_TASSERT(position(finder_1) == 6)
	SEQAN_TASSERT(find(finder_1, pattern_1))
	SEQAN_TASSERT(position(finder_1) == 7)
	SEQAN_TASSERT(!find(finder_1, pattern_1))


	String<char> haystack_2 = "yyyXXaba";
	String<char> needle_2 = "yyyaba";
	Finder<String<char> > finder_2(haystack_2);
	Pattern<String<char>, TPatternSpec> pattern_2(needle_2, -2);
	SEQAN_TASSERT(find(finder_2, pattern_2))
	SEQAN_TASSERT(position(finder_2) == 5)
	SEQAN_TASSERT(infix(finder_2) == "yyyXXa");
	SEQAN_TASSERT(find(finder_2, pattern_2))
	SEQAN_TASSERT(position(finder_2) == 7)
	SEQAN_TASSERT(infix(finder_2) == "yyyXXaba");
	SEQAN_TASSERT(!find(finder_2, pattern_2))


	String<char> haystack_3 = "testtexttext";
	String<char> needle_3 = "mismatch";
	Finder<String<char> > finder_3(haystack_3);
	Pattern<String<char>, TPatternSpec> pattern_3(needle_3, -2);
	SEQAN_TASSERT(!find(finder_3, pattern_3))


	String<char> haystack_4 = "testtext";
	String<char> needle_4 = "a longer mismatch";
	Finder<String<char> > finder_4(haystack_4);
	Pattern<String<char>, TPatternSpec> pattern_4(needle_4, -2);
	SEQAN_TASSERT(!find(finder_4, pattern_4))


	String<char> haystack_5 = "exactmatching";
	String<char> needle_5 = "exact";
	Finder<String<char> > finder_5(haystack_5);
	Pattern<String<char>, TPatternSpec> pattern_5(needle_5, 0);
	SEQAN_TASSERT(find(finder_5, pattern_5))
	SEQAN_TASSERT(position(finder_5) == 4)
	SEQAN_TASSERT(!find(finder_5, pattern_5))


	String<char> haystack_6 = "this is a text that is a bit longer than one machine word of 32 or 64 bits. AAXYX";
	String<char> needle_6 =   "this is a text that is a bit longer than one machine word of 32 or 64 bits. XYX";
	Finder<String<char> > finder_6(haystack_6);
	Pattern<String<char>, TPatternSpec> pattern_6(needle_6, -2);
	SEQAN_TASSERT(find(finder_6, pattern_6))
	SEQAN_TASSERT(infix(finder_6) == "this is a text that is a bit longer than one machine word of 32 or 64 bits. AAX");
	SEQAN_TASSERT(find(finder_6, pattern_6))
	SEQAN_TASSERT(infix(finder_6) == "this is a text that is a bit longer than one machine word of 32 or 64 bits. AAXYX");
	SEQAN_TASSERT(!find(finder_6, pattern_6))

}


//////////////////////////////////////////////////////////////////////////////

int main() 
{
	SEQAN_TREPORT("TEST BEGIN")

	Test_OnlineAlg<Simple>();	
	Test_OnlineAlg<Horspool>();	
	Test_OnlineAlg<ShiftAnd>();
	Test_OnlineAlg<ShiftOr>();
	Test_OnlineAlg<BndmAlgo>();
	Test_OnlineAlg<BFAM<Oracle> >();
	Test_OnlineAlg<BFAM<Trie> >();

	Test_OnlineAlgWildcards<WildShiftAnd>();

	Test_OnlineAlgMulti<AhoCorasick>(false);
	//Test_OnlineAlgMulti<MultipleShiftAnd>(false);  //leakt 
	//Test_OnlineAlgMulti<SetHorspool>();		//kompiliert nicht
	Test_OnlineAlgMulti<WuManber>(true);
	Test_OnlineAlgMulti<MultiBFAM<Oracle> >(true);
	Test_OnlineAlgMulti<MultiBFAM<Trie> >(true);

	Test_Approx_Prefix_EditDist<DPSearch<Score<>, FindPrefix> >();
	Test_Approx_Prefix_EditDist<Myers<FindPrefix> >();

	Test_Approx();



	
//	debug::verifyCheckpoints("projects/library/seqan/find/find_myers_ukkonen.h");

	debug::verifyCheckpoints("projects/library/seqan/find/find_wild_shiftand.h");
	debug::verifyCheckpoints("projects/library/seqan/find/find_horspool.h");
	debug::verifyCheckpoints("projects/library/seqan/find/find_base.h");
	debug::verifyCheckpoints("projects/library/seqan/find/find_shiftand.h");
	debug::verifyCheckpoints("projects/library/seqan/find/find_shiftor.h");
	debug::verifyCheckpoints("projects/library/seqan/find/find_bndm.h");
	debug::verifyCheckpoints("projects/library/seqan/find/find_bom.h");
	//debug::verifyCheckpoints("projects/library/seqan/find/find_quasar.h");
	debug::verifyCheckpoints("projects/library/seqan/find/find_ahocorasick.h");
	debug::verifyCheckpoints("projects/library/seqan/find/find_multiple_shiftand.h");
	debug::verifyCheckpoints("projects/library/seqan/find/find_set_horspool.h");
	debug::verifyCheckpoints("projects/library/seqan/find/find_wumanber.h");
	debug::verifyCheckpoints("projects/library/seqan/find/find_abndm.h");
	debug::verifyCheckpoints("projects/library/seqan/find/find_pex.h");

	SEQAN_TREPORT("TEST END")
	return 0;
}

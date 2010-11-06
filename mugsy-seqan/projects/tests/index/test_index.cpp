#include <seqan/platform.h>

#include <iostream>
#include <fstream>

//#define SEQAN_DEBUG
//#define SEQAN_TEST

#include "test_index_creation.h"
#include <seqan/align.h>
#include <seqan/find.h>
#include <seqan/index.h>


using namespace std;
using namespace seqan;


//////////////////////////////////////////////////////////////////////////////



void testBuild()
{
		typedef String<char> TText;
		typedef StringSet< TText, Owner<ConcatDirect<> > > TMulti;

		String<char> gen1, gen2, gen3;
		cout << open(gen1, "corpus/NC_000117.txt");
		cout << open(gen2, "corpus/NC_002620.txt");
		cout << open(gen3, "corpus/NC_007429.txt");

        Index<TMulti> esa;
		appendValue(indexText(esa), gen1);
		appendValue(indexText(esa), gen2);
		appendValue(indexText(esa), gen3);


		indexRequire(esa, ESA_SA());
		indexRequire(esa, ESA_LCP());
		indexRequire(esa, ESA_BWT());
		indexRequire(esa, ESA_ChildTab());

        save(esa, "corpus/chlamydia");
}

template <typename TIter>
inline void _printNode(TIter const it) 
{
		cout << countOccurrences(it) << "\t";
		cout << representative(it) << "\t";
//		cout << "parentEdgeLabel:" << parentEdgeLabel(it);
		cout << endl;
}

void testMultiIndex()
{
		typedef String<Dna5> TText;
		typedef StringSet< TText, Owner<> > TMulti;

		String<Dna5> t[6];
		//t[0] = "caterpillar";
		//t[1] = "catwoman";
		//t[2] = "pillow";
		//t[3] = "willow";
		//t[4] = "ill";
		//t[5] = "wow";

		t[0] = "caggctcgcgt";
		t[1] = "caggaacg";
		t[2] = "tcgttg";
		t[3] = "tggtcg";
		t[4] = "agg";
		t[5] = "ctg";

		t[0] = "ac";
		t[1] = "ac";
//		t[2] = "aatt";

        Index<TMulti> esa;
		for(unsigned i=0; i<2; ++i)
			appendValue(indexText(esa), t[i]);

		// efficient dfs iterator (hiding edges with empty labels)
		{
			cout << "BottomUp without empty edges" << endl;
			Iter<Index<TMulti>, VSTree< BottomUp<> > > it(esa);
			while (!atEnd(it)) {
				_printNode(it);
				goNext(it);
			}
		}

		// efficient dfs iterator
		{
			cout << endl << "BottomUp with empty edges" << endl;
			Iter<Index<TMulti>, VSTree< BottomUp<PostorderEmptyEdges> > > it(esa);
			while (!atEnd(it)) {
				_printNode(it);
				goNext(it);
			}
		}

		// topdown dfs iterator (hiding edges with empty labels)
		{
			cout << endl << "TopDown postorder without empty edges" << endl;
			Iter<Index<TMulti>, VSTree< TopDown<ParentLinks<Postorder> > > > it(esa);
			while (goDown(it)) ;
			while (!atEnd(it)) {
				_printNode(it);
				goNext(it);
			}
		}

		// topdown dfs iterator
		{
			cout << endl << "TopDown postorder with empty edges" << endl;
			Iter<Index<TMulti>, VSTree< TopDown<ParentLinks<PostorderEmptyEdges> > > > it(esa);
			while (goDown(it)) ;
			while (!atEnd(it)) {
				_printNode(it);
				goNext(it);
			}
		}

		// topdown dfs iterator (hiding edges with empty labels)
		{
			cout << endl << "TopDown preorder without empty edges" << endl;
			Iter<Index<TMulti>, VSTree< TopDown<ParentLinks<Preorder> > > > it(esa);
			goDown(it,'c');
			while (!atEnd(it)) {
				_printNode(it);
				goNext(it);
			}
		}

		// topdown dfs iterator
		{
			cout << endl << "TopDown preorder with empty edges" << endl;
			Iter<Index<TMulti>, VSTree< TopDown<ParentLinks<PreorderEmptyEdges> > > > it(esa);
			while (!atEnd(it)) {
				_printNode(it);
				goNext(it);
			}
		}

		// topdown iterator w/o parent links (hiding edges with empty labels)
		{
			cout << endl << "TopDown with empty edges" << endl;
			Iter<Index<TMulti>, VSTree< TopDown<HideEmptyEdges> > > it(esa);
			_printNode(it);
			while (goDown(it))
				_printNode(it);
		}

		// topdown iterator w/o parent links
		{
			cout << endl << "TopDown with empty edges" << endl;
			Iter<Index<TMulti>, VSTree< TopDown<EmptyEdges> > > it(esa);
			_printNode(it);
			while (goDown(it))
				_printNode(it);
		}
/*
		indexRequire(esa, ESA_SA());
		indexRequire(esa, ESA_BWT());
		for(int i=0; i<length(indexRawSA(esa)); ++i)
			cout << saAt(i,esa) << " " << bwtAt(i,esa) << "    " << suffix(t[getValueI1(saAt(i,esa))], getValueI2(saAt(i,esa))) << endl;

//		resize(indexLCP(esa), length(indexRawText(esa)));
//		createLCPTableExt(indexLCP(esa), indexText(esa), indexSA(esa), Kasai());
		indexRequire(esa, ESA_LCP());
		for(int i=0; i<length(indexRawSA(esa)); ++i)
			cout << lcpAt(i,esa) << "    " << suffix(t[getValueI1(saAt(i,esa))], getValueI2(saAt(i,esa))) << endl;

		for(int i=0; i<length(indexRawSA(esa)); ++i)
			cout << saAt(i,esa) << " = " << indexRawSA(esa)[i] << "    " << endl;
		for(int i=0; i<length(indexRawSA(esa)); ++i)
			cout << bwtAt(i,esa) << " = " << indexBWT(esa).tab[i] << "    " << endl;
		for(int i=0; i<length(indexRawSA(esa)); ++i)
			cout << lcpAt(i,esa) << " = " << indexLCP(esa)[i] << "    " << endl;
*/

/*

		resize(sa, length(indexRawText(esa)));
		createSuffixArrayExt(sa, indexText(esa), Skew7());

		for(int i=0; i<length(indexRawText(esa)); ++i)
			cout << indexRawText(esa)[i] << "    ";

		String<unsigned> lcp;
		resize(lcp, length(indexRawText(esa)));
		createLCPTableExt(lcp, indexText(esa), sa, Kasai());
*/
}

template <typename TIndex1, typename TIndex2>
bool compareTreeIterators(TIndex1 &index1, TIndex2 &index2)
{
	Iter<TIndex1, VSTree< TopDown< ParentLinks<Preorder> > > > it1(index1);
	Iter<TIndex2, VSTree< TopDown< ParentLinks<Preorder> > > > it2(index2);

	while (!atEnd(it1) && !atEnd(it2)) {
		if (representative(it1) != representative(it2)) 
		{
			cout << representative(it1) << " != " << representative(it2) << endl;
			return false;
		}
		goNext(it1);
		goNext(it2);
	}

	if (!(atEnd(it1) && atEnd(it1)))
		return false;

	return true;
}

template <typename TIndexSpec1, typename TIndexSpec2>
bool compareIndices()
{
	bool equal = true;
	{
		CharString text("mississippi");
		Index<CharString, TIndexSpec1> index1(text);
		Index<CharString, TIndexSpec2> index2(text);
		if (!compareTreeIterators(index1, index2)) {
			cout << text << " indices differ" << endl;
			equal = false;
		}
	}
	{
		DnaString text("acaaacatat");
		Index<DnaString, TIndexSpec1> index1(text);
		Index<DnaString, TIndexSpec2> index2(text);
		if (!compareTreeIterators(index1, index2)) {
			cout << text << " indices differ" << endl;
			equal = false;
		}
	}
	{
		StringSet<CharString> t;
		resize(t, 6);
		t[0] = "caterpillar";
		t[1] = "catwoman";
		t[2] = "pillow";
		t[3] = "willow";
		t[4] = "ill";
		t[5] = "wow";
		Index<StringSet<CharString>, TIndexSpec1> index1(t);
		Index<StringSet<CharString>, TIndexSpec2> index2(t);
		if (!compareTreeIterators(index1, index2)) {
			cout << " indices differ" << endl;
			equal = false;
		}
	}
	{
		StringSet<DnaString> t;
		resize(t, 6);
		t[0] = "caggctcgcgt";
		t[1] = "caggaacg";
		t[2] = "tcgttg";
		t[3] = "tggtcg";
		t[4] = "agg";
		t[5] = "ctg";
		Index<StringSet<DnaString>, TIndexSpec1> index1(t);
		Index<StringSet<DnaString>, TIndexSpec2 > index2(t);
		if (!compareTreeIterators(index1, index2)) {
			cout << " indices differ" << endl;
			equal = false;
		}
	}
	return equal;
}

template <typename TIndexSpec>
void testSTreeIterators()
{
		typedef Index<String<char>, TIndexSpec> TIndex;

		String<char> text("acaaacatatz");
//		String<char> text("AAAAAGGGGG");
		TIndex index(text);
		Iter<TIndex, VSTree< TopDown< ParentLinks<Preorder> > > > it(index);
		Iter<TIndex, VSTree< TopDown<> > > itNoLinks(it);	// test conversion
		//Iter<TIndex, VSTree< BottomUp<> > > it(index);

//		while (goDown(it));
		while (!atEnd(it)) {
//			cout << countOccurrences(it) << "\t";
			cout << representative(it) << "\t";
			cout << "parentEdgeLabel: " << parentEdgeLabel(it); // << " " << value(it).node << "  " << value(it).range;
			cout << endl;
			goNext(it);
		}
			cout << endl;
		_dump(index);
/*		goBegin(it);
		while (!atEnd(it)) {
			cout << countOccurrences(it) << "\t";
			cout << representative(it) << "\t";
			cout << "parentEdgeLabel: " << parentEdgeLabel(it) << " " << value(it).node << "  " << value(it).range;
			cout << endl;
			goNext(it);
		}
		_dump(index);
*/}


void testMaxRepeats()
{
//		typedef String<char, External<> > TText;
		typedef String<char> TText;

        Index<TText> esa;
//        open(esa, "corpus/NC_000117.txt");
		indexText(esa) = "HALLOBALLOHALLEBALLO";

		FILE* dotFile = fopen("stree.dot","w");
		write(dotFile, esa, DotDrawing());
		fclose(dotFile);

        Iterator< Index<TText>, MaxRepeats >::Type it(esa, 3);
		typedef MaxRepeat< Index<TText> > TRepeat;
        while (!atEnd(it)) {
			cout << representative(it) << ":";
			Iterator<TRepeat, MaxRepeatOccurrences>::Type mit(it);
			while (!atEnd(mit)) {
				cout << "\t" << *mit;
				++mit;
			}
			cout << endl;
            ++it;
        }
}


void testMultiMEMs()
{
		typedef String<char> TText;
		typedef StringSet< TText, Owner<ConcatDirect<> > > TMulti;

        Index<TMulti> esa;

		String<Dna5> t[6];
		t[0] = "caterpillar";
		t[1] = "catwoman";
		t[2] = "pillow";
		t[3] = "willow";
		t[4] = "ill";
		t[5] = "wow";

		FILE* dotFile = fopen("stree.dot","w");
		write(dotFile, esa, DotDrawing());
		fclose(dotFile);

        Iterator< Index<TMulti>, MultiMEMs >::Type it(esa, 3);
		typedef MultiMEM< Index<TMulti> > TMultiMEM;
        while (!atEnd(it)) {
			cout << representative(it) << ":";
			Iterator<TMultiMEM>::Type mit(it);
			while (!atEnd(mit)) {
				cout << "\t" << *mit;
				++mit;
			}
			cout << endl;
            ++it;
        }
}

template <typename TIteratorSpec>
void testSuperMaxRepeats()
{
//		typedef String<char, External<> > TText;
		typedef String<char> TText;

        Index<TText> esa;
//        open(esa, "corpus/NC_000117.txt");
		indexText(esa) = "HALLOBALLOHALLEBALLO";

        typename Iterator< Index<TText>, TIteratorSpec >::Type it(esa);
        while (!atEnd(it)) {
			cout << representative(it) << ":";
			for(typename Size<Index<TText> >::Type i = 0; i < countOccurrences(it); ++i)
				cout << "\t" << getOccurrences(it)[i];
			cout << endl;
            ++it;
        }
}


void testMUMs()
{
		typedef String<char> TText;
		typedef StringSet< TText, Owner<ConcatDirect<> > > TMulti;
		typedef Index<TMulti, Index_ESA<> > TIndex;

		String<char> t[3];

		t[0] = "fefhalloballo";
		t[1] = "halloballefser";
		t[2] = "grballoballo";


        TIndex esa;
		for(int i = 0; i < 3; ++i)
			appendValue(indexText(esa), t[i]);			// add sequences to multiple index

		Iterator<TIndex, MUMs>::Type  it(esa, 3);		// set minimum MUM length to 3
		String< SAValue<TIndex>::Type > occs;			// temp. string storing the hit positions

		cout << resetiosflags(ios::left);
		while (!atEnd(it)) 
		{
			cout << representative(it) << ":";
			occs = getOccurrences(it);					// gives hit positions (seqNo,seqOfs)
			orderOccurrences(occs);						// order them by seqNo
			
			for(unsigned i = 0; i < length(occs); ++i)
				cout << "\t" << getValueI2(occs[i]);
			cout << endl;
			cout << alignment(it) << endl;

			++it;
		}
}


template <typename TAlgorithmSpec>
void testFind()
{
		String<unsigned int> pos;

	//____________________________________________________________________________
	// Test1 - small needle

		String<char> haystack("Dies ist ein Haystack. Ja, das ist wirklich einer!");
		Index<String<char> > index(haystack);

		Finder<Index< String<char> >, TAlgorithmSpec> finder(index);

		String<char> needle1("ist");
		seqan::Pattern<String<char> > pattern(needle1);	

		while (find(finder, pattern))
			append(pos,position(finder));

		SEQAN_TASSERT(pos[0] == 5);
		SEQAN_TASSERT(pos[1] == 31);

	//____________________________________________________________________________
	// Test2 - large needle

		haystack = "abcdefghijklmnopqrstuvwxyzabcdefghijklmnopqrstuvwxyzabcdefgaabcdef";
		clear(index);
		clear(finder);

		needle1 = "abcdefghijklmnopqrstuvwxyzabcdefg";
		setNeedle(pattern, needle1);

		clear(pos);
		while (find(finder, pattern))
			append(pos,position(finder));

		SEQAN_TASSERT(pos[1] == 0);
		SEQAN_TASSERT(pos[0] == 26);
}

void dumpStatus(bool status) {
	if (status)
		cout << "OK" << endl;
	else
		cout << "FAILED!" << endl;
	}


void testShapes();
void testQGramIndex();
bool testIndexCreation();

int main()
{
	SEQAN_TREPORT("TEST BEGIN")

	cout << "===================================" << endl;
	cout << "----Basic Suffix Tree iterators----" << endl;
	testSTreeIterators<Index_Wotd<> >();
	testSTreeIterators<Index_Wotd<WotdOriginal> >();
	testSTreeIterators<Index_ESA<> >();
	cout << "===================================" << endl;
	cout << "----Compare Indices----------------" << endl;
	dumpStatus(compareIndices<Index_ESA<>, Index_Wotd<> >());
	cout << "===================================" << endl;
	cout << "----Suffix Array based Finder------" << endl;
	testFind<ESA_FIND_MLR>();
	cout << "===================================" << endl;
	cout << "----Multiple Sequence Index--------" << endl;
	testMultiIndex();
	cout << "===================================" << endl;
	cout << "----MUMs---------------------------" << endl;
	testMUMs();
	cout << "===================================" << endl;
	cout << "----Maximal Repeats----------------" << endl;
	testMaxRepeats();
	cout << "===================================" << endl;
	cout << "----Supermaximal Repeats-----------" << endl;
	testSuperMaxRepeats<SuperMaxRepeats>();
	cout << "===================================" << endl;
	cout << "----Supermaximal Repeats---fast----" << endl;
	testSuperMaxRepeats<SuperMaxRepeatsFast>();
	cout << "===================================" << endl;
	cout << endl;

	cout << "===================================" << endl;
	cout << "----QGram Index--------------------" << endl;
	testQGramIndex();
	cout << "===================================" << endl;
	cout << "----Shapes-------------------------" << endl;
	testShapes();
	cout << "===================================" << endl;
	cout << endl;

	cout << "===================================" << endl;
	cout << "----SA, LCP, and ChildTab test-----" << endl;
	testIndexCreation();
	cout << "===================================" << endl;
//	testBuild();

	SEQAN_TREPORT("TEST END")
	return 0;
}
